import os
import glob
import subprocess
import GeneralUtilityMethods
from OptimizeProtein import yasara

# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

yasara.info.mode = 'txt'


# The original OptProt script was designed as a type of main script that is responsible for determining which
# calculations will be made on which pdbs, based on reading in and parsing the options file.
# It includes preparatory steps for calculating protein structural and aggregation propensity properties,
# using Yasara functions. It also includes calls to FoldX to perform repair of the structure data.
# This version of OptProt does exactly the same things as the original. It is simply refactored to be object-oriented.
class OptProt(object):
    # declare all global variables
    __r_path__ = ''
    __foldx_path__ = ''
    __agadir_path__ = ''
    __qsub_path__ = ''
    __start_path__ = ''
    __scripts_path__ = ''
    __execute_python_and_path_to_script_fwdslash__ = ''
    __execute_r_and_path_to_script_fwdslash__ = ''
    __proteinChains_and_paths_to_r_foldx_agadir__ = ''
    __proteinChains_and_paths_to_r_foldx_agadir_and_charge__ = ''
    __repair_job_prefix__ = ''

    def __init__(self, start_path, scripts_path, r_path, foldx_path, agadir_path, qsub_path):
        global __start_path__
        global __scripts_path__
        global __r_path__
        global __foldx_path__
        global __agadir_path__
        global __qsub_path__
        global __execute_python_and_path_to_script_fwdslash__
        global __execute_r_and_path_to_script_fwdslash__
        global __repair_job_prefix__
        __start_path__ = start_path
        __scripts_path__ = scripts_path
        __r_path__ = r_path
        __foldx_path__ = foldx_path
        __agadir_path__ = agadir_path
        __qsub_path__ = qsub_path
        __repair_job_prefix__ = 'RPjob_'
        self.__single_space__ = ' '
        __execute_python_and_path_to_script_fwdslash__ = 'python' + self.__single_space__ + __scripts_path__ + '/'
        __execute_r_and_path_to_script_fwdslash__ = 'R <' + self.__single_space__ + __scripts_path__ + '/'
        self.__pdb_list__ = []
        self.__command__ = ''
        self.__charge__ = ''
        self.__proteinChains__ = ''
        self._print_absolute_path_to('R', __r_path__)
        self._print_absolute_path_to('FoldX', __foldx_path__)
        self._print_absolute_path_to('TANGO', __agadir_path__)
        self._print_absolute_path_to('Qsub', __qsub_path__)

    # Extracts instructions from the option file.
    # This includes which computations to run & which proteins to run them on in 2 parts:
    # 1. list of pdb names
    # 2. computation names ("command"), charge, protein chains
    def parse_option_file(self, __option_file__):
        self.__parse_optionfile_for_pdblist(__option_file__)
        self.__parse_optionfile_for_computations_charge_proteinChains(__option_file__)
        global __proteinChains_and_paths_to_r_foldx_agadir__
        global __proteinChains_and_paths_to_r_foldx_agadir_and_charge__
        __proteinChains_and_paths_to_r_foldx_agadir__ = self.__single_space__ + self.__proteinChains__ + \
            self.__single_space__ + __r_path__ + self.__single_space__ + __foldx_path__ + self.__single_space__ + \
            __agadir_path__
        __proteinChains_and_paths_to_r_foldx_agadir_and_charge__ = __proteinChains_and_paths_to_r_foldx_agadir__ + \
            self.__single_space__ + self.__charge__

    # Tidies up each pdb file using Yasara functions.
    # Converts each pdb to fasta protein sequence via pdb2fasta.py.
    # Runs Agadirwrapper on each fasta protein sequence via agadir.py.
    # Runs FoldX repair on each pdb via repair.py.
    def run_yasara_agadir_repair(self):
        global __repair_job_prefix__
        if not os.path.exists('Results'):
            os.makedirs('Results')
        for pdb in self.__pdb_list__:
            pdb_name = pdb.split('.')[0]
            self._build_results_directory_tree_for_each(pdb)  # also changes current directory into Results/pdb
            self._run_yasara_to_organise_pdb(pdb, pdb_name)
            self._copy_pdb_foldx_agadir_files_to_new_subdirectories(pdb)
            self._run_agadir()
            self._run_repair_on_grid_engine(pdb_name)
        GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__repair_job_prefix__, 'PDBs to be repaired')

    def perform_selected_computations(self):
        for pdb in self.__pdb_list__:
            self._build_results_directory_tree_for_each(pdb)  # also changes current directory into Results/pdb
            self._compute_stretchplot(pdb)
            if not os.path.exists('Runs/' + self.__command__):
                self._build_directory_tree_for_computations()
            self._compute_commands()
            os.chdir(__start_path__)

    # # # # # # # # # PRIVATE METHODS # # # # # # # # #

    # # # # # # Called by parse_option_file() # # # # #

    # 1. Assigns to object variables the list of pdb names specified in the option file.
    def __parse_optionfile_for_pdblist(self, __option_file__):
        self.__pdb_list__ = []
        for line in __option_file__:
            if '#' in line:
                continue
            if 'PDBs:' in line:
                if line.split(':')[-1].strip(';\n') == "All":
                    pdb_paths = glob.glob('./PDBs/*.pdb')
                    for pdb_path in pdb_paths:
                        pdb_temp = pdb_path.split('/')[-1]
                        self.__pdb_list__.append(pdb_temp)
                pdb_string = line.split(':')[-1].strip(';\n')
                if ',' in pdb_string:
                    pdbs_temp = pdb_string.split(',')
                    for pdb_temp in pdbs_temp:
                        self.__pdb_list__.append(pdb_temp)
                else:
                    self.__pdb_list__.append(pdb_string)
        print 'PDBs to analyse:\t\t' + ",\t".join(self.__pdb_list__)

    # 2. Assigns to object variables the computation names ("command"), charge, protein chains specified
    # in the option file.
    def __parse_optionfile_for_computations_charge_proteinChains(self, __option_file__):
        for line in __option_file__:
            if '#' in line:
                continue
            if 'Command:' in line:
                self.__command__ = line.split(':')[-1].strip(';\n')
            if 'Charge:' in line:
                self.__charge__ = line.split(':')[-1].strip(';\n')
            if 'ProteinChains:' in line:
                self.__proteinChains__ = []
                protein_chains_string = line.split(':')[-1].strip(';\n')
                if ',' in protein_chains_string:
                    protein_chains_temp = protein_chains_string.split(',')
                    for protein_chain_temp in protein_chains_temp:
                        self.__proteinChains__.append(protein_chain_temp)
                    self.__proteinChains__ = "_".join(self.__proteinChains__)
                else:
                    self.__proteinChains__ = protein_chains_string
        print 'Command to be executed:\t\t' + self.__command__
        print 'Protein chains to be considered:\t' + self.__proteinChains__

    # # # # Called by run_yasara_agadir_repair() # # # #

    # Yasara uses the term "molecules" to describe protein chains,
    # so the variable names in this method are mol instead of protein_chain.
    def _run_yasara_to_organise_pdb(self, pdb, pdb_name):
        yasara.run('DelObj all')
        yasara.run('LoadPDB' + self.__single_space__ + __start_path__ + '/PDBs/' + pdb)
        yasara.run('DelRes !Protein')
        temp_mols = yasara.run('ListMol All,Format=MOLNAME')
        for mol in temp_mols:
            yasara.run('RenumberRes all and Mol' + self.__single_space__ + mol + ',First=1')
        yasara.run(
            'SavePDB 1,' + __start_path__ + '/Results/' + pdb_name + '/PDBs/' + pdb + ',Format=PDB,Transform=Yes')

    def _copy_pdb_foldx_agadir_files_to_new_subdirectories(self, pdb):
        cp_start_path = 'cp' + self.__single_space__ + __start_path__
        subprocess.call(cp_start_path + '/PDBs/' + pdb + ' ./PDBs/.', shell=True)
        subprocess.call(cp_start_path + '/PDBs/' + pdb + ' ./Repair/.', shell=True)
        subprocess.call(cp_start_path + '/SourceFiles/FoldXFiles/* ./Repair/.', shell=True)
        subprocess.call(cp_start_path + '/SourceFiles/AgadirFiles/* ./Agadir/.', shell=True)

    def _run_repair_on_grid_engine(self, pdb_name):
        global __repair_job_prefix__
        global __foldx_path__
        global __scripts_path__
        repair_python_script = 'repair.py'
        self._print_OptProt_calling_script(repair_python_script)
        grid_engine_job_name = __repair_job_prefix__ + pdb_name
        no_queue = ''
        no_max_memory = ''
        no_cluster = ''
        python_script_with_path_and_qsub = __scripts_path__ + '/' + repair_python_script + \
                                self.__single_space__ + __qsub_path__
        GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                   __foldx_path__, python_script_with_path_and_qsub)
        subprocess.call(__qsub_path__ + 'qsub job.q', shell=True)
        os.chdir(__start_path__)

    def _run_agadir(self):
        pdb2fasta_python_script = 'pdb2fasta.py'
        agadir_python_script = 'agadir.py'
        self._print_OptProt_calling_script(pdb2fasta_python_script)
        subprocess.call(__execute_python_and_path_to_script_fwdslash__ + pdb2fasta_python_script, shell=True)
        self._print_OptProt_calling_script(agadir_python_script)
        subprocess.call(__execute_python_and_path_to_script_fwdslash__ + agadir_python_script, shell=True)

    # # # # Called by perform_selected_computations() # # # #

    def _compute_stretchplot(self, pdb):
        stretchplot_python_script = 'stretchplot.py'
        stretchplot_r_script = 'stretchplot.R'
        self._print_OptProt_calling_script(stretchplot_python_script)
        subprocess.call(
            __execute_python_and_path_to_script_fwdslash__ + stretchplot_python_script + self.__single_space__ +
            __proteinChains_and_paths_to_r_foldx_agadir__, shell=True)
        self._print_OptProt_calling_script(stretchplot_r_script)
        subprocess.call(
            __execute_r_and_path_to_script_fwdslash__ + stretchplot_r_script + self.__single_space__ + '--no-save',
            shell=True)

    def _build_directory_tree_for_computations(self):
        if self.__command__ == 'All':
            self._build_directory_tree_for_stabilize_solubis_solubismild_cysscan()
        elif self.__command__ != 'Stretchplot':
            self._create_Runs_subfolder_for(self.__command__)

    def _build_directory_tree_for_stabilize_solubis_solubismild_cysscan(self):
        self._create_Runs_subfolder_for('Stabilize', 'Solubis', 'SolubisMild', 'CysScan')

    def _create_Runs_subfolder_for(self, *args):
        for command in args:
            os.makedirs('Runs/' + command)

    def _compute_commands(self):
        if self.__command__ == 'All':
            self._compute_stabilize_solubis_solubismild_cysscan()
        else:
            self._compute(self.__command__)

    def _compute_stabilize_solubis_solubismild_cysscan(self):
        self._compute('Stabilize', 'Solubis', 'SolubisMild', 'CysScan')

    def _compute(self, *args):
        for command in args:
            python_script = self._convert_command_name_to_python_script_name(command)
            self._print_OptProt_calling_script(python_script)
            if self.__command__ == 'Supercharge' or self.__command__ == 'DelPos' or self.__command__ == 'Indiv':
                subprocess.call(
                    __execute_python_and_path_to_script_fwdslash__ + python_script +
                    __proteinChains_and_paths_to_r_foldx_agadir_and_charge__, shell=True)
            else:
                subprocess.call(
                    __execute_python_and_path_to_script_fwdslash__ + python_script +
                    __proteinChains_and_paths_to_r_foldx_agadir__, shell=True)

    def _convert_command_name_to_python_script_name(self, command):
        python_script_name = command
        if self.__command__ == 'Stabilize':
            python_script_name = 'stabilize'
        elif self.__command__ == 'Solubis':
            python_script_name = 'solubis'
        elif self.__command__ == 'SolubisMild':
            python_script_name = 'solubis_mild'
        elif self.__command__ == 'SolubisDouble':
            python_script_name = 'solubis_double'
        elif self.__command__ == 'Supercharge':
            python_script_name = 'supercharge'
        return python_script_name + '.py'

    # # # # # # # # # UTILITY METHODS # # # # # # # # #

    # Note that this method also changes the current directory to the new Results/pdb directory
    def _build_results_directory_tree_for_each(self, pdb):
        pdb_name = pdb.split('.')[0]
        if not os.path.exists('Results/' + pdb_name):
            os.makedirs('Results/' + pdb_name)
        os.chdir('Results/' + pdb_name)
        self._create_folder_for('Repair', 'Agadir', 'Runs', 'PDBs', 'Fasta')

    def _create_folder_for(self, *args):
        for folder_name in args:
            if not os.path.exists(folder_name):
                os.makedirs(folder_name)

    def _print_OptProt_calling_script(self, python_script):
        print 'OptProt.py calling ' + python_script + '.......'

    def _print_absolute_path_to(self, target, path_to_target):
        print 'Absolute path to ' + target + ':' + '\t\t' + path_to_target


# pydevd.stoptrace()
