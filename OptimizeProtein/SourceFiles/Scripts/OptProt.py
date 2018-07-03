from OptimizeProtein import yasara
import os
import glob
import subprocess
from GeneralUtilityMethods import GUM
from Agadir import Agadir
from Solubis import Solubis
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

yasara.info.mode = 'txt'


# The original OptProt script was designed as a type of main script that is responsible for determining which
# calculations will be made on which pdbs, based on reading in and parsing the options file.
# It includes preparatory steps for calculating protein structural and aggregation propensity properties,
# using Yasara functions. It also includes calls to FoldX to perform repair of the structure data.
# This version of OptProt does exactly the same things as the original. It is simply refactored to be object-oriented.
# The methods work fine when a single pdb is specified in the Options_file, but has not been tested where more than 1
# pdb is specified. It is likely that this latter requirement would be resolved with a few changes. To begin that
# process, I have added a new dictionary data structure pdb_protein_chains_dict.)
class OptProt(object):

    def __init__(self, start_path, scripts_path, r_path, foldx_path, agadir_path, qsub_path):
        self.start_path = start_path
        self.scripts_path = scripts_path
        self.r_path = r_path
        self.foldx_path = foldx_path
        self.agadir_path = agadir_path
        self.qsub_path = qsub_path
        self.repair_job_prefix = 'RPjob_'
        self.space = ' '
        self.execute_python_and_path_to_script_fwdslash = 'python' + self.space + self.scripts_path + '/'
        self.execute_r_and_path_to_script_fwdslash = 'R <' + self.space + self.scripts_path + '/'
        self.pdb_list = []
        self.command = ''
        self.charge = ''
        self.underscore_separated_proteinChains = ''
        self.pdb_protein_chains_dict = {}
        self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir = ''
        self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir_and_charge = ''
        self._print_absolute_path_to('R', self.r_path)
        self._print_absolute_path_to('FoldX', self.foldx_path)
        self._print_absolute_path_to('TANGO', self.agadir_path)
        self._print_absolute_path_to('Qsub', self.qsub_path)

    # Extracts instructions from the option file.
    # This includes which computations to run & which proteins to run them on in 2 parts:
    # 1. list of pdb names
    # 2. computation names ("command"), charge, protein chains
    def parse_option_file(self, option_file):
        self._parse_optionfile_for_pdblist(option_file)
        self._parse_optionfile_for_computations_charge_proteinChains(option_file)
        self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir = ''.join(self.underscore_separated_proteinChains) + self.space + self.r_path + \
                                                         self.space + self.foldx_path + \
                                                         self.space + self.agadir_path
        self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir_and_charge = self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir + \
                                                                    self.space + self.charge

    # Tidies up each pdb file using Yasara functions.
    # Converts each pdb to fasta protein sequence via pdb2fasta.py.
    # Runs Agadirwrapper on each fasta protein sequence via agadir.py.
    # Runs FoldX repair on each pdb via repair.py.
    def run_yasara_agadir_repair(self):
        if not os.path.exists('Results'):
            os.makedirs('Results')
        for pdb in self.pdb_list:
            pdb_name = pdb.split('.')[0]
            self._build_results_directory_tree_for_each(pdb)
            # current directory Results/pdb
            # self._run_yasara_to_organise_pdb(pdb, pdb_name)
            self._copy_pdb_foldx_agadir_files_to_new_subdirectories(pdb)
            relative_path_of_pdb_to_read = './PDBs/'
            relative_path_for_new_fasta_folder_to_write = './'
            pdb_name_chain_fasta_dict = GUM.extract_pdb_name_fasta_and_chains_from_pdb(pdb, relative_path_of_pdb_to_read)
            GUM.write_fasta_to_folder(pdb_name_chain_fasta_dict, relative_path_for_new_fasta_folder_to_write)
            agadir_results_path = self.start_path + '/Results/' + pdb_name
            agadir_instance = Agadir(agadir_results_path)
            agadir_instance.run_agadir_with_fasta_files('./')
            self._run_repair_on_grid_engine(pdb_name)
        message_to_print = 'PDBs to be repaired'
        GUM.wait_for_grid_engine_job_to_complete(self.repair_job_prefix, message_to_print)

    def perform_selected_computations(self):
        for pdb in self.pdb_list:
            self._build_results_directory_tree_for_each(pdb)  # also changes current directory into Results/pdb
            self._compute_stretchplot(pdb)
            if not os.path.exists('Runs/' + self.command):
                self._build_directory_tree_for_computations()
            self._compute_commands()
            os.chdir(self.start_path)

    # # # # # # # # # PRIVATE METHODS # # # # # # # # #

    # # # # # # Called by parse_option_file() # # # # #

    # 1. Assigns to object variables the list of pdb names specified in the option file.
    def _parse_optionfile_for_pdblist(self, option_file):
        for line in option_file:
            if '#' in line:
                continue
            if 'PDBs:' in line:
                if line.split(':')[-1].strip(';\n') == "All":
                    pdb_paths = glob.glob('./PDBs/*.pdb')
                    for pdb_path in pdb_paths:
                        pdb_temp = pdb_path.split('/')[-1]
                        self.pdb_list.append(pdb_temp)
                pdb_string = line.split(':')[-1].strip(';\n')
                if ',' in pdb_string:
                    pdbs_temp = pdb_string.split(',')
                    for pdb_temp in pdbs_temp:
                        self.pdb_list.append(pdb_temp)
                else:
                    self.pdb_list.append(pdb_string)
        print 'PDBs to analyse:\t\t' + ",\t".join(self.pdb_list)

    # TAKE NOTE: This method has a bug for assigning option_file-specified protein chains to a single global variable
    # (called underscore_separated_proteinChains). It would fail in the scenario where the option_file specifies more
    # than 1 pdb.
    def _parse_optionfile_for_computations_charge_proteinChains(self, option_file):
        for line in option_file:
            if '#' in line:
                continue
            if 'Command:' in line:
                self.command = line.split(':')[-1].strip(';\n').strip()
            if 'Charge:' in line:
                self.charge = line.split(':')[-1].strip(';\n').strip()
            if 'ProteinChains:' in line:
                self.underscore_separated_proteinChains = []
                protein_chains_option = line.split(':')[-1].strip(';\n')
                if protein_chains_option == 'All':
                    for pdb in self.pdb_list:
                        protein_chains_list = GUM.extract_all_chains_from_pdb(pdb, './PDBs/')
                        self.underscore_separated_proteinChains = '_'.join(protein_chains_list)
                        self.pdb_protein_chains_dict = {pdb: self.underscore_separated_proteinChains}
                elif ',' in protein_chains_option:
                    protein_chains_list = protein_chains_option.split(',')
                    protein_chains_list_stripped = []
                    for protein_chain in protein_chains_list:
                        protein_chains_list_stripped.append(protein_chain.strip())
                    self.underscore_separated_proteinChains = '_'.join(protein_chains_list_stripped)
                else:
                    self.underscore_separated_proteinChains = protein_chains_option
        print 'Command to be executed:\t\t' + self.command
        print 'Protein chains to be considered:\t' + self.underscore_separated_proteinChains

    # # # # Called by run_yasara_agadir_repair() # # # #

    # Yasara uses the term "molecules" to describe protein chains,
    # so the variable names in this method are "mol" instead of "protein_chain.
    def _run_yasara_to_organise_pdb(self, pdb, pdb_name):
        yasara.run('DelObj all')
        yasara.run('LoadPDB' + self.space + self.start_path + '/PDBs/' + pdb)
        yasara.run('DelRes !Protein')
        temp_mols = yasara.run('ListMol All,Format=MOLNAME')
        for mol in temp_mols:
            yasara.run('RenumberRes all and Mol' + self.space + mol + ',First=1')
        yasara.run(
            'SavePDB 1,' + self.start_path + '/Results/' + pdb_name + '/PDBs/' + pdb + ',Format=PDB,Transform=Yes')

    # The current working directory is Results/<pdb_name>0
    def _copy_pdb_foldx_agadir_files_to_new_subdirectories(self, pdb):
        cp_start_path = 'cp' + self.space + self.start_path
        subprocess.call(cp_start_path + '/PDBs/' + pdb + ' ./PDBs/.', shell=True)
        subprocess.call(cp_start_path + '/PDBs/' + pdb + ' ./Repair/.', shell=True)
        subprocess.call(cp_start_path + '/SourceFiles/FoldXFiles/* ./Repair/.', shell=True)
        subprocess.call(cp_start_path + '/SourceFiles/AgadirFiles/* ./Agadir/.', shell=True)

    def _run_repair_on_grid_engine(self, pdb_name):
        repair_python_script = 'repair.py'
        self._print_OptProt_calling_script(repair_python_script)
        grid_engine_job_name = self.repair_job_prefix + pdb_name
        no_queue = ''
        no_max_memory = ''
        no_cluster = ''
        using_runscript = False
        python_script_with_path_and_qsub = self.scripts_path + '/' + repair_python_script + \
                                           self.space + self.qsub_path
        GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                   using_runscript, self.foldx_path, python_script_with_path_and_qsub)
        subprocess.call(self.qsub_path + 'qsub job.q', shell=True)
        os.chdir(self.start_path)

    # # # # Called by perform_selected_computations() # # # #

    def _compute_stretchplot(self, pdb):
        stretchplot_python_script = 'stretchplot.py'
        stretchplot_r_script = 'stretchplot.R'
        self._print_OptProt_calling_script(stretchplot_python_script)
        subprocess.call(
            self.execute_python_and_path_to_script_fwdslash + stretchplot_python_script + self.space +
            self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir, shell=True)
        self._print_OptProt_calling_script(stretchplot_r_script)
        subprocess.call(
            self.execute_r_and_path_to_script_fwdslash + stretchplot_r_script + self.space + '--no-save',
            shell=True)

    def _build_directory_tree_for_computations(self):
        if self.command == 'All':
            self._build_directory_tree_for_stabilize_solubis_solubismild_cysscan()
        elif self.command != 'Stretchplot':
            self._create_Runs_subfolder_for(self.command)

    def _build_directory_tree_for_stabilize_solubis_solubismild_cysscan(self):
        self._create_Runs_subfolder_for('Stabilize', 'Solubis', 'SolubisMild', 'CysScan')

    def _create_Runs_subfolder_for(self, *args):
        for command in args:
            os.makedirs('Runs/' + command)

    def _compute_commands(self):
        if self.command == 'All':
            self._compute_stabilize_solubis_solubismild_cysscan()
        else:
            self._compute(self.command)

    def _compute_stabilize_solubis_solubismild_cysscan(self):
        self._compute('Stabilize', 'Solubis', 'SolubisMild', 'CysScan')

    def _compute(self, *args):
        for command in args:
            python_script = self._convert_command_name_to_python_script_name(command)
            if python_script == 'solubis.py':
                solubis_instance = Solubis(self.underscore_separated_proteinChains)
                solubis_instance.run_FoldX_BuildModel_all_APRs_GKs_scanning_mutations()
                solubis_instance.run_FoldX_AnalyseComplex_And_Agadir()
                solubis_instance.write_summary_solubis_file()
                self._print_OptProt_calling_script(python_script + ' (but via creating an object)')
            elif self.command == 'Supercharge' or self.command == 'DelPos' or self.command == 'Indiv':
                subprocess.call(
                    self.execute_python_and_path_to_script_fwdslash + python_script + self.space +
                    self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir_and_charge, shell=True)
            else:
                subprocess.call(
                    self.execute_python_and_path_to_script_fwdslash + python_script + self.space +
                    self.underscore_separated_proteinChains_and_paths_to_r_foldx_agadir, shell=True)

    def _convert_command_name_to_python_script_name(self, command):
        python_script_name = command
        if self.command == 'Stabilize':
            python_script_name = 'stabilize'
        elif self.command == 'Solubis':
            python_script_name = 'solubis'
        elif self.command == 'SolubisMild':
            python_script_name = 'solubis_mild'
        elif self.command == 'SolubisDouble':
            python_script_name = 'solubis_double'
        elif self.command == 'Supercharge':
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
