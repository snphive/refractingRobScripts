import os
import glob
import subprocess
import time
from OptimizeProtein import yasara

# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

yasara.info.mode = 'txt'


# The original OptProt script was designed as a type of main script that is responsible for determining which
# calculations will be made on which pdbs, based on reading in and parsing the options file.
# It includes steps that are necessary precursors to any of our protein structural and aggregation propensity
# calculations. This includes utilising yasara to clean up the pdb file.
class OptProt:
    # declare all global variables
    __r_path__ = ''
    __foldx_path__ = ''
    __agadir_path__ = ''
    __qsub_path__ = ''
    __start_path__ = ''
    __scripts_path__ = ''
    __execute_python_and_path_to_script_fwdslash__ = ''
    __execute_r_and_path_to_script_fwdslash__ = ''
    __molecules_and_paths_to_r_foldx_agadir__ = ''
    __molecules_and_paths_to_r_foldx_agadir_and_charge__ = ''

    def __init__(self, start_path, scripts_path):
        global __start_path__
        global __scripts_path__
        global __execute_python_and_path_to_script_fwdslash__
        global __execute_r_and_path_to_script_fwdslash__
        __start_path__ = start_path
        __scripts_path__ = scripts_path
        self.__single_space__ = ' '
        __execute_python_and_path_to_script_fwdslash__ = 'python' + self.__single_space__ + __scripts_path__ + '/'
        __execute_r_and_path_to_script_fwdslash__ = 'R <' + self.__single_space__ + __scripts_path__ + '/'
        self.__pdb_list__ = []
        self.__command__ = ''
        self.__charge__ = ''
        self.__molecules__ = ''

    # Extracts instructions from the option file.
    # This includes which computations to run & which proteins to run them on in 3 parts:
    # 1. list of pdb names
    # 2. computation names ("command"), charge, protein chains ("molecules")
    # 3. absolute paths to all required software and to the grid engine executables
    def parse_option_file(self, __option_file__):
        self.__parse_optionfile_for_pdblist(__option_file__)
        self.__parse_optionfile_for_computations_charge_mols(__option_file__)
        self.__parse_optionfile_for_paths_to_r_foldx_agadir_qsub(__option_file__)
        global __molecules_and_paths_to_r_foldx_agadir__
        global __molecules_and_paths_to_r_foldx_agadir_and_charge__
        __molecules_and_paths_to_r_foldx_agadir__ = self.__single_space__ + self.__molecules__ + self.__single_space__ \
                                                    + __r_path__ + self.__single_space__ + __foldx_path__ + \
                                                    self.__single_space__ + __agadir_path__
        __molecules_and_paths_to_r_foldx_agadir_and_charge__ = __molecules_and_paths_to_r_foldx_agadir__ + \
                                                               self.__single_space__ + self.__charge__

    # Tidies up each pdb file using Yasara functions.
    # Converts each pdb to fasta protein sequence via pdb2fasta.py.
    # Runs Agadirwrapper on each fasta protein sequence via agadir.py.
    # Runs FoldX repair on each pdb via repair.py.
    def run_yasara_agadir_repair(self):
        if not os.path.exists('Results'):
            os.makedirs('Results')
        for PDB in self.__pdb_list__:
            pdb_name = PDB.split('.')[0]
            self._build_results_directory_tree_for_each(PDB)
            # current directory is the pdb subfolder of the Results folder
            self._run_yasara_to_organise_pdb(PDB, pdb_name)
            self._copy_pdb_foldx_agadir_files_to_new_subdirectories(PDB)
            self._run_agadir()
            self._run_repair_on_grid_engine(pdb_name)

    def perform_selected_computations(self):
        self._wait_for_repair_to_complete()
        for PDB in self.__pdb_list__:
            self._build_results_directory_tree_for_each(PDB)
            self._compute_stretchplot(PDB)
            if not os.path.exists('Runs/' + self.__command__):
                self._build_directory_tree_for_computations()
            self._compute_commands()
            os.chdir(__start_path__)

    # # # # # # # # # PRIVATE METHODS # # # # # # # # #

    # # # # # # Called by parse_option_file() # # # # #

    # 1. returns list of pdb names from option file.
    def __parse_optionfile_for_pdblist(self, __option_file__):
        self.__pdb_list__ = []
        for line in __option_file__:
            if '#' in line:
                continue
            if 'PDBs:' in line:
                if line.split(':')[-1].strip(';\n') == "All":
                    pdbpaths = glob.glob('./PDBs/*.pdb')
                    for pdbpath in pdbpaths:
                        pdb_temp = pdbpath.split('/')[-1]
                        self.__pdb_list__.append(pdb_temp)
                pdb_string = line.split(':')[-1].strip(';\n')
                if ',' in pdb_string:
                    pdbs_temp = pdb_string.split(',')
                    for pdb_temp in pdbs_temp:
                        self.__pdb_list__.append(pdb_temp)
                else:
                    self.__pdb_list__.append(pdb_string)
        print 'PDBs to analyse:\t\t' + ",\t".join(self.__pdb_list__)

    # 2. returns computation names ("command"), charge, protein chains ("molecules") from option file.
    def __parse_optionfile_for_computations_charge_mols(self, __option_file__):
        for line in __option_file__:
            if '#' in line:
                continue
            if 'Command:' in line:
                self.__command__ = line.split(':')[-1].strip(';\n')
            if 'Charge:' in line:
                self.__charge__ = line.split(':')[-1].strip(';\n')
            if 'Mols:' in line:
                self.__molecules__ = []
                MolString = line.split(':')[-1].strip(';\n')
                if ',' in MolString:
                    Mols_temp = MolString.split(',')
                    for Mols_temp in Mols_temp:
                        self.__molecules__.append(Mols_temp)
                    self.__molecules__ = "_".join(self.__molecules__)
                else:
                    self.__molecules__ = MolString
        print 'Command to be executed:\t\t' + self.__command__
        print 'Molecules to be considered:\t' + self.__molecules__

    # 3. returns absolute paths to all required software and to the grid engine executables from option file.
    def __parse_optionfile_for_paths_to_r_foldx_agadir_qsub(self, __option_file__):
        global __r_path__
        global __foldx_path__
        global __agadir_path__
        global __qsub_path__
        for line in __option_file__:
            if 'R_Path' in line:
                __r_path__ = line.split(':')[-1].strip(';\n')
            if 'FoldX_Path' in line:
                __foldx_path__ = line.split(':')[-1].strip(';\n')
            if 'Agadir_Path' in line:
                __agadir_path__ = line.split(':')[-1].strip(';\n')
            if 'Qsub_Path' in line:
                __qsub_path__ = line.split(':')[-1].strip(';\n')
        self._print_Absolute_path_to('R', __r_path__)
        self._print_Absolute_path_to('FoldX', __foldx_path__)
        self._print_Absolute_path_to('TANGO', __agadir_path__)
        self._print_Absolute_path_to('Qsub', __qsub_path__)

    # # # # Called by run_yasara_agadir_repair() # # # #

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
        repair_python_script = 'repair.py'
        self._print_OptProt_calling_script(repair_python_script)
        g = open('./job.q', 'w')
        g.write('#!/bin/bash\n')
        g.write('#$ -N RPjob_' + pdb_name + '\n')
        g.write('#$ -V\n')
        g.write('#$ -cwd\n')
        g.write('source ~/.bash_profile\n')
        g.write(
            __execute_python_and_path_to_script_fwdslash__ + repair_python_script + self.__single_space__ + __qsub_path__ + '\n')
        g.close()
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

    def _wait_for_repair_to_complete(self):
        check_qstat = subprocess.Popen(__qsub_path__ + 'qstat', stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
        while 'RPjob_' in output_qstat:
            print 'Waiting for PDBs to be repaired'
            time.sleep(10)
            check_qstat = subprocess.Popen(__qsub_path__ + 'qstat', stdout=subprocess.PIPE)
            output_qstat = check_qstat.stdout.read()

    def _compute_stretchplot(self, pdb):
        stretchplot_python_script = 'stretchplot.py'
        stretchplot_r_script = 'stretchplot.R'
        self._print_OptProt_calling_script(stretchplot_python_script)
        subprocess.call(
            __execute_python_and_path_to_script_fwdslash__ + stretchplot_python_script + self.__single_space__ +
            __molecules_and_paths_to_r_foldx_agadir__, shell=True)
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
                    __molecules_and_paths_to_r_foldx_agadir_and_charge__, shell=True)
            else:
                subprocess.call(
                    __execute_python_and_path_to_script_fwdslash__ + python_script +
                    __molecules_and_paths_to_r_foldx_agadir__, shell=True)

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
