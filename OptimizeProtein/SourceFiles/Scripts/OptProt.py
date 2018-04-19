import os
import glob
import subprocess
import time
from OptimizeProtein import yasara
# import pydevd
# pydevd.settrace('localhost', port=51234, stdoutToServer=True, stderrToServer=True)

yasara.info.mode = 'txt'


class OptProt:
    # declare all global variables
    __r_path__ = ''
    __foldx_path__ = ''
    __agadir_path__ = ''
    __qsub_path__ = ''
    __start_path__ = ''
    __scripts_path__ = ''
    __execute_python_and_path_to_script_fwdslash__ = ''
    __molecules_and_paths_to_r_foldx_agadir__ = ''
    __molecules_and_paths_to_r_foldx_agadir_and_charge__ = ''

    def __init__(self, start_path, scripts_path):
        global __start_path__
        global __scripts_path__
        global __execute_python_and_path_to_script_fwdslash__
        __start_path__ = start_path
        __scripts_path__ = scripts_path
        __execute_python_and_path_to_script_fwdslash__ = 'python ' + __scripts_path__ + '/'
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
        self.__build_pdb_list_from(__option_file__)
        self.__parse_cmd_charge_mols_from(__option_file__)
        self.__parse_paths_for_r_foldx_agadir_qsub_from(__option_file__)
        global __molecules_and_paths_to_r_foldx_agadir__
        global __molecules_and_paths_to_r_foldx_agadir_and_charge__
        __molecules_and_paths_to_r_foldx_agadir__ = ' ' + self.__molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + \
            ' ' + __agadir_path__
        __molecules_and_paths_to_r_foldx_agadir_and_charge__ = __molecules_and_paths_to_r_foldx_agadir__ + ' ' + \
            self.__charge__

    # Called by parse_option_file() for:
    # 1. list of pdb names
    def __build_pdb_list_from(self, __option_file__):
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

    # Called by parse_option_file() for:
    # 2. computation names ("command"), charge, protein chains ("molecules")
    def __parse_cmd_charge_mols_from(self, __option_file__):
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

    # Called by parse_option_file() for:
    # 3. absolute paths to all required software and to the grid engine executables
    def __parse_paths_for_r_foldx_agadir_qsub_from(self, __option_file__):
        global __r_path__
        global __foldx_path__
        global __agadir_path__
        global __qsub_path__
        for line in __option_file__:
            if 'R_Path:' in line:
                __r_path__ = line.split(':')[-1].strip(';\n')
            if 'FoldX_Path:' in line:
                __foldx_path__ = line.split(':')[-1].strip(';\n')
            if 'Agadir_Path:' in line:
                __agadir_path__ = line.split(':')[-1].strip(';\n')
            if 'Qsub_Path' in line:
                __qsub_path__ = line.split(':')[-1].strip(';\n')
        print 'Absolute path to R:\t\t' + __r_path__
        print 'Absolute path to FoldX:\t\t' + __foldx_path__
        print 'Absolute path to TANGO:\t\t' + __agadir_path__
        print 'Absolute path to Qsub:\t\t' + __qsub_path__

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
            # current directory is now inside Results folder
            self._run_yasara_to_organise_pdb(PDB, pdb_name)
            self._copy_pdb_foldx_agadir_files_to_new_subdirectories(PDB)
            self._print_OptProt_calling_script('pdb2fasta.py')
            subprocess.call('python ' + __scripts_path__ + '/pdb2fasta.py', shell=True)
            self._print_OptProt_calling_script('agadir.py')
            subprocess.call('python ' + __scripts_path__ + '/agadir.py', shell=True)
            self._run_repair_on_grid_engine(pdb_name)

    def _run_repair_on_grid_engine(self, pdb_name):
        repair_python_script = 'repair.py'
        self._print_OptProt_calling_script(repair_python_script)
        g = open('./job.q', 'w')
        g.write('#!/bin/bash\n')
        g.write('#$ -N RPjob_' + pdb_name + '\n')
        g.write('#$ -V\n')
        g.write('#$ -cwd\n')
        g.write('source ~/.bash_profile\n')
        g.write(__execute_python_and_path_to_script_fwdslash__ + repair_python_script + ' ' + __qsub_path__ + '\n')
        g.close()
        subprocess.call(__qsub_path__ + 'qsub job.q', shell=True)
        os.chdir(__start_path__)

    def _copy_pdb_foldx_agadir_files_to_new_subdirectories(self, PDB):
        cp_start_path = 'cp ' + __start_path__
        subprocess.call(cp_start_path + '/PDBs/' + PDB + ' ./PDBs/.', shell=True)
        subprocess.call(cp_start_path + '/PDBs/' + PDB + ' ./Repair/.', shell=True)
        subprocess.call(cp_start_path + '/SourceFiles/FoldXFiles/* ./Repair/.', shell=True)
        subprocess.call(cp_start_path + '/SourceFiles/AgadirFiles/* ./Agadir/.', shell=True)

    def wait_for_repair_to_complete(self):
        check_qstat = subprocess.Popen(__qsub_path__ + 'qstat', stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
        while 'RPjob_' in output_qstat:
            print 'Waiting for PDBs to be repaired'
            time.sleep(10)
            check_qstat = subprocess.Popen(__qsub_path__ + 'qstat', stdout=subprocess.PIPE)
            output_qstat = check_qstat.stdout.read()

    def perform_selected_computations(self):
        for PDB in self.__pdb_list__:
            self._build_results_directory_tree_for_each(PDB)
            self._compute_stretchplot(PDB)
            if not os.path.exists('Runs/' + self.__command__):
                self._build_directory_tree_for_computations()
            self._compute_commands()
            os.chdir(__start_path__)

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

    def _compute_stretchplot(self, PDB):
        self._print_OptProt_calling_script('stretchplot.py')
        subprocess.call(
            'python ' + __scripts_path__ + '/stretchplot.py ' + __molecules_and_paths_to_r_foldx_agadir__,
            shell=True)
        self._print_OptProt_calling_script('stretchplot.R')
        subprocess.call('R < ' + __scripts_path__ + '/stretchplot.R --no-save', shell=True)

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
                    __molecules_and_paths_to_r_foldx_agadir_and_charge__,
                    shell=True)
            else:
                subprocess.call(
                    __execute_python_and_path_to_script_fwdslash__ + python_script + __molecules_and_paths_to_r_foldx_agadir__,
                    shell=True)

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

    def _run_yasara_to_organise_pdb(self, PDB, pdb_name):
        yasara.run('DelObj all')
        yasara.run('LoadPDB ' + __start_path__ + '/PDBs/' + PDB)
        yasara.run('DelRes !Protein')
        tempMols = yasara.run('ListMol All,Format=MOLNAME')
        for mol in tempMols:
            yasara.run('RenumberRes all and Mol ' + mol + ',First=1')
        yasara.run(
            'SavePDB 1,' + __start_path__ + '/Results/' + pdb_name + '/PDBs/' + PDB + ',Format=PDB,Transform=Yes')

    def _print_OptProt_calling_script(self, python_script):
        print 'OptProt.py calling ' + python_script + '....'


# pydevd.stoptrace()
