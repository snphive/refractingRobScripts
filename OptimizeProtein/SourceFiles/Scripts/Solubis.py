# Get Tango results, parse stretches, get resnums, create fastafiles for agadirwrapper, get difference total tango
from OptimizeProtein import yasara
import glob
import os
import subprocess
import yaml
from GeneralUtilityMethods import GUM
from Agadir import Agadir
from FoldX import FoldX


yasara.info.mode = 'txt'


class Solubis(object):

    def __init__(self, underscore_separated_protein_chains):
        with open("/switchlab/group/shazib/OptimizeProteinShazibCopy/SourceFiles/Scripts/pathsAndDictionaries.yaml",
                  'r') as stream:
            try:
                paths_and_dictionaries = yaml.load(stream)
                self.r_path = paths_and_dictionaries['ROOT']['R_path']
                self.foldx_path = paths_and_dictionaries['ROOT']['FoldX_path']
                self.agadir_path = paths_and_dictionaries['ROOT']['Agadir_path']
                self.results_path = paths_and_dictionaries['ROOT']['results_path']
                self.aa_dict_1to3 = paths_and_dictionaries['ROOT']['aa_dict_1to3']
                self.aa_dict_3to1 = paths_and_dictionaries['ROOT']['aa_dict_3to1']
                self.gatekeepers = paths_and_dictionaries['ROOT']['gatekeepers']
            except yaml.YAMLError as exc:
                print(exc)
        self.protein_chains_list = underscore_separated_protein_chains.split('_')
        self.solubis_job_prefix = 'SB_'
        self.analyze_complex_job_prefix = 'AC_'
        self.Int_AnalyComp = 'Interaction_AnalyseComplex_'
        self.gatekeepers = ['R']  # for testing purposes only (to speed up solubis)
        self.name_of_repair_pdb = self._get_pdb_name_from_repair_folder()
        self.RepairPDB_pdb_name = 'RepairPDB_' + self.name_of_repair_pdb

        print self.name_of_repair_pdb
        self.results_pdb_Runs_Solubis_path = ''
        self.all_mutants_in_Solubis_path = ''
        self.individual_list_all_mutants = open('individual_list.txt', 'w')  # NOT SURE WHAT THIS LIST IS USED FOR

    def _get_pdb_name_from_repair_folder(self):
        path_to_repair_pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
        return path_to_repair_pdb.split('/')[-1].split('.')[0].split('_')[-1]

    def _change_dir_to_pdb_Run_Solubis(self):
        self.results_pdb_path = self.results_path + '/' + self.name_of_repair_pdb
        self.results_pdb_Runs_Solubis_path = self.results_pdb_path + '/Runs/Solubis'
        os.chdir(self.results_pdb_Runs_Solubis_path)

    def run_FoldX_BuildModel_all_APRs_GKs_scanning_mutations(self):
        self._change_dir_to_pdb_Run_Solubis()
        foldx_instance = FoldX()
        for protein_chain in self.protein_chains_list:
            tangowindow_file = open(self.results_pdb_path + '/Agadir/' + self.name_of_repair_pdb + '_' +
                                    protein_chain + '/PSX_tangowindow.out', 'r').readlines()
            for line in tangowindow_file[1:]:
                apr_sequence_and_index = Agadir.get_APR_sequence_and_index_from_line_in_tangowindow_file(line)
                for i, amino_acid in enumerate(apr_sequence_and_index[0]):
                    mutation_position = apr_sequence_and_index[1] + i + 1
                    for gatekeeper in self.gatekeepers:
                        mutant_name = amino_acid + protein_chain + str(mutation_position) + gatekeeper
                        self._build_fasta_agadir_directory_structure_for_mutant(mutant_name)
                        self._copy_Agadir_options_FoldX_files_repair_pdb_to_mutant_folder(mutant_name,
                                                                                          self.name_of_repair_pdb)
                        foldx_instance.prepare_for_FoldX_BuildModel(mutant_name, self.name_of_repair_pdb)
                        self.individual_list_all_mutants.write(mutant_name + ';\n')
                        grid_engine_job_name = self.solubis_job_prefix + mutant_name
                        python_script = ''
                        self._run_runscript_on_grid_engine(grid_engine_job_name, python_script)
                        os.chdir('./..')  # maybe pass the absolute path here ?
        self.individual_list_all_mutants.close()  # NOT SURE WHAT THIS LIST IS USED FOR

    def run_FoldX_AnalyseComplex_And_Agadir(self):
        GUM.wait_for_grid_engine_job_to_complete(self.solubis_job_prefix, 'all Solubis jobs to finish')
        self.all_mutants_in_Solubis_path = sorted(glob.glob('./*'))
        foldx_instance = FoldX()
        for mutant_in_Solubis_path in self.all_mutants_in_Solubis_path:
            if os.path.isdir(mutant_in_Solubis_path):
                mutant_folder_name = mutant_in_Solubis_path.split('/')[-1]
                os.chdir(mutant_folder_name)
                self._make_backup_of_runscript()
                foldx_instance.prepare_for_FoldX_AnalyseComplex(self.RepairPDB_pdb_name)
                print os.getcwd()
                # Note: _1_0 is random choice, other two are same.
                # cwd is Results/pdbname/Runs/Solubis/mutant_name
                relative_path_of_pdb_to_read = './'
                relative_path_for_new_fasta_folder_to_write = './'
                pdb_name_chain_fasta_dict = GUM.extract_pdb_name_fasta_and_chains_from_pdb(self.RepairPDB_pdb_name +
                                                                            '_1_0.pdb', relative_path_of_pdb_to_read)
                GUM.write_fasta_to_folder(pdb_name_chain_fasta_dict, relative_path_for_new_fasta_folder_to_write)
                grid_engine_job_name = self.analyze_complex_job_prefix + mutant_folder_name
                agadir_results_path = os.getcwd()
                python_script_with_path = self.results_pdb_path + '/../../SourceFiles/Scripts/run_agadir.py ' + \
                                          agadir_results_path
                self._run_runscript_on_grid_engine(grid_engine_job_name, python_script_with_path)
                os.chdir('./..')  # is this Runs/Solubis folder ?
                cwd = os.getcwd()  # for debugging
            else:
                print mutant_in_Solubis_path

    def write_summary_solubis_file(self):
        GUM.wait_for_grid_engine_job_to_complete(self.analyze_complex_job_prefix, 'all AnalyseComplex jobs to finish')
        os.chdir(self.results_pdb_path)
        solubis_summary_file = open('SummarySolubis.txt', 'w')
        solubis_summary_file.write('Mutation\tProteinChahin\tddG\tdTANGO\tComplexSum\t')
        has_I_AC_fxout_file_for_repair_pdb = False
        if os.path.isfile('./Repair/' + self.Int_AnalyComp + self.RepairPDB_pdb_name + '.fxout'):  # use absolute path?
            has_I_AC_fxout_file_for_repair_pdb = True
            interact_ac_pdb_fxout_file = open('./Repair/' + self.Int_AnalyComp + self.RepairPDB_pdb_name + '.fxout').readlines()
            protein_chain_complexes = []
            for line in interact_ac_pdb_fxout_file[9:]:
                protein_chain_1 = line.split('\t')[1]
                protein_chain_2 = line.split('\t')[2]
                protein_chain_complex = 'Complex_' + protein_chain_1 + '_' + protein_chain_2
                protein_chain_complexes.append(protein_chain_complex)
                tab_separated_protein_chain_complexes = "\t".join(protein_chain_complexes)
            solubis_summary_file.write(tab_separated_protein_chain_complexes + '\n')
        else:
            solubis_summary_file.write('\n')

        all_agadir_outputs_wt_pdb_chains = glob.glob('./Agadir/*')  # use the absolute path here
        TangoWT = 0
        for agadir_outputs_wt_pdb_chain in all_agadir_outputs_wt_pdb_chains:  # agadir_outputs_pdb_chains
            if os.path.isdir(agadir_outputs_wt_pdb_chain):
                print agadir_outputs_wt_pdb_chain
                f = open(agadir_outputs_wt_pdb_chain + '/PSX_globaltotal.out', 'r').readlines()
                TangoWT += float(f[1].split()[2])

        print 'Tango GlobalTotal score for WT protein = ' + str(TangoWT)
        os.chdir(self.results_pdb_Runs_Solubis_path)
        for solubis_mutant_folder_path in self.all_mutants_in_Solubis_path:
            if os.path.isdir(solubis_mutant_folder_path):
                f = open(solubis_mutant_folder_path + '/Average_BuildModel_' + self.RepairPDB_pdb_name + '.fxout', 'r').readlines()
                # # THIS IS THE FOLDX STABILITY VALUE # #
                ddG = f[9].split()[2]
                agadir_output_files_for_this_mutant_chain_paths = glob.glob(solubis_mutant_folder_path + '/Agadir/*')
                TangoMut = 0
                for agadir_output_files_for_this_mutant_chain_path in agadir_output_files_for_this_mutant_chain_paths:
                    if os.path.isdir(agadir_output_files_for_this_mutant_chain_path):
                        # redundant ? path_agad_list is obtained line above by getting everything in Agadir folder
                        # def get_globalTango_from_PSX_globaltotalout_file(path):
                        f = open(agadir_output_files_for_this_mutant_chain_path +
                                 '/PSX_globaltotal.out', 'r').readlines()  # globalTango_file
                        TangoMut += float(f[1].split()[2])
                        path_agad = agadir_output_files_for_this_mutant_chain_path
                        # ./Runs/Solubis/AH92R/Agadir/RepairPDB_Ab82b0sLigand_H
                print 'Tango GlobalTotal score for mutant protein = ' + str(TangoMut)
                print path_agad
                # don't understand. Goes through all chains but only stores last one to open the tango file below
                f = open(path_agad + '/PSX_globaltotal.out', 'r').readlines()  # opening the globalTango file again ?
                print f[1]  # which it then just prints and doesn't use ? - is this leftover from debugging
                if has_I_AC_fxout_file_for_repair_pdb:
                    repair_pdb_name_1_ = self.RepairPDB_pdb_name + '_1_'
                    wt_repair_pdb_name_1_ = 'WT_' + repair_pdb_name_1_
                    path_I_AC_repair_pdb_name_1_ = solubis_mutant_folder_path + '/' + self.Int_AnalyComp + repair_pdb_name_1_
                    path_I_AC_wt_repair_pdb_name_1_ = solubis_mutant_folder_path + '/' + self.Int_AnalyComp + wt_repair_pdb_name_1_
                    _0_1_2_fxout = ['0.fxout', '1.fxout', '2.fxout']
                    # for fxout in _0_1_2_fxout:
                    #     _get_interaction_energies_from_I_AC_file(path_IAC_repair_pdb_name_1_ + fxout)
                    # complex_Mut1 = _get_interaction_energies_from_I_AC_file(path_IAC_repair_pdb_name_1_ + _0_1_2_fxout[0])
                    f = open(path_I_AC_repair_pdb_name_1_ + _0_1_2_fxout[0], 'r').readlines()
                    complex_Mut1 = []
                    for line in f[9:]:
                        complex_Mut1.append(float(line.split()[5]))
                    f = open(path_I_AC_repair_pdb_name_1_ + _0_1_2_fxout[1], 'r').readlines()
                    complex_Mut2 = []
                    for line in f[9:]:
                        complex_Mut2.append(float(line.split()[5]))
                    f = open(path_I_AC_repair_pdb_name_1_ + _0_1_2_fxout[2], 'r').readlines()
                    complex_Mut3 = []
                    for line in f[9:]:
                        complex_Mut3.append(float(line.split()[5]))
                    f = open(path_I_AC_wt_repair_pdb_name_1_ + _0_1_2_fxout[0], 'r').readlines()
                    complex_WT1 = []
                    for line in f[9:]:
                        complex_WT1.append(float(line.split()[5]))
                    f = open(path_I_AC_wt_repair_pdb_name_1_ + _0_1_2_fxout[1], 'r').readlines()
                    complex_WT2 = []
                    for line in f[9:]:
                        complex_WT2.append(float(line.split()[5]))
                    f = open(path_I_AC_wt_repair_pdb_name_1_ + _0_1_2_fxout[2], 'r').readlines()
                    complex_WT3 = []
                    for line in f[9:]:
                        complex_WT3.append(float(line.split()[5]))

                    dInteractionEnergies_per_complex_list = []
                    dInteractionEnergies_all_complexes_summed = 0
                    for x, line in enumerate(f[9:]):
                        average_interaction_energy_per_complex_WT = float((complex_WT1[x] + complex_WT2[x] + complex_WT3[x]) / 3)
                        average_interaction_energy_per_complex_Mut = float((complex_Mut1[x] + complex_Mut2[x] + complex_Mut3[x]) / 3)
                        dInteractionEnergies_for_each_complex = average_interaction_energy_per_complex_Mut - average_interaction_energy_per_complex_WT
                        dInteractionEnergies_all_complexes_summed += dInteractionEnergies_for_each_complex
                        dInteractionEnergies_per_complex_list.append(str(dInteractionEnergies_for_each_complex))
                    dInteractionEnergies_per_complex = "\t".join(dInteractionEnergies_per_complex_list)
                else:
                    dInteractionEnergies_per_complex = ""
                    dInteractionEnergies_all_complexes_summed = 0
                dTango = float(TangoMut) - float(TangoWT)
                mutation_name = solubis_mutant_folder_path.split('/')[-1]
                protein_chain = mutation_name[1]  # Does this give same as above two lines?
                row_of_values_per_mutation = mutation_name + '\t' + protein_chain + '\t' + ddG + '\t' + str(dTango) + '\t' + \
                                     str(dInteractionEnergies_all_complexes_summed) + '\t' + \
                                     dInteractionEnergies_per_complex + '\n'
                solubis_summary_file.write(row_of_values_per_mutation)
        solubis_summary_file.close()

    # def _get_interaction_energies_from_I_AC_file(fxout_file_name):
    #     fxout_file = open(fxout_file_name, 'r').readlines()
    #     interaction_energies = []
    #     for line in fxout_file[9:]:
    #         interaction_energy.append(float(line.split()[5]))
    #     return interaction_energies

    def _build_fasta_agadir_directory_structure_for_mutant(self, mutant_name):
        if not os.path.exists(mutant_name):
            os.makedirs(mutant_name)
        os.chdir(mutant_name)
        if not os.path.exists('Fasta'):
            os.makedirs('Fasta')
        if not os.path.exists('Agadir'):
            os.makedirs('Agadir')
        if not os.path.exists('Agadir') or not os.path.exists('Fasta'):
            print 'Something is wrong - Fasta or Agadir folders were not created'

    def _copy_Agadir_options_FoldX_files_repair_pdb_to_mutant_folder(self, mutant_name, pdb_name):
        if os.getcwd().split('/')[-1] != mutant_name:
            os.chdir(mutant_name)
        if not os.path.exists('Agadir/Options.txt'):
            subprocess.call('cp ' + self.results_pdb_path + '/../../SourceFiles/AgadirFiles/* ./Agadir/.', shell=True)
        if os.path.exists(self.results_pdb_path + '/Repair/RepairPDB_' + pdb_name + '.pdb'):
            subprocess.call('cp ' + self.results_pdb_path + '/Repair/RepairPDB_' + pdb_name + '.pdb .', shell=True)
            subprocess.call('cp ' + self.results_pdb_path + '/../../SourceFiles/FoldXFiles/* .', shell=True)
        else:
            print 'Something is wrong - RepairPDB_ pdb does not exist for: ' + pdb_name

    def _run_runscript_on_grid_engine(self, grid_engine_job_name, python_script_with_path):
        no_queue = ''
        no_max_memory = ''
        no_cluster = ''
        using_runscript = True
        GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                             using_runscript, self.foldx_path, python_script_with_path)
        subprocess.call('qsub job.q', shell=True)

    def _make_backup_of_runscript(self):
        subprocess.call('cp runscript.txt runscript_build.txt', shell=True)
        subprocess.call('rm runscript.txt', shell=True)
