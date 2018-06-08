# Get Tango results, parse stretches, get resnums, create fastafiles for agadirwrapper, get difference total tango
import os
import glob
import subprocess
import GeneralUtilityMethods
from OptimizeProtein import yasara
import yaml


yasara.info.mode = 'txt'


class Solubis(object):

    def __init__(self, protein_chains):
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
        self.protein_chains = protein_chains
        self.solubis_job_prefix = 'SB_'
        self.analyze_complex_job_prefix = 'AC_'
        self.Int_AnalyComp = 'Interaction_AnalyseComplex_'
        self.gatekeepers = ['R']  # for testing purposes only (to speed up solubis)
        self.self.pdb_name_of_repaired_pdb = self._get_pdb_name_from_repaired_pdb_folder(self)
        print self.self.pdb_name_of_repaired_pdb
        self.results_pdb_Runs_Solubis_path = ''
        self.run_FoldX_BuildModel_all_APRs_GKs_scanning_mutations(self)

    def _get_pdb_name_from_repaired_pdb_folder(self):
        path_to_repaired_pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
        return path_to_repaired_pdb.split('/')[-1].split('.')[0].split('_')[-1]

    def _change_dir_to_pdb_Run_Solubis(self):
        self.results_pdb_path = self.results_path + '/' + self.self.pdb_name_of_repaired_pdb
        self.results_pdb_Runs_Solubis_path = self.self.results_pdb_path + '/Runs/Solubis'
        os.chdir(self.results_pdb_Runs_Solubis_path)

def run_FoldX_BuildModel_all_APRs_GKs_scanning_mutations(self):
    self._change_dir_to_pdb_Run_Solubis(self)
    individual_list_all_mutants = open('individual_list.txt', 'w')  # NOT SURE WHAT THIS LIST IS USED FOR
    for protein_chain in self.protein_chains:
        ##### FOR EACH CHAIN, GET ITS TANGO WINDOW RESULTS , FROM  WHICH IT SHOULD GET EACH APR SEQUENCE AND START POSITION #####
        tangowindow_file = open(self.results_pdb_path + '/Agadir/' + self.pdb_name_of_repaired_pdb + '_' + protein_chain + '/PSX_tangowindow.out', 'r').readlines()
        for line in tangowindow_file[1:]:
            apr_sequence = line.split()[-4]
            apr_index = int(float(line.split()[-6]))
            #### FOR EACH APR, MUTATE EVERY RESIDUE TO EVERY GATEKEEPER AND RUN BUILD MODEL TO CREATE REPAIRED MUTANT FILES ####
            for i, amino_acid in enumerate(apr_sequence):
                mutation_position = apr_index + i + 1
                for gatekeeper in self.gatekeepers:
                    # def _build_mutant_fasta_agadir_directory_structure(gatekeeper)
                    mutant_name = amino_acid + protein_chain + str(mutation_position) + gatekeeper
                    if not os.path.exists(mutant_name):
                        os.makedirs(mutant_name)
                    os.chdir(mutant_name)
                    if not os.path.exists('Fasta'):
                        os.makedirs('Fasta')
                    if not os.path.exists('Agadir'):
                        os.makedirs('Agadir')
                    if not os.path.exists('Agadir/Options.txt'):
                        subprocess.call('cp ' + self.results_pdb_path + '/../../SourceFiles/AgadirFiles/* ./Agadir/.',
                                        shell=True)
                    if os.path.exists(self.results_pdb_path + '/Repair/RepairPDB_' + self.pdb_name_of_repaired_pdb + '.pdb'):
                        subprocess.call('cp ' + self.results_pdb_path + '/Repair/RepairPDB_' + self.pdb_name_of_repaired_pdb + '.pdb .', shell=True)
                        subprocess.call('cp ' + self.results_pdb_path + '/../../SourceFiles/FoldXFiles/* .', shell=True)
                    else:
                        print 'Something is wrong'
                    # def _run_FoldX_BuildModel(mutation_name):
                    path_to_runscript = './'
                    repaired_pdbs = 'RepairPDB_' + self.pdb_name_of_repaired_pdb + '.pdb'
                    show_sequence_detail = False
                    action = '<BuildModel>#,individual_list.txt'
                    print_networks = False
                    calculate_stability = False
                    GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, repaired_pdbs,
                                                                       show_sequence_detail, action, print_networks,
                                                                       calculate_stability)
                    individual_list_for_this_mutant_only = open('individual_list.txt', 'w')
                    individual_list_for_this_mutant_only.write(mutant_name + ';\n')
                    individual_list_for_this_mutant_only.close()

                    grid_engine_job_name = self.solubis_job_prefix + mutant_name
                    no_queue = ''
                    no_max_memory = ''
                    no_cluster = ''
                    using_runscript = True
                    no_python_script = ''
                    GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                               using_runscript, self.foldx_path, no_python_script)
                    individual_list_all_mutants.write(mutant_name + ';\n')
                    subprocess.call('qsub job.q', shell=True)
                    os.chdir('./..')  # maybe pass the absolute path here
    individual_list_all_mutants.close()  # NOT SURE WHAT THIS LIST IS USED FOR
    GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(self.solubis_job_prefix, 'all Solubis jobs to finish')
    self.run_FoldX_AnalyseComplex_(self)

    def run_FoldX_AnalyseComplex(self):
        self.all_mutants_in_Solubis_path = sorted(glob.glob('./*'))
        repair_pdb_name = 'RepairPDB_' + self.pdb_name_of_repaired_pdb
        repair_pdb_name_1_ = repair_pdb_name + '_1_'
        wt_repair_pdb_name_1_ = 'WT_' + repair_pdb_name_1_
        _0_1_2_pdbs = ['0.pdb,', '1.pdb,', '2.pdb,']

        for mutant_in_Solubis_path in self.all_mutants_in_Solubis_path:
            if os.path.isdir(mutant_in_Solubis_path):
                mutant_folder_name = mutant_in_Solubis_path.split('/')[-1]
                os.chdir(mutant_folder_name)
                # def _make_backup_of_runscript():
                subprocess.call('cp runscript.txt runscript_build.txt', shell=True)
                subprocess.call('rm runscript.txt', shell=True)
                path_to_runscript = './'
                pdbs_to_analyse = repair_pdb_name_1_ + _0_1_2_pdbs[0] + \
                                  repair_pdb_name_1_ + _0_1_2_pdbs[1] + \
                                  repair_pdb_name_1_ + _0_1_2_pdbs[2] + \
                                  wt_repair_pdb_name_1_ + _0_1_2_pdbs[0] + \
                                  wt_repair_pdb_name_1_ + _0_1_2_pdbs[1] + \
                                  wt_repair_pdb_name_1_ + _0_1_2_pdbs[2]
                show_sequence_detail = False
                action = '<AnalyseComplex>#'
                print_networks = False
                calculate_stability = False
                GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, pdbs_to_analyse,
                                                                   show_sequence_detail, action, print_networks,
                                                                   calculate_stability)

                print os.getcwd()
                # Note: _1_0 is random choice, other two are same.
                # cwd is Results/pdbname/Runs/Solubis/mutant_name
                GeneralUtilityMethods.GUM.extract_fasta_from_pdb(repair_pdb_name + '_1_0.pdb', './')

                grid_engine_job_name = self.analyze_complex_job_prefix + mutant_folder_name
                no_queue = ''
                no_max_memory = ''
                no_cluster = ''
                using_runscript = True
                python_script_with_path = self.results_pdb_path + '/../../SourceFiles/Scripts/agadir.py'
                GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                           using_runscript, self.foldx_path, python_script_with_path)

                subprocess.call('qsub job.q', shell=True)
                os.chdir('./..')  # is this Runs/Solubis folder ?
            else:
                print mutant_in_Solubis_path

        message_to_print = 'all AnalyseComplex jobs to finish'
        GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(self.analyze_complex_job_prefix, message_to_print)
        self.write_summary_solubis_file(self)

    def write_summary_solubis_file(self):
        os.chdir(self.results_pdb_path)
        solubis_summary_file = open('SummarySolubis.txt', 'w')
        solubis_summary_file.write('Mutation\tProteinChahin\tddG\tdTANGO\tComplexSum\t')
        has_I_AC_fxout_file_for_repaired_pdb = False
        if os.path.isfile('./Repair/' + self.Int_AnalyComp + self.repair_pdb_name + '.fxout'):  # perhaps use absolute path?
            has_I_AC_fxout_file_for_repaired_pdb = True
            interact_ac_pdb_fxout_file = open('./Repair/' + self.Int_AnalyComp + self.repair_pdb_name + '.fxout').readlines()
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
                f = open(solubis_mutant_folder_path + '/Average_BuildModel_' + self.repair_pdb_name + '.fxout', 'r').readlines()
                ddG = f[9].split()[2]  # THIS IS THE FOLDX STABILITY VALUE
                all_agadir_output_files_for_this_mutant_chain_paths = glob.glob(solubis_mutant_folder_path + '/Agadir/*')
                TangoMut = 0
                for agadir_output_files_for_this_mutant_chain_path in all_agadir_output_files_for_this_mutant_chain_paths:
                    if os.path.isdir(agadir_output_files_for_this_mutant_chain_path): # isn't this redundant cos path_agad_list is obtained the line above by getting everything in the Agadir folder.
                        # def get_globalTango_from_PSX_globaltotalout_file(path):
                        f = open(agadir_output_files_for_this_mutant_chain_path + '/PSX_globaltotal.out', 'r').readlines()  # globalTango_file
                        TangoMut += float(f[1].split()[2])
                        path_agad = agadir_output_files_for_this_mutant_chain_path  # ./Runs/Solubis/AH92R/Agadir/RepairPDB_Ab82b0sLigand_H
                print 'Tango GlobalTotal score for mutant protein = ' + str(TangoMut)
                print path_agad  # don't understand this. It goes through all the chains but only stores the last one to open the tango file on line below
                f = open(path_agad + '/PSX_globaltotal.out', 'r').readlines()  # opening the globalTango file again ?
                print f[1]  # which it then just prints and doesn't use ? - is this leftover from debugging
                if has_I_AC_fxout_file_for_repaired_pdb:
                    path_I_AC_repaired_pdb_name_1_ = solubis_mutant_folder_path + '/' + self.Int_AnalyComp + self.repair_pdb_name_1_
                    path_I_AC_wt_repaired_pdb_name_1_ = solubis_mutant_folder_path + '/' + self.Int_AnalyComp + self.wt_repair_pdb_name_1_
                    _0_1_2_fxout = ['0.fxout', '1.fxout', '2.fxout']
                    # for fxout in _0_1_2_fxout:
                    #     _get_interaction_energies_from_I_AC_file(path_IAC_repaired_pdb_name_1_ + fxout)
                    # complex_Mut1 = _get_interaction_energies_from_I_AC_file(path_IAC_repaired_pdb_name_1_ + _0_1_2_fxout[0])
                    f = open(path_I_AC_repaired_pdb_name_1_ + _0_1_2_fxout[0], 'r').readlines()
                    complex_Mut1 = []
                    for line in f[9:]:
                        complex_Mut1.append(float(line.split()[5]))
                    f = open(path_I_AC_repaired_pdb_name_1_ + _0_1_2_fxout[1], 'r').readlines()
                    complex_Mut2 = []
                    for line in f[9:]:
                        complex_Mut2.append(float(line.split()[5]))
                    f = open(path_I_AC_repaired_pdb_name_1_ + _0_1_2_fxout[2], 'r').readlines()
                    complex_Mut3 = []
                    for line in f[9:]:
                        complex_Mut3.append(float(line.split()[5]))
                    f = open(path_I_AC_wt_repaired_pdb_name_1_ + _0_1_2_fxout[0], 'r').readlines()
                    complex_WT1 = []
                    for line in f[9:]:
                        complex_WT1.append(float(line.split()[5]))
                    f = open(path_I_AC_wt_repaired_pdb_name_1_ + _0_1_2_fxout[1], 'r').readlines()
                    complex_WT2 = []
                    for line in f[9:]:
                        complex_WT2.append(float(line.split()[5]))
                    f = open(path_I_AC_wt_repaired_pdb_name_1_ + _0_1_2_fxout[2], 'r').readlines()
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
