# Get Tango results, parse stretches, get resnums, create fastafiles for agadirwrapper, get difference total tango
import os
import sys
import glob
import subprocess
import GeneralUtilityMethods
from OptimizeProtein import yasara
import yaml

protein_chains = sys.argv[1].split('_')

yasara.info.mode = 'txt'

r_path = ''
foldx_path = ''
agadir_path = ''
results_path = ''
results_pdb_path = ''
aa_dict_1to3 = {}
aa_dict_3to1 = {}
gatekeepers = []
__solubis_job_prefix__ = 'SB_'
__analyze_complex_job_prefix__ = 'AC_'


with open("/switchlab/group/shazib/OptimizeProteinShazibCopy/SourceFiles/Scripts/pathsAndDictionaries.yaml",
          'r') as stream:
    try:

        paths_and_dictionaries = yaml.load(stream)
        r_path = paths_and_dictionaries['ROOT']['R_path']
        foldx_path = paths_and_dictionaries['ROOT']['FoldX_path']
        agadir_path = paths_and_dictionaries['ROOT']['Agadir_path']
        results_path = paths_and_dictionaries['ROOT']['results_path']
        aa_dict_1to3 = paths_and_dictionaries['ROOT']['aa_dict_1to3']
        aa_dict_3to1 = paths_and_dictionaries['ROOT']['aa_dict_3to1']
        gatekeepers = paths_and_dictionaries['ROOT']['gatekeepers']

    except yaml.YAMLError as exc:
        print(exc)


path_repaired_pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]

# def _get_pdb_name_from_repaired_pdb_folder():
# repaired_pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
# return repaired_pdb.split('/')[-1].split('.')[0].split('_')[-1]

pdb_name = path_repaired_pdb.split('/')[-1].split('.')[0].split('_')[-1] # removes the path and the repair prefix too
gatekeepers = ['R']  # for testing purposes only (to speed up solubis)
print pdb_name
# agadirs = []   can be deleted - NOT USED
globaltotal_Tango_all_chains_sum = 0
# def change_dir_to_pdb_Run_Solubis():
results_pdb_path = results_path + '/' + pdb_name
results_pdb_Runs_Solubis_path = results_pdb_path + '/Runs/Solubis'
os.chdir(results_pdb_Runs_Solubis_path)

##################################  FOR EACH CHAIN, FIND APR STRETCHES IN THIS REPAIRED PDB FILE  ##################################

individual_list_all_mutants = open('individual_list.txt', 'w')  # NOT SURE WHAT THIS LIST IS USED FOR
index_program = 0

for protein_chain in protein_chains:
################################## FOR EACH CHAIN, DIG OUT THE TANGO WINDOW SCORE ##################################
    tangowindow_file = open(results_pdb_path + '/Agadir/' + pdb_name + '_' + protein_chain + '/PSX_tangowindow.out', 'r').readlines()
############### FOR EACH LINE IN TANGO WINDOW DIG OUT THE APR and ITS POSITION IN THE FULL LENGTH PROTEIN ##########
    for line in tangowindow_file[1:]:
        # getTangoScore(line): try return line.split()[-2]
        apr = line.split()[-4]
        apr_index = int(float(line.split()[-6]))
#### FOR EACH RESIDUE IN THE APR, CREATE A MUTATION FOR EVERY GATEKEEPER. CREATE THE DIRECTORY STRUCTURE FOR EACH GK
        # MUTATION, AND USE BUILDMODEL TO MUTATE THE REPAIRED PDB TO PRODUCE A MUTANT PDB - ####
        for i, amino_acid in enumerate(apr):
            mutation_position = str(apr_index + i + 1)  # NOT SURE WHY THIS NEEDS TO BE MADE INTO A STRING
            for gatekeeper in gatekeepers:
                mutation_name = amino_acid + protein_chain + mutation_position + gatekeeper
                if not os.path.exists(mutation_name):
                    os.makedirs(mutation_name)
                os.chdir(mutation_name)
                if not os.path.exists('Fasta'):
                    os.makedirs('Fasta')
                if not os.path.exists('Agadir'):
                    os.makedirs('Agadir')
                if not os.path.exists('Agadir/Options.txt'):
                    subprocess.call('cp ' + results_pdb_path + '/../../SourceFiles/AgadirFiles/* ./Agadir/.',
                                    shell=True)
                if os.path.exists(results_pdb_path + '/Repair/RepairPDB_' + pdb_name + '.pdb'):
                    subprocess.call('cp ' + results_pdb_path + '/Repair/RepairPDB_' + pdb_name + '.pdb .', shell=True)
                    subprocess.call('cp ' + results_pdb_path + '/../../SourceFiles/FoldXFiles/* .', shell=True)
                else:
                    print 'Something is wrong'

                # def _mutate_pdb_with_FoldX_BuildModel():
                path_to_runscript = './'
                repaired_pdbs = 'RepairPDB_' + pdb_name + '.pdb'
                show_sequence_detail = False
                action = '<BuildModel>#,individual_list.txt'
                print_networks = False
                calculate_stability = False
                GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, repaired_pdbs,
                                                                   show_sequence_detail, action, print_networks,
                                                                   calculate_stability)

                individual_list_for_this_mutant_only = open('individual_list.txt', 'w')
                individual_list_for_this_mutant_only.write(mutation_name + ';\n')
                individual_list_for_this_mutant_only.close()

                grid_engine_job_name = __solubis_job_prefix__ + mutation_name
                no_queue = ''
                no_max_memory = ''
                no_cluster = ''
                using_runscript = True
                no_python_script = ''
                GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                           using_runscript, foldx_path, no_python_script)
                individual_list_all_mutants.write(mutation_name + ';\n')
                subprocess.call('qsub job.q', shell=True)
                os.chdir('./..')

individual_list_all_mutants.close()  # NOT SURE WHAT THIS LIST IS USED FOR
GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__solubis_job_prefix__, 'all Solubis jobs to finish')

##################################  PERFORM ANALYSE COMPLEX  ##################################

############# GET ALL THE SUBDIRECTORIES, THESE ARE THE DIRECTORIES MADE ABOVE FOR THE BUILDMODEL ALGORITHM, WHICH WILL
####### CONTAIN THE NEWLY BUILT REPAIRED PDBS FOR THE GATEKEEPER MUTANTS ################
# ListOfSolubisMutationFolderNames
dirs = sorted(glob.glob('./*'))  # get list of all solubis mutation names
repaired_pdb_name = 'RepairPDB_' + pdb_name
print 'Name:\t' + repaired_pdb_name
# SolubisMutationResultsFolder
for path in dirs:
    if os.path.isdir(path):
        os.chdir(path.split('/')[-1])
        mutation = path.split('/')[-1]
        subprocess.call('cp runscript.txt runscript_build.txt', shell=True)
        subprocess.call('rm runscript.txt', shell=True)

        path_to_runscript = './'
        common_prefix = repaired_pdb_name + '_1_'
        WT_common_prefix = 'WT_' + common_prefix
        pdbs_to_analyse = common_prefix + '0.pdb,' + \
                          common_prefix + '1.pdb,' + \
                          common_prefix + '2.pdb,' + \
                          WT_common_prefix + '0.pdb,' + \
                          WT_common_prefix + '1.pdb,' + \
                          WT_common_prefix + '2.pdb,'
        show_sequence_detail = False
        action = '<AnalyseComplex>#'
        print_networks = False
        calculate_stability = False
        GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, pdbs_to_analyse,
                                                           show_sequence_detail, action, print_networks,
                                                           calculate_stability)

################## THIS PART OF CODE SEEMS TO SUGGEST THAT THE REPAIR ALGORITHM CAN CHANGE THE AMINO ACID SEQUENCE ??! ##################

# def convert_repairedPdb_to_FASTA():
# not sure why _1_0 is the pdb of choice
        repaired_1_0_pdb = repaired_pdb_name + '_1_0.pdb'
        print os.getcwd()
        fasta_destination_path = os.getcwd() + 'Fasta/'
        # pdb2fasta_instance.convert_pdb_to_FASTA(repaired_pdb, fasta_destination_path)
        pdb_file = open(repaired_1_0_pdb).readlines()
        atomlines = []
        mols = []
        for line in pdb_file:
            if 'ATOM' == line[0:4]:
                mol = line[21]
                atomlines.append(line)
                if mol not in mols:
                    mols.append(mol)
        for mol in mols:
            fastalist = []
            resnum = '0'
            for line in atomlines:
                if line[21] == mol and resnum != line[22:26].strip(' '):
                    resnum = line[22:26].strip(' ')
                    aa = line[17:20]
                    fastalist.append(aa_dict_3to1[aa])
            fasta = "".join(fastalist)
            print repaired_pdb_name + '_' + mol
            print fasta
            # _write_repaired_pdb_name_in_fasta_file():
            repaired_pdb_fasta_file = open('Fasta/' + repaired_pdb_name + '_' + mol + '.fasta', 'w')
            repaired_pdb_fasta_file.write('>' + repaired_pdb_name + '_' + mol + '\n')
            repaired_pdb_fasta_file.write(fasta)
            repaired_pdb_fasta_file.close()

        grid_engine_job_name = __analyze_complex_job_prefix__ + mutation
        no_queue = ''
        no_max_memory = ''
        no_cluster = ''
        using_runscript = True
        python_script_with_path = results_pdb_path + '/../../SourceFiles/Scripts/agadir.py'
        GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                   using_runscript, foldx_path, python_script_with_path)

        subprocess.call('qsub job.q', shell=True)
        os.chdir('./..')
    else:
        print path

message_to_print = 'all AnalyseComplex jobs to finish'
GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__analyze_complex_job_prefix__, message_to_print)

################### WRITE A SUMMARY FILE WITH YOUR SOLUBIS RESULTS - not sure this summary text file actually contains anything ##################
############### FIRST IT LOOKS TO SEE IF ANALYSE COMPLEX FILE EXISTS (IN REPAIR FOLDER), I.E. I ASSUME THIS WAS CREATED ABOVE #########
os.chdir(results_pdb_path)
solubis_summary_file = open('SummarySolubis.txt', 'w')
solubis_summary_file.write('Mutation\tMol\tddG\tdTANGO\tComplexSum\t')
################### THE FIRST THING THE SUMMARY FILE WILL HAVE WRITTEN TO IT IS the list of Interactions between protein chains from the interaction_analysecomplex_pdb file)###################
analyse_complex = False
if os.path.isfile('./Repair/Interaction_AnalyseComplex_' + repaired_pdb_name + '.fxout'):  # WHAT IS THE ABSOLUTE PATH HERE ? THIS SEEMS TO BE IN Repair SUBFOLDER
    analyse_complex = True
    interact_ac_pdb_fxout_file = open('./Repair/Interaction_AnalyseComplex_' + repaired_pdb_name + '.fxout').readlines()
    interactions = []
    for line in interact_ac_pdb_fxout_file[9:]:
        pieces = line.split('\t')
        protein_chain_1 = line.split('\t')[1]  # does this work to be same as mol1 and mol2 ??
        protein_chain_2 = line.split('\t')[2]  # does this work to be same as mol1 and mol2 ??
        protein_chain_complex = protein_chain_1 + '_' + protein_chain_2
        interactions.append(protein_chain_complex)
    interaction_list = []
    for i in interactions:
        interaction_list.append('Complex_' + i)
    interaction_string = "\t".join(interaction_list)
    solubis_summary_file.write(interaction_string + '\n')
else:
    solubis_summary_file.write('\n')
dirs = sorted(glob.glob('./Runs/Solubis/*'))

agadir_outputs_wt_pdb_chains = glob.glob('./Agadir/*')  # does this ignore Options.txt ? Just the pdb_protein_chain folders?
TangoWT = 0
for agadir_outputs_wt_pdb_chain in agadir_outputs_wt_pdb_chains:  # agadir_outputs_pdb_chains
    if os.path.isdir(agadir_outputs_wt_pdb_chain):
        print agadir_outputs_wt_pdb_chain
        f = open(agadir_outputs_wt_pdb_chain + '/PSX_globaltotal.out', 'r').readlines()
        pieces = f[1].split()
        TangoMol = float(pieces[2])
        TangoWT = TangoWT + TangoMol
print 'TangoTotalWT = ' + str(TangoWT)

for path in dirs:  # ./Runs/Solubis/AH92R
    if os.path.isdir(path):  # returns True if path is an existing directory (? subdirectory of this current dir ?)
        # mut = path.split('/')[-1]  # AH92R
        # mol = mut[1]  # H
        f = open(path + '/Average_BuildModel_' + repaired_pdb_name + '.fxout', 'r').readlines()  # av_buildModel_repairedPdb_fxout_file
        ddG = f[9].split()[2]  # 13.3781
#  GO THROUGH THE AGADIR RESULTS FOLDERS FOR EACH CHAIN AND EXTRACT THE GLOBAL TANGO SCORE
        agadir_outputs_mutant_pdb_chains = glob.glob(path + '/Agadir/*')  # everything in Runs/Solubis/AH92R/Agadir/ (they are strings that include the relative path)
        TangoMut = 0  # how can this be 1630.285, surely it should be 0 !!?
        for agadir_outputs_mutant_pdb_chain in agadir_outputs_mutant_pdb_chains:
            if os.path.isdir(i): # isn't this redundant seeing as how path_agad_list is obtained the line above by getting everything in the Agadir folder.
                # def get_globalTango_from_PSX_globaltotalout_file(path)
                f = open(agadir_outputs_mutant_pdb_chain + '/PSX_globaltotal.out', 'r').readlines()  # globalTango_file
                pieces = f[1].split()
                TangoMol = float(pieces[2])  # ?
                TangoMut = TangoMut + TangoMol # ? 1630.. ?
                path_agad = agadir_outputs_mutant_pdb_chain  # ./Runs/Solubis/AH92R/Agadir/RepairPDB_Ab82b0sLigand_H  NOT SURE WHAT IT IS SHOWING CHAIN P !!?? It's because Agadir is outputting results for all the chains.. why ??
        print 'TangoTotalMut = ' + str(TangoMut)  # does a number need to be parsed to a string in order to print it out??
        print 'TangoTotalMut not as string... = ' + TangoMut  # does a number need to be parsed to a string in order to print it out??
        print path_agad
        f = open(path_agad + '/PSX_globaltotal.out', 'r').readlines()  # opening the globalTango file again ?
        print f[1]
        # mol = mut[1]
        if analyse_complex:
            common_prefix_AC = path + '/Interaction_AnalyseComplex_' + repaired_pdb_name + '_1_'
            common_prefix_AC_WT = path + '/Interaction_AnalyseComplex_WT_' + repaired_pdb_name + '_1_'
            f = open(common_prefix_AC + '0.fxout', 'r').readlines()
            complex_Mut1 = []
            for line in f[9:]:
                complex_Mut1.append(float(line.split()[5]))
            f = open(common_prefix_AC + '1.fxout', 'r').readlines()
            complex_Mut2 = []
            for line in f[9:]:
                complex_Mut2.append(float(line.split()[5]))
            f = open(common_prefix_AC + '2.fxout', 'r').readlines()
            complex_Mut3 = []
            for line in f[9:]:
                complex_Mut3.append(float(line.split()[5]))
            f = open(common_prefix_AC_WT + '0.fxout', 'r').readlines()
            complex_WT1 = []
            for line in f[9:]:
                complex_WT1.append(float(line.split()[5]))
            f = open(common_prefix_AC_WT + '1.fxout', 'r').readlines()
            complex_WT2 = []
            for line in f[9:]:
                complex_WT2.append(float(line.split()[5]))
            f = open(common_prefix_AC_WT + '2.fxout', 'r').readlines()
            complex_WT3 = []
            for line in f[9:]:
                complex_WT3.append(float(line.split()[5]))
            ComplexList = []
            ComplexSum = 0
            for x, line in enumerate(f[9:]):
                complex_WT = float((complex_WT1[x] + complex_WT2[x] + complex_WT3[x]) / 3)
                complex_Mut = float((complex_Mut1[x] + complex_Mut2[x] + complex_Mut3[x]) / 3)
                Complex = complex_Mut - complex_WT
                ComplexSum = ComplexSum + Complex
                ComplexList.append(str(Complex))
            Complex = "\t".join(ComplexList)
        else:
            Complex = ""
            ComplexSum = 0
        dTango = float(TangoMut) - float(TangoWT)
        mut = path.split('/')[-1]
        mol = mut[1]
        mol2 = path.split('/')[-1][1]
        sumstring = mut + '\t' + mol + '\t' + ddG + '\t' + str(dTango) + '\t' + str(ComplexSum) + '\t' + Complex + '\n'
        solubis_summary_file.write(sumstring)
solubis_summary_file.close()
