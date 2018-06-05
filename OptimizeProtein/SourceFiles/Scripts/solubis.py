# Get Tango results, parse stretches, get resnums, create fastafiles for agadirwrapper, get difference total tango
import os
import sys
import glob
import subprocess
import GeneralUtilityMethods
from OptimizeProtein import yasara
import yaml

protein_chains = sys.argv[1].split('_')  # maybe get this through the constructor when you have made it into a class
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
__I_AC__ = 'Interaction_AnalyseComplex_'

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
# def _get_pdb_name_from_repaired_pdb_folder():
# repaired_pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
# return repaired_pdb.split('/')[-1].split('.')[0].split('_')[-1]
path_repaired_pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
pdb_name = path_repaired_pdb.split('/')[-1].split('.')[0].split('_')[-1] # removes the path and the repair prefix too

gatekeepers = ['R']  # for testing purposes only (to speed up solubis)
print pdb_name
globaltotal_Tango_all_chains_sum = 0
# def change_dir_to_pdb_Run_Solubis():
results_pdb_path = results_path + '/' + pdb_name
results_pdb_Runs_Solubis_path = results_pdb_path + '/Runs/Solubis'
os.chdir(results_pdb_Runs_Solubis_path)
individual_list_all_mutants = open('individual_list.txt', 'w')  # NOT SURE WHAT THIS LIST IS USED FOR
index_program = 0
############ RUN BUILDMODEL ON EACH PROTEIN CHAIN OF THIS PDB, 3 TIMES ###############
for protein_chain in protein_chains:
##### FOR EACH CHAIN, GET ITS TANGO WINDOW RESULTS , FROM  WHICH IT SHOULD GET EACH APR SEQUENCE AND START POSITION #####
    tangowindow_file = open(results_pdb_path + '/Agadir/' + pdb_name + '_' + protein_chain + '/PSX_tangowindow.out', 'r').readlines()
    for line in tangowindow_file[1:]:
        apr_sequence = line.split()[-4]
        apr_index = int(float(line.split()[-6]))
#### FOR EACH APR, MUTATE EVERY RESIDUE TO EVERY GATEKEEPER AND RUN BUILD MODEL TO CREATE REPAIRED MUTANT FILES ####
        for i, amino_acid in enumerate(apr_sequence):
            mutation_position = apr_index + i + 1
            for gatekeeper in gatekeepers:
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
                    subprocess.call('cp ' + results_pdb_path + '/../../SourceFiles/AgadirFiles/* ./Agadir/.',
                                    shell=True)
                if os.path.exists(results_pdb_path + '/Repair/RepairPDB_' + pdb_name + '.pdb'):
                    subprocess.call('cp ' + results_pdb_path + '/Repair/RepairPDB_' + pdb_name + '.pdb .', shell=True)
                    subprocess.call('cp ' + results_pdb_path + '/../../SourceFiles/FoldXFiles/* .', shell=True)
                else:
                    print 'Something is wrong'
                # def _run_FoldX_BuildModel(mutation_name):
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
                individual_list_for_this_mutant_only.write(mutant_name + ';\n')
                individual_list_for_this_mutant_only.close()

                grid_engine_job_name = __solubis_job_prefix__ + mutant_name
                no_queue = ''
                no_max_memory = ''
                no_cluster = ''
                using_runscript = True
                no_python_script = ''
                GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                           using_runscript, foldx_path, no_python_script)
                individual_list_all_mutants.write(mutant_name + ';\n')
                subprocess.call('qsub job.q', shell=True)

                os.chdir('./..')  # maybe pass the absolute path here

individual_list_all_mutants.close()  # NOT SURE WHAT THIS LIST IS USED FOR
GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__solubis_job_prefix__, 'all Solubis jobs to finish')

###################  RUN FOLDX ANALYSE COMPLEX on the outputs of BUILDMODEL ##################################
#  os.chdir(results_pdb_Runs_Solubis_path)  # this is the way to show the absolute path for glob glob otherwise it collates the whole path
all_solubis_mutant_folder_paths = sorted(glob.glob('./*'))
repaired_pdb_name = 'RepairPDB_' + pdb_name
repaired_pdb_name_1_ = repaired_pdb_name + '_1_'
wt_repaired_pdb_name_1_ = 'WT_' + repaired_pdb_name_1_
_0_1_2_pdbs = ['0.pdb,', '1.pdb,', '2.pdb,']

for solubis_mutant_folder_path in all_solubis_mutant_folder_paths:
    if os.path.isdir(solubis_mutant_folder_path):
        mutant_folder_name = solubis_mutant_folder_path.split('/')[-1]
        os.chdir(mutant_folder_name)
        subprocess.call('cp runscript.txt runscript_build.txt', shell=True)  # rename runscript to runscript_build.txt - WHY?
        subprocess.call('rm runscript.txt', shell=True)

        path_to_runscript = './'
        repaired_pdb_name_1_ = repaired_pdb_name + '_1_'
        wt_repaired_pdb_name_1_ = 'WT_' + repaired_pdb_name_1_
        pdbs_to_analyse = repaired_pdb_name_1_ + _0_1_2_pdbs[0] + \
                          repaired_pdb_name_1_ + _0_1_2_pdbs[1] + \
                          repaired_pdb_name_1_ + _0_1_2_pdbs[2] + \
                          wt_repaired_pdb_name_1_ + _0_1_2_pdbs[0] + \
                          wt_repaired_pdb_name_1_ + _0_1_2_pdbs[1] + \
                          wt_repaired_pdb_name_1_ + _0_1_2_pdbs[2]
        show_sequence_detail = False
        action = '<AnalyseComplex>#'
        print_networks = False
        calculate_stability = False
        GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, pdbs_to_analyse,
                                                           show_sequence_detail, action, print_networks,
                                                           calculate_stability)

## THIS PART OF CODE SEEMS TO SUGGEST THAT THE REPAIR ALGORITHM CAN CHANGE THE AMINO ACID SEQUENCE .. is that right ?? ###

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
            # _write_repaired_pdb_name_in_fasta_file():  HOWEVER, IT NAMES IT BY REPAIRED PDB BUT GETS SEQUENCE FROM REPAIRED PDB 1_0_ .. WHAT's going on here?
            repaired_pdb_fasta_file = open('Fasta/' + repaired_pdb_name + '_' + mol + '.fasta', 'w')
            repaired_pdb_fasta_file.write('>' + repaired_pdb_name + '_' + mol + '\n')
            repaired_pdb_fasta_file.write(fasta)
            repaired_pdb_fasta_file.close()

        grid_engine_job_name = __analyze_complex_job_prefix__ + mutant_folder_name
        no_queue = ''
        no_max_memory = ''
        no_cluster = ''
        using_runscript = True
        python_script_with_path = results_pdb_path + '/../../SourceFiles/Scripts/agadir.py'
        GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                   using_runscript, foldx_path, python_script_with_path)

        subprocess.call('qsub job.q', shell=True)
        os.chdir('./..')  # is this Runs/Solubis folder ?
        whatdir = os.getcwd()  # check here
    else:
        print solubis_mutant_folder_path

message_to_print = 'all AnalyseComplex jobs to finish'
GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__analyze_complex_job_prefix__, message_to_print)

#### WRITE A SUMMARY FILE WITH YOUR SOLUBIS RESULTS - not sure this summary text file actually contains anything ####
os.chdir(results_pdb_path)
solubis_summary_file = open('SummarySolubis.txt', 'w')
solubis_summary_file.write('Mutation\tMol\tddG\tdTANGO\tComplexSum\t')
has_interaction_AC_fxout_file_for_repaired_pdb = False
if os.path.isfile('./Repair/' + __I_AC__ + repaired_pdb_name + '.fxout'):  # perhaps use absolute path?
    has_interaction_AC_fxout_file_for_repaired_pdb = True
    interact_ac_pdb_fxout_file = open('./Repair/' + __I_AC__ + repaired_pdb_name + '.fxout').readlines()
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

whatiscurrentwd = os.getcwd()
## I'm not sure this works even in the orginal code, because the folder in results/pdbname. For the similar for and if on line 134 the cwd is runs/solubis
for solubis_mutant_folder_path in all_solubis_mutant_folder_paths:
    if os.path.isdir(solubis_mutant_folder_path):  # returns True if path is an existing directory (? subdirectory of this current dir ?)
        # mut = path.split('/')[-1]  # AH92R
        # mol = mut[1]  # H
        f = open(solubis_mutant_folder_path + '/Average_BuildModel_' + repaired_pdb_name + '.fxout', 'r').readlines()  # av_buildModel_repairedPdb_fxout_file
        ddG = f[9].split()[2]  # 13.3781

        all_agadir_output_files_for_this_mutant_chain_paths = glob.glob(solubis_mutant_folder_path + '/Agadir/*')
        TangoMut = 0
        TangoMut2 = 0  # check this works
        for agadir_output_files_for_this_mutant_chain_path in all_agadir_output_files_for_this_mutant_chain_paths:
            if os.path.isdir(agadir_output_files_for_this_mutant_chain_path): # isn't this redundant cos path_agad_list is obtained the line above by getting everything in the Agadir folder.
                # def get_globalTango_from_PSX_globaltotalout_file(path):
                f = open(agadir_output_files_for_this_mutant_chain_path + '/PSX_globaltotal.out', 'r').readlines()  # globalTango_file
                pieces = f[1].split()  # can delete if line 243 works
                TangoMol = float(pieces[2])  # can delete if line 243 works
                TangoMut = TangoMut + TangoMol  # can delete if line 243 works
                TangoMut2 += float(f[1].split()[2])  # check this works
                path_agad = agadir_output_files_for_this_mutant_chain_path  # ./Runs/Solubis/AH92R/Agadir/RepairPDB_Ab82b0sLigand_H
        print 'Tango GlobalTotal score for mutant protein = ' + str(TangoMut)
        print path_agad  # don't understand this. It goes through all the chains but only stores the last one to open the tango file on line below
        f = open(path_agad + '/PSX_globaltotal.out', 'r').readlines()  # opening the globalTango file again ?
        print f[1]  # which it then just prints and doesn't use ???
        # mol = mut[1]
        if has_interaction_AC_fxout_file_for_repaired_pdb:
            path_IAC_repaired_pdb_name_1_ = solubis_mutant_folder_path + '/' + __I_AC__ + repaired_pdb_name_1_
            path_IAC_wt_repaired_pdb_name_1_ = solubis_mutant_folder_path + '/' + __I_AC__ + wt_repaired_pdb_name_1_
            _0_1_2_fxout = ['0.fxout', '1.fxout', '2.fxout']
            f = open(path_IAC_repaired_pdb_name_1_ + _0_1_2_fxout[0], 'r').readlines()
            complex_Mut1 = []
            for line in f[9:]:
                complex_Mut1.append(float(line.split()[5]))
            f = open(path_IAC_repaired_pdb_name_1_ + _0_1_2_fxout[1], 'r').readlines()
            complex_Mut2 = []
            for line in f[9:]:
                complex_Mut2.append(float(line.split()[5]))
            f = open(path_IAC_repaired_pdb_name_1_ + _0_1_2_fxout[2], 'r').readlines()
            complex_Mut3 = []
            for line in f[9:]:
                complex_Mut3.append(float(line.split()[5]))
            f = open(path_IAC_wt_repaired_pdb_name_1_ + _0_1_2_fxout[0], 'r').readlines()
            complex_WT1 = []
            for line in f[9:]:
                complex_WT1.append(float(line.split()[5]))
            f = open(path_IAC_wt_repaired_pdb_name_1_ + _0_1_2_fxout[1], 'r').readlines()
            complex_WT2 = []
            for line in f[9:]:
                complex_WT2.append(float(line.split()[5]))
            f = open(path_IAC_wt_repaired_pdb_name_1_ + _0_1_2_fxout[2], 'r').readlines()
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
        mut = solubis_mutant_folder_path.split('/')[-1]
        mol = mut[1]
        mol2 = solubis_mutant_folder_path.split('/')[-1][1] # DOes this give same as above two lines?
        sumstring = mut + '\t' + mol + '\t' + ddG + '\t' + str(dTango) + '\t' + str(ComplexSum) + '\t' + Complex + '\n'
        solubis_summary_file.write(sumstring)
solubis_summary_file.close()
