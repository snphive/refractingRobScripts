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
qsub_path = ''  # don't think this is needed so may remove it.. (from repair.py too)
results_directory = ''
results_pdb_directory = ''
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
        qsub_path = paths_and_dictionaries['ROOT']['Qsub_path']
        results_directory = paths_and_dictionaries['ROOT']['results_path']
        aa_dict_1to3 = paths_and_dictionaries['ROOT']['aa_dict_1to3']
        aa_dict_3to1 = paths_and_dictionaries['ROOT']['aa_dict_3to1']
        gatekeepers = paths_and_dictionaries['ROOT']['gatekeepers']

    except yaml.YAMLError as exc:
        print(exc)

# starting_directory = os.getcwd()
pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
gatekeepers = ['R']  # for testing purposes only (to speed up solubis)
pdb_name = pdb.split('/')[-1].split('.')[0].split('_')[-1]
results_pdb_directory = results_directory + '/' + pdb_name
results_pdb_Runs_Solubis_directory = results_pdb_directory + '/Runs/Solubis'
print pdb_name
agadirs = []
totals = []
total_pdb = 0

os.chdir(results_pdb_Runs_Solubis_directory)

indiv = open('individual_list.txt', 'w')
index_program = 0
for protein_chain in protein_chains:
    f = open(results_pdb_directory + '/Agadir/' + pdb_name + '_' + protein_chain + '/PSX_globaltotal.out', 'r').readlines()
    header_for_PSX_globaltotal_table = f[0].split()
    old_fasta_lines = open(results_pdb_directory + '/Fasta/' + pdb_name + '_' + protein_chain + '.fasta', 'r').readlines()
    old_fasta = []
    if len(old_fasta_lines) == 2:
        for a in old_fasta_lines[1]:
            old_fasta.append(a)
    else:
        for old_fasta_line in old_fasta_lines[1:]:
            for a in old_fasta_line:
                old_fasta.append(a)
    old_fasta = "".join(old_fasta)
    for x, program in enumerate(header_for_PSX_globaltotal_table):
        if 'TANGO' == program:
            index_program = x
    total = f[1].split()[index_program]
    totals.append(total)
    total_pdb = total_pdb + float(total)
    # get stretches
    f = open(results_pdb_directory + '/Agadir/' + pdb_name + '_' + protein_chain + '/PSX_tangowindow.out', 'r').readlines()
    stretches = []
    for line in f[1:]:
        new_stretches = []
        pieces = line.split()
        score = pieces[-2]
        length = pieces[-1]
        stretch_with_gk = pieces[2:5]
        stretch_with_gk = "".join(stretch_with_gk)
        stretch = pieces[-4]
        index_stretch = old_fasta.find(stretch)
        for x, aa in enumerate(stretch):
            index_mutation = str(index_stretch + x + 1)
            for gatekeeper in gatekeepers:
                new_stretch = list(stretch)
                new_stretch[x] = gatekeeper
                new_stretch = "".join(new_stretch)
                new_fasta = old_fasta.replace(stretch, new_stretch)
                mutation = aa + protein_chain + index_mutation + gatekeeper
                if not os.path.exists(mutation):
                    os.makedirs(mutation)
                os.chdir(mutation)
                if not os.path.exists('Fasta'):
                   os.makedirs('Fasta')
                if not os.path.exists('Agadir'):
                    os.makedirs('Agadir')
                if not os.path.exists('Agadir/Options.txt'):
                    subprocess.call('cp ' + results_pdb_directory + '/../../SourceFiles/AgadirFiles/* ./Agadir/.',
                                shell=True)
                if os.path.exists(results_pdb_directory + '/Repair/RepairPDB_' + pdb_name + '.pdb'):
                    subprocess.call('cp ' + results_pdb_directory + '/Repair/RepairPDB_' + pdb_name + '.pdb .', shell=True)
                    subprocess.call('cp ' + results_pdb_directory + '/../../SourceFiles/FoldXFiles/* .', shell=True)
                else:
                    print 'Something is wrong'
                repaired_pdbs = 'RepairPDB_' + pdb_name + '.pdb'
                runscript_foldx_command = '<BuildModel>#,individual_list.txt'
                path_to_runscript = './'
                should_print_networks = False
                GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, repaired_pdbs,
                                                                   runscript_foldx_command, should_print_networks)

                h = open('individual_list.txt', 'w')
                h.write(mutation + ';\n')
                h.close()

                grid_engine_job_name = __solubis_job_prefix__ + mutation
                no_queue = ''
                no_max_memory = ''
                no_cluster = ''
                no_python_script = ''
                GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                           foldx_path, no_python_script)

                indiv.write(mutation + ';\n')
                subprocess.call('qsub job.q', shell=True)
                os.chdir('./..')

indiv.close()
GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__solubis_job_prefix__, 'all Solubis jobs to finish')

# ListOfSolubisMutationResultsFolders
dirs = sorted(glob.glob('./*'))
name_agad = pdb_name
print 'Name agad:\t' + name_agad
pdb_name = 'RepairPDB_' + pdb_name
print 'Name:\t' + pdb_name
# SolubisMutationResultsFolder
for path in dirs:
    if os.path.isdir(path):
        os.chdir(path.split('/')[-1])
        mutation = path.split('/')[-1]
        subprocess.call('cp runscript.txt runscript_build.txt', shell=True)
        subprocess.call('rm runscript.txt', shell=True)

        pdbs_to_analyse = pdb_name + '_1_0.pdb,' + pdb_name + '_1_1.pdb,' + pdb_name + '_1_2.pdb,WT_' + pdb_name + \
                          '_1_0.pdb,WT_' + pdb_name + '_1_1.pdb,WT_' + pdb_name + '_1_2.pdb,'
        runscript_foldx_command = '<AnalyseComplex>#'
        path_to_runscript = './'
        should_print_networks = False
        GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, pdbs_to_analyse,
                                                           runscript_foldx_command, should_print_networks)

        pdb = pdb_name + '_1_0.pdb'
        print os.getcwd()
        f = open(pdb).readlines()
        atomlines = []
        mols = []
        for line in f:
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
            print pdb_name + '_' + mol
            print fasta
            f = open('Fasta/' + pdb_name + '_' + mol + '.fasta', 'w')
            f.write('>' + pdb_name + '_' + mol + '\n')
            f.write(fasta)
            f.close()

        grid_engine_job_name = __analyze_complex_job_prefix__ + mutation
        no_queue = ''
        no_max_memory = ''
        no_cluster = ''
        python_script_with_path = results_pdb_directory + '/../../SourceFiles/Scripts/agadir.py'
        GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, no_queue, no_max_memory, no_cluster,
                                                   foldx_path, python_script_with_path)

        subprocess.call('qsub job.q', shell=True)
        os.chdir('./..')
    else:
        print path

message_to_print = 'all AnalyseComplex jobs to finish'
GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__analyze_complex_job_prefix__, message_to_print)

os.chdir(results_pdb_directory)
g = open('SummarySolubis.txt', 'w')
g.write('Mutation\tMol\tddG\tdTANGO\tComplexSum\t')
analyseComplex = False
if os.path.isfile('./Repair/Interaction_AnalyseComplex_' + pdb_name + '.fxout'):
    analyseComplex = True
    h = open('./Repair/Interaction_AnalyseComplex_' + pdb_name + '.fxout').readlines()
    Interactions = []
    for line in h[9:]:
        pieces = line.split('\t')
        mol1 = pieces[1]
        mol2 = pieces[2]
        molcomp = mol1 + '_' + mol2
        Interactions.append(molcomp)
    interactionlist = []
    for i in Interactions:
        interactionlist.append('Complex_' + i)
    interactionstring = "\t".join(interactionlist)
    g.write(interactionstring + '\n')
else:
    g.write('\n')
dirs = sorted(glob.glob('./Runs/Solubis/*'))

path_agad_list = glob.glob('./Agadir/*')
TangoWT = 0
for i in path_agad_list:
    if os.path.isdir(i):
        print i
        f = open(i + '/PSX_globaltotal.out', 'r').readlines()
        pieces = f[1].split()
        TangoMol = float(pieces[2])
        TangoWT = TangoWT + TangoMol
print 'TangoTotalWT = ' + str(TangoWT)

for path in dirs:
    if os.path.isdir(path):
        mut = path.split('/')[-1]
        mol = mut[1]
        f = open(path + '/Average_BuildModel_' + pdb_name + '.fxout', 'r').readlines()
        ddG = f[9].split()[2]
        path_agad_list = glob.glob(path + '/Agadir/*')
        TangoMut = 0
        for i in path_agad_list:
            if os.path.isdir(i):
                f = open(i + '/PSX_globaltotal.out', 'r').readlines()
                pieces = f[1].split()
                TangoMol = float(pieces[2])
                TangoMut = TangoMut + TangoMol
                path_agad = i
        print 'TangoTotalMut = ' + str(TangoMut)
        print path_agad
        f = open(path_agad + '/PSX_globaltotal.out', 'r').readlines()
        print f[1]
        mol = mut[1]
        if analyseComplex:
            f = open(path + '/Interaction_AnalyseComplex_' + pdb_name + '_1_0.fxout', 'r').readlines()
            complex_Mut1 = []
            for line in f[9:]:
                complex_Mut1.append(float(line.split()[5]))
            f = open(path + '/Interaction_AnalyseComplex_' + pdb_name + '_1_1.fxout', 'r').readlines()
            complex_Mut2 = []
            for line in f[9:]:
                complex_Mut2.append(float(line.split()[5]))
            f = open(path + '/Interaction_AnalyseComplex_' + pdb_name + '_1_2.fxout', 'r').readlines()
            complex_Mut3 = []
            for line in f[9:]:
                complex_Mut3.append(float(line.split()[5]))
            f = open(path + '/Interaction_AnalyseComplex_WT_' + pdb_name + '_1_0.fxout', 'r').readlines()
            complex_WT1 = []
            for line in f[9:]:
                complex_WT1.append(float(line.split()[5]))
            f = open(path + '/Interaction_AnalyseComplex_WT_' + pdb_name + '_1_1.fxout', 'r').readlines()
            complex_WT2 = []
            for line in f[9:]:
                complex_WT2.append(float(line.split()[5]))
            f = open(path + '/Interaction_AnalyseComplex_WT_' + pdb_name + '_1_2.fxout', 'r').readlines()
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
        sumstring = mut + '\t' + mol + '\t' + ddG + '\t' + str(dTango) + '\t' + str(ComplexSum) + '\t' + Complex + '\n'
        g.write(sumstring)
g.close()
