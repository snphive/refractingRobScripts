# Get Tango results, parse stretches, get resnums, create fastafiles for agadirwrapper, get difference total tango
import os
import sys
import glob
import subprocess
import time
from OptimizeProtein import yasara
import yaml

Mols = sys.argv[1].split('_')
# R_Path = sys.argv[2] # this is not used so why is it here?
# FoldX_Path = sys.argv[3]
# Agadir_Path = sys.argv[4] # this is not used so why is it here?

yasara.info.mode = 'txt'

r_path = ''
foldx_path = ''
agadir_path = ''
qsub_path = ''  # don't think this is needed so may remove it.. (from repair.py too)
starting_directory = ''
results_directory = ''
aa_dict_1to3 = {}
aa_dict_3to1 = {}
gatekeepers = []

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

# aa_dict_1to3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
# 				'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
# 				'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}
# aa_dict_3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'H1S': 'H',
# 				'H2S': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
# 				'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


# starting_directory = os.getcwd()
pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
# gatekeepers = ['R', 'P', 'K', 'E', 'D']
gatekeepers = ['R']  # commented out full list of gatekeepers and using just 1 residue to speed up my test runs
pdb_name = pdb.split('/')[-1].split('.')[0].split('_')[-1]
starting_directory = results_directory + '/' + pdb_name
print pdb_name
agadirs = []
totals = []
total_pdb = 0
os.chdir('Runs/Solubis')
indiv = open('individual_list.txt', 'w')
index_program = 0
for mol in Mols:
    f = open(starting_directory + '/Agadir/' + pdb_name + '_' + mol + '/PSX_globaltotal.out', 'r').readlines()
    header_for_PSX_globaltotal_table = f[0].split()
    old_fasta_lines = open(starting_directory + '/Fasta/' + pdb_name + '_' + mol + '.fasta', 'r').readlines()
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
    f = open(starting_directory + '/Agadir/' + pdb_name + '_' + mol + '/PSX_tangowindow.out', 'r').readlines()
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
            for gate in gatekeepers:
                new_stretch = list(stretch)
                new_stretch[x] = gate
                new_stretch = "".join(new_stretch)
                new_fasta = old_fasta.replace(stretch, new_stretch)
                mutation = aa + mol + index_mutation + gate
                if not os.path.exists(mutation):
                    os.makedirs(mutation)
                os.chdir(mutation)
                if not os.path.exists('Fasta'):
                    os.makedirs('Fasta')
                if not os.path.exists('Agadir'):
                    os.makedirs('Agadir')
                if not os.path.exists('Agadir/Options.txt'):
                    subprocess.call('cp ' + starting_directory + '/../../SourceFiles/AgadirFiles/* ./Agadir/.',
                                    shell=True)
                if os.path.exists(starting_directory + '/Repair/RepairPDB_' + pdb_name + '.pdb'):
                    subprocess.call('cp ' + starting_directory + '/Repair/RepairPDB_' + pdb_name + '.pdb .', shell=True)
                    subprocess.call('cp ' + starting_directory + '/../../SourceFiles/FoldXFiles/* .', shell=True)
                else:
                    print 'Something is wrong'
                f = open('runscript.txt', 'w')
                f.write('<TITLE>FOLDX_runscript;\n')
                f.write('<JOBSTART>#;\n')
                f.write('<PDBS>RepairPDB_' + pdb_name + '.pdb;\n')
                f.write('<BATCH>#;\n')
                f.write('<COMMANDS>FOLDX_commandfile;\n')
                f.write('<BuildModel>#,individual_list.txt;\n')
                f.write('<END>#;\n')
                f.write('<OPTIONS>FOLDX_optionfile;\n')
                f.write('<Temperature>298;\n')
                f.write('<IonStrength>0.05;\n')
                f.write('<ph>7;\n')
                f.write('<moveNeighbours>true;\n')
                f.write('<VdWDesign>2;\n')
                f.write('<numberOfRuns>3;\n')
                f.write('<OutPDB>#;\n')
                f.write('<END>#;\n')
                f.write('<JOBEND>#;\n')
                f.write('<ENDFILE>#;\n')
                f.close()
                h = open('individual_list.txt', 'w')
                h.write(mutation + ';\n')
                h.close()
                g = open('./job.q', 'w')
                g.write('#!/bin/bash\n')
                g.write('#$ -N SB_' + mutation + '\n')
                g.write('#$ -V\n')
                g.write('#$ -cwd\n')
                g.write('source ~/.bash_profile\n')
                g.write(foldx_path + ' -runfile runscript.txt\n')
                g.close()
                indiv.write(mutation + ';\n')
                subprocess.call('qsub job.q', shell=True)
                os.chdir('./..')

indiv.close()
check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'SB_' in output_qstat:
    print 'Waiting for all Solubis jobs to finish'
    time.sleep(10)
    check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
    output_qstat = check_qstat.stdout.read()
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
        f = open('runscript.txt', 'w')
        f.write('<TITLE>FOLDX_runscript;\n')
        f.write('<JOBSTART>#;\n')
        f.write(
            '<PDBS>' + pdb_name + '_1_0.pdb,' + pdb_name + '_1_1.pdb,' + pdb_name + '_1_2.pdb,WT_' + pdb_name + '_1_0.pdb,WT_'
            + pdb_name + '_1_1.pdb,WT_' + pdb_name + '_1_2.pdb,;\n')
        f.write('<BATCH>#;\n')
        f.write('<COMMANDS>FOLDX_commandfile;\n')
        f.write('<AnalyseComplex>#;\n')
        f.write('<END>#;\n')
        f.write('<OPTIONS>FOLDX_optionfile;\n')
        f.write('<Temperature>298;\n')
        f.write('<IonStrength>0.05;\n')
        f.write('<ph>7;\n')
        f.write('<moveNeighbours>true;\n')
        f.write('<VdWDesign>2;\n')
        f.write('<numberOfRuns>3;\n')
        f.write('<OutPDB>#;\n')
        f.write('<END>#;\n')
        f.write('<JOBEND>#;\n')
        f.write('<ENDFILE>#;\n')
        f.close()
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
        g = open('./job.q', 'w')
        g.write('#!/bin/bash\n')
        g.write('#$ -N AC_' + mutation + '\n')
        g.write('#$ -V\n')
        g.write('#$ -cwd\n')
        g.write('source ~/.bash_profile\n')
        g.write(foldx_path + ' -runfile runscript.txt\n')
        g.write('python ' + starting_directory + '/../../SourceFiles/Scripts/agadir.py\n')
        g.close()
        subprocess.call('qsub job.q', shell=True)
        os.chdir('./..')
    else:
        print path

check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'AC_' in output_qstat:
    print 'Waiting for all AnalyseComplex jobs to finish'
    time.sleep(10)
    check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
    output_qstat = check_qstat.stdout.read()
os.chdir(starting_directory)
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
