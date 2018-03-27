import os
import glob
import subprocess
import time
import yasara

yasara.info.mode = 'txt'

Start = os.getcwd()
if not os.path.exists('Results'):
    os.makedirs('Results')

Scripts = Start + '/SourceFiles/Scripts'

# declare all global variables
__option_file__ = ''
__pdb_list__ = []
__command__ = ''
__charge__ = ''
__molecules__ = ''
__r_path__ = ''
__foldx_path__ = ''
__agadir_path__ = ''
__qsub_path__ = ''


def parse_option_file():

    global __pdb_list__
    global __command__
    global __molecules__
    global __r_path__
    global __foldx_path__
    global __agadir_path__
    global __qsub_path__
    global __option_file__
    __option_file__ = open('./OptProt_Options.txt', 'r').readlines()

    for line in __option_file__:
        if '#' in line:
            continue
        if 'PDBs:' in line:
            if line.split(':')[-1].strip(';\n') == "All":
                pdbpaths = glob.glob('./PDBs/*.pdb')
                for pdbpath in pdbpaths:
                    pdb_temp = pdbpath.split('/')[-1]
                    __pdb_list__.append(pdb_temp)
            pdb_string = line.split(':')[-1].strip(';\n')
            if ',' in pdb_string:
                pdbs_temp = pdb_string.split(',')
                for pdb_temp in pdbs_temp:
                    __pdb_list__.append(pdb_temp)
            else:
                __pdb_list__.append(pdb_string)
        if 'Command:' in line:
            __command__ = line.split(':')[-1].strip(';\n')
        if 'Charge:' in line:
            __charge__ = line.split(':')[-1].strip(';\n')
        if 'Mols:' in line:
            __molecules__ = []
            MolString = line.split(':')[-1].strip(';\n')
            if ',' in MolString:
                Mols_temp = MolString.split(',')
                for Mols_temp in Mols_temp:
                    __molecules__.append(Mols_temp)
                __molecules__ = "_".join(__molecules__)
            else:
                __molecules__ = MolString
        if 'R_Path:' in line:
            __r_path__ = line.split(':')[-1].strip(';\n')
        if 'FoldX_Path:' in line:
            __foldx_path__ = line.split(':')[-1].strip(';\n')
        if 'Agadir_Path:' in line:
            __agadir_path__ = line.split(':')[-1].strip(';\n')
        if 'Qsub_Path' in line:
            __qsub_path__ = line.split(':')[-1].strip(';\n')
    print 'PDBs to analyse:\t\t' + ",\t".join(__pdb_list__)
    print 'Command to be executed:\t\t' + __command__
    print 'Molecules to be considered:\t' + __molecules__
    print 'Absolute path to R:\t\t' + __r_path__
    print 'Absolute path to FoldX:\t\t' + __foldx_path__
    print 'Absolute path to TANGO:\t\t' + __agadir_path__
    print 'Absolute path to Qsub:\t\t' + __qsub_path__


parse_option_file(__option_file__)

for PDB in __pdb_list__:
    name = PDB.split('.')[0]

    if not os.path.exists('Results/' + name):
        os.makedirs('Results/' + name)
    os.chdir('Results/' + name)
    if not os.path.exists('Repair'):
        os.makedirs('Repair')
    if not os.path.exists('Agadir'):
        os.makedirs('Agadir')
    if not os.path.exists('Runs'):
        os.makedirs('Runs')
    if not os.path.exists('PDBs'):
        os.makedirs('PDBs')
    if not os.path.exists('Fasta'):
        os.makedirs('Fasta')
    yasara.run('DelObj all')
    yasara.run('LoadPDB ' + Start + '/PDBs/' + PDB)
    yasara.run('DelRes !Protein')
    tempMols = yasara.run('ListMol All,Format=MOLNAME')

    for mol in tempMols:
        yasara.run('RenumberRes all and Mol ' + mol + ',First=1')

    yasara.run('SavePDB 1,' + Start + '/Results/' + name + '/PDBs/' + PDB + ',Format=PDB,Transform=Yes')
    subprocess.call('cp ' + Start + '/PDBs/' + PDB + ' ./PDBs/.', shell=True)
    subprocess.call('cp ' + Start + '/PDBs/' + PDB + ' ./Repair/.', shell=True)
    subprocess.call('cp ' + Start + '/SourceFiles/FoldXFiles/* ./Repair/.', shell=True)
    subprocess.call('cp ' + Start + '/SourceFiles/AgadirFiles/* ./Agadir/.', shell=True)
    print 'pdb2fasta.py'
    subprocess.call('python ' + Scripts + '/pdb2fasta.py', shell=True)
    print 'agadir.py'
    subprocess.call('python ' + Scripts + '/agadir.py', shell=True)
    print 'repair.py'

    g = open('./job.q', 'w')
    g.write('#!/bin/bash\n')
    g.write('#$ -N RPjob_' + name + '\n')
    g.write('#$ -V\n')
    g.write('#$ -cwd\n')
    g.write('source ~/.bash_profile\n')
    g.write('python ' + Scripts + '/repair.py ' + __qsub_path__ + '\n')
    g.close()
    subprocess.call(__qsub_path__ + 'qsub job.q', shell=True)
    os.chdir(Start)

check_qstat = subprocess.Popen(__qsub_path__ + 'qstat', stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()

while 'RPjob_' in output_qstat:
    print 'Waiting for PDBs to be repaired'
    time.sleep(10)
    check_qstat = subprocess.Popen(__qsub_path__ + 'qstat', stdout=subprocess.PIPE)
    output_qstat = check_qstat.stdout.read()

for PDB in __pdb_list__:
    name = PDB.split('.')[0]

    if not os.path.exists('Results/' + name):
        os.makedirs('Results/' + name)
    os.chdir('Results/' + name)
    if not os.path.exists('Repair'):
        os.makedirs('Repair')
    if not os.path.exists('Agadir'):
        os.makedirs('Agadir')
    if not os.path.exists('Runs'):
        os.makedirs('Runs')
    if not os.path.exists('PDBs'):
        os.makedirs('PDBs')
    if not os.path.exists('Fasta'):
        os.makedirs('Fasta')

    print 'stretchplot.py'
    subprocess.call(
        'python ' + Scripts + '/stretchplot.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__,
        shell=True)
    print 'stretchplot.R'
    subprocess.call('R < ' + Scripts + '/stretchplot.R --no-save', shell=True)

    if not os.path.exists('Runs/' + __command__):
        if __command__ != 'Stretchplot' and __command__ != 'All':
            os.makedirs('Runs/' + __command__)
        if __command__ == 'All':
            os.makedirs('Runs/Stabilize')
            os.makedirs('Runs/Solubis')
            os.makedirs('Runs/SolubisMild')
            os.makedirs('Runs/CysScan')
    if __command__ == 'Stabilize' or __command__ == 'All':
        print 'stabilize.py'
        subprocess.call(
            'python ' + Scripts + '/stabilize.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__,
            shell=True)
    if __command__ == 'Solubis' or __command__ == 'All':
        print 'solubis.py'
        subprocess.call(
            'python ' + Scripts + '/solubis.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__,
            shell=True)
    if __command__ == 'SolubisMild' or __command__ == 'All':
        print 'solubis_mild.py'
        subprocess.call(
            'python ' + Scripts + '/solubis_mild.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__,
            shell=True)
    if __command__ == 'CysScan' or __command__ == 'All':
        print 'CysScan.py'
        subprocess.call(
            'python ' + Scripts + '/CysScan.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__,
            shell=True)
    if __command__ == 'SolubisDouble':
        print 'solubis_double.py'
        subprocess.call(
            'python ' + Scripts + '/solubis_double.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__,
            shell=True)
    if __command__ == 'Supercharge':
        print 'supercharge.py'
        subprocess.call(
            'python ' + Scripts + '/supercharge.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__ + ' ' + __charge__,
            shell=True)
    if __command__ == 'DelPos':
        print 'DelPos.py'
        subprocess.call(
            'python ' + Scripts + '/DelPos.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__ + ' ' + __charge__,
            shell=True)
    if __command__ == 'Indiv':
        print 'Indiv.py'
        subprocess.call(
            'python ' + Scripts + '/Indiv.py ' + __molecules__ + ' ' + __r_path__ + ' ' + __foldx_path__ + ' ' + __agadir_path__ + ' ' + __charge__,
            shell=True)
    subprocess.call('R < ' + Scripts + '/massplots.R --no-save', shell=True)
    # subprocess.call('python '+Scripts+'/clean.py',shell=True)
    os.chdir(Start)
