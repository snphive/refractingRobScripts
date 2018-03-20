# Create individuals list that mutates every residue to every residue in the entire PDB.
import os,sys,glob,subprocess,time,yasara

Mols = sys.argv[1].split('_')
R_Path = sys.argv[2]
FoldX_Path = sys.argv[3]
Agadir_Path = sys.argv[4]

startingDir = os.getcwd()
AAlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

yasara.info.mode = 'txt'

aa_dict_1to3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}
aa_dict_3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','H1S':'H','H2S':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}

pdb = glob.glob('./Repair/RepairPDB*.pdb')[0]
print pdb
print os.getcwd()
name = pdb.split('/')[-1].split('.')[0].split('_')[-1]
rawpdb = pdb.split('/')[-1]
f = open('./Repair/SequenceDetail_RepairPDB_'+name+'.fxout').readlines()
os.chdir('Runs/Stabilize')
rescombos = []
for line in f[9:]:
	pieces = line.split('\t')
	if len(pieces) == 1:
		continue
	if len(pieces) == 0:
		continue
	if pieces[1] == 'GAP':
		continue
	if pieces[1] not in aa_dict_3to1.keys():
		continue
	res3 = pieces[1]
	res1 = aa_dict_3to1[res3]
	mol = pieces[2]
	if mol not in Mols:
		continue
	resnum = pieces[3]
	rescombo = [res1,mol,resnum]
	rescombo = "".join(rescombo)
	rescombos.append(rescombo)
indiv = open('individual_list.txt','w')
for x, aa in enumerate(rescombos):
	for newaa in AAlist:
		if newaa == aa[0]:
			continue
		templist = [aa,newaa]
		newaa = "".join(templist)
		indiv.write(newaa+';\n')
		if not os.path.exists(newaa):
			os.makedirs(newaa)
		os.chdir(newaa)
		if not os.path.exists('Fasta'):
			os.makedirs('Fasta')
		if not os.path.exists('Agadir'):
			os.makedirs('Agadir')
		if not os.path.exists('Agadir/Options.txt'):
			subprocess.call('cp '+startingDir+'/../../SourceFiles/AgadirFiles/* ./Agadir/.',shell=True)
		if os.path.exists(startingDir+'/Repair/RepairPDB_'+name+'.pdb'):
			if not os.path.exists('RepairPDB_'+name+'.pdb'):
				subprocess.call('cp '+startingDir+'/Repair/RepairPDB_'+name+'.pdb .',shell=True)
				subprocess.call('cp '+startingDir+'/../../SourceFiles/FoldXFiles/* .',shell=True)
		else:
			print 'Something is wrong'
		f = open('runscript.txt','w')
		f.write('<TITLE>FOLDX_runscript;\n')
		f.write('<JOBSTART>#;\n')
		f.write('<PDBS>RepairPDB_'+name+'.pdb;\n')
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
		h = open('individual_list.txt','w')
		h.write(newaa+';\n')
		h.close()
		g = open('./job.q','w')
		g.write('#!/bin/bash\n')
		g.write('#$ -N ST_'+newaa+'\n')
		g.write('#$ -V\n')
		g.write('#$ -cwd\n')
		g.write('source ~/.bash_profile\n')
		g.write(FoldX_Path+' -runfile runscript.txt\n')
		g.close()
		subprocess.call('qsub job.q',shell=True)
		os.chdir('./..')
indiv.close()
check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'ST_' in output_qstat:
	print 'Waiting for all Stabilize jobs to finish'
	time.sleep(10)
	check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
	output_qstat = check_qstat.stdout.read()

dirs = sorted(glob.glob('./*'))
name_agad = name
print 'Name agad:\t'+name_agad
name = 'RepairPDB_'+name
print 'Name:\t'+name
for path in dirs:
	if os.path.isdir(path):
		os.chdir(path.split('/')[-1])
		mutation = path.split('/')[-1]
		subprocess.call('cp runscript.txt runscript_build.txt',shell=True)
		subprocess.call('rm runscript.txt',shell=True)
		f = open('runscript.txt','w')
		f.write('<TITLE>FOLDX_runscript;\n')
		f.write('<JOBSTART>#;\n')
		f.write('<PDBS>'+name+'_1_0.pdb,'+name+'_1_1.pdb,'+name+'_1_2.pdb,WT_'+name+'_1_0.pdb,WT_'+name+'_1_1.pdb,WT_'+name+'_1_2.pdb,;\n')
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
		pdb = name+'_1_0.pdb'
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
				if line[21] == mol and resnum!=line[22:26].strip(' '):
					resnum = line[22:26].strip(' ')
					aa = line[17:20]
					fastalist.append(aa_dict_3to1[aa])
			fasta = "".join(fastalist)
			print name+'_'+mol
			print fasta
			f = open('Fasta/'+name+'_'+mol+'.fasta','w')
			f.write('>'+name+'_'+mol+'\n')
			f.write(fasta)
			f.close()
		g = open('./job.q','w')
		g.write('#!/bin/bash\n')
		g.write('#$ -N AC_'+mutation+'\n')
		g.write('#$ -V\n')
		g.write('#$ -cwd\n')
		g.write('source ~/.bash_profile\n')
		g.write(FoldX_Path+' -runfile runscript.txt\n')
		g.write('python '+startingDir+'/../../SourceFiles/Scripts/agadir.py\n')
		g.close()
		subprocess.call('qsub job.q',shell=True)
		os.chdir('./..')
	else:
		print path

check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'AC_' in output_qstat:
	print 'Waiting for all AnalyseComplex jobs to finish'
	time.sleep(10)
	check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
	output_qstat = check_qstat.stdout.read()
os.chdir(startingDir)
g = open('SummaryStabilize.txt','w')
g.write('Mutation\tMol\tddG\tdTANGO\tComplexSum\t')
analyseComplex = False
if os.path.isfile('./Repair/Interaction_AnalyseComplex_'+name+'.fxout'): 
	analyseComplex = True
	h = open('./Repair/Interaction_AnalyseComplex_'+name+'.fxout').readlines()
	Interactions = []
	for line in h[9:]:
		pieces = line.split('\t')
		mol1 = pieces[1]
		mol2 = pieces[2]
		molcomp = mol1+'_'+mol2
		Interactions.append(molcomp)
	interactionlist = []
	for i in Interactions:
		interactionlist.append('Complex_'+i)
	interactionstring = "\t".join(interactionlist)
	g.write(interactionstring+'\n')
else:
	g.write('\n')
dirs = sorted(glob.glob('./Runs/Stabilize/*'))

path_agad_list = glob.glob('./Agadir/*')
TangoWT = 0
for i in path_agad_list:
	if os.path.isdir(i):
		print i
		f = open(i+'/PSX_globaltotal.out','r').readlines()
		pieces = f[1].split()
		TangoMol = float(pieces[2])
		TangoWT = TangoWT + TangoMol
print 'TangoTotalWT = '+str(TangoWT)


for path in dirs:
	if os.path.isdir(path):
		mut = path.split('/')[-1]
		mol = mut[1]
		f = open(path+'/Average_BuildModel_'+name+'.fxout','r').readlines()
		ddG = f[9].split()[2]
		path_agad_list = glob.glob(path+'/Agadir/*')
		TangoMut = 0
		for i in path_agad_list:
			if os.path.isdir(i):
				f = open(i+'/PSX_globaltotal.out','r').readlines()
				pieces = f[1].split()
				TangoMol = float(pieces[2])
				TangoMut = TangoMut + TangoMol
				path_agad = i
		print 'TangoTotalMut = '+str(TangoMut)
		print path_agad
		f = open(path_agad+'/PSX_globaltotal.out','r').readlines()
		print f[1]
		mol = mut[1]
		if analyseComplex:
			f = open(path+'/Interaction_AnalyseComplex_'+name+'_1_0.fxout','r').readlines()
			complex_Mut1 = []
			for line in f[9:]:
					complex_Mut1.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_'+name+'_1_1.fxout','r').readlines()
			complex_Mut2 = []
			for line in f[9:]:
					complex_Mut2.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_'+name+'_1_2.fxout','r').readlines()
			complex_Mut3 = []
			for line in f[9:]:
					complex_Mut3.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name+'_1_0.fxout','r').readlines()
			complex_WT1 = []
			for line in f[9:]:
					complex_WT1.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name+'_1_1.fxout','r').readlines()
			complex_WT2 = []
			for line in f[9:]:
					complex_WT2.append(float(line.split()[5]))
			f = open(path+'/Interaction_AnalyseComplex_WT_'+name+'_1_2.fxout','r').readlines()
			complex_WT3 = []
			for line in f[9:]:
					complex_WT3.append(float(line.split()[5]))
			ComplexList=[]
			ComplexSum = 0
			for x,line in enumerate(f[9:]):
				complex_WT = float((complex_WT1[x]+complex_WT2[x]+complex_WT3[x])/3)
				complex_Mut = float((complex_Mut1[x]+complex_Mut2[x]+complex_Mut3[x])/3)
				Complex = complex_Mut-complex_WT
				ComplexSum = ComplexSum + Complex
				ComplexList.append(str(Complex))
			Complex = "\t".join(ComplexList)
		else:
			Complex = ""
			ComplexSum = 0
		dTango = float(TangoMut)-float(TangoWT)
		sumstring = mut+'\t'+mol+'\t'+ddG+'\t'+str(dTango)+'\t'+str(ComplexSum)+'\t'+Complex+'\n'
		g.write(sumstring)
g.close()




