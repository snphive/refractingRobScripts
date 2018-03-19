import os,sys,glob,subprocess,time,yasara

yasara.info.mode = 'txt'

Start = os.getcwd()
if not os.path.exists('Results'):
	os.makedirs('Results')

Scripts = Start+'/SourceFiles/Scripts'

OptionFile = open('./OptProt_Options.txt','r').readlines()

for line in OptionFile:
	if '#' in line:
		continue
	if 'PDBs:' in line:
		PDBs = []
		if line.split(':')[-1].strip(';\n') == "All":
			pdbpaths = glob.glob('./PDBs/*.pdb')
			for pdbpath in pdbpaths:
				pdb_temp = pdbpath.split('/')[-1]
				PDBs.append(pdb_temp)
		PDBString = line.split(':')[-1].strip(';\n')
		if ',' in PDBString:
			pdbs_temp = PDBString.split(',')
			for pdb_temp in pdbs_temp:
				PDBs.append(pdb_temp)
		else:
			PDBs.append(PDBString)
	if 'Command:' in line:
		Command = line.split(':')[-1].strip(';\n')
	if 'Charge:' in line:
		Charge = line.split(':')[-1].strip(';\n')
	if 'Mols:' in line:
		Mols = []
		MolString = line.split(':')[-1].strip(';\n')
		if ',' in MolString:
			Mols_temp = MolString.split(',')
			for Mols_temp in Mols_temp:
				Mols.append(Mols_temp)
			Mols = "_".join(Mols)
		else:
			Mols = MolString
	if 'R_Path:' in line:
		R_Path = line.split(':')[-1].strip(';\n')
	if 'FoldX_Path:' in line:
		FoldX_Path = line.split(':')[-1].strip(';\n')
	if 'Agadir_Path:' in line:
		Agadir_Path = line.split(':')[-1].strip(';\n')
	if 'Qsub_Path' in line:
		Qsub_Path = line.split(':')[-1].strip(';\n')
print 'PDBs to analyse:\t\t'+",\t".join(PDBs)
print 'Command to be executed:\t\t'+Command
print 'Molecules to be considered:\t'+Mols
print 'Absolute path to R:\t\t'+R_Path
print 'Absolute path to FoldX:\t\t'+FoldX_Path
print 'Absolute path to TANGO:\t\t'+Agadir_Path
print 'Absolute path to Qsub:\t\t'+Qsub_Path

for PDB in PDBs:
	name = PDB.split('.')[0]
	if not os.path.exists('Results/'+name):
		os.makedirs('Results/'+name)
	os.chdir('Results/'+name)
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
	yasara.run('LoadPDB '+Start+'/PDBs/'+PDB)
	yasara.run('DelRes !Protein')
	tempMols = yasara.run('ListMol All,Format=MOLNAME')
	for mol in tempMols:
		yasara.run('RenumberRes all and Mol '+mol+',First=1')
	yasara.run('SavePDB 1,'+Start+'/Results/'+name+'/PDBs/'+PDB+',Format=PDB,Transform=Yes')
	subprocess.call('cp '+Start+'/PDBs/'+PDB+' ./PDBs/.',shell=True)
	subprocess.call('cp '+Start+'/PDBs/'+PDB+' ./Repair/.',shell=True)
	subprocess.call('cp '+Start+'/SourceFiles/FoldXFiles/* ./Repair/.',shell=True)
	subprocess.call('cp '+Start+'/SourceFiles/AgadirFiles/* ./Agadir/.',shell=True)
	print 'pdb2fasta.py'
	subprocess.call('python '+Scripts+'/pdb2fasta.py',shell=True)
	print 'agadir.py'
	subprocess.call('python '+Scripts+'/agadir.py',shell=True)
	print 'repair.py'
	g = open('./job.q','w')
	g.write('#!/bin/bash\n')
	g.write('#$ -N RPjob_'+name+'\n')
	g.write('#$ -V\n')
	g.write('#$ -cwd\n')
	g.write('source ~/.bash_profile\n')
	g.write('python '+Scripts+'/repair.py '+Qsub_Path+'\n')
	g.close()
	subprocess.call(Qsub_Path+'qsub job.q',shell=True)
	os.chdir(Start)

check_qstat = subprocess.Popen(Qsub_Path+'qstat',stdout=subprocess.PIPE)
output_qstat = check_qstat.stdout.read()
while 'RPjob_' in output_qstat:
	print 'Waiting for PDBs to be repaired'
	time.sleep(10)
	check_qstat = subprocess.Popen(Qsub_Path+'qstat',stdout=subprocess.PIPE)
	output_qstat = check_qstat.stdout.read()

for PDB in PDBs:
	name = PDB.split('.')[0]
	if not os.path.exists('Results/'+name):
		os.makedirs('Results/'+name)
	os.chdir('Results/'+name)
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
	subprocess.call('python '+Scripts+'/stretchplot.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path,shell=True)
	print 'stretchplot.R'
	subprocess.call('R < '+Scripts+'/stretchplot.R --no-save',shell=True)
	if not os.path.exists('Runs/'+Command):
		if Command != 'Stretchplot' and Command != 'All':
			os.makedirs('Runs/'+Command)
		if Command == 'All':
			os.makedirs('Runs/Stabilize')
			os.makedirs('Runs/Solubis')
			os.makedirs('Runs/SolubisMild')
			os.makedirs('Runs/CysScan')
	if Command == 'Stabilize' or Command == 'All':
		print 'stabilize.py'
		subprocess.call('python '+Scripts+'/stabilize.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path,shell=True)
	if Command == 'Solubis' or Command == 'All':
		print 'solubis.py'
		subprocess.call('python '+Scripts+'/solubis.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path,shell=True)
	if Command == 'SolubisMild' or Command == 'All':
		print 'solubis_mild.py'
		subprocess.call('python '+Scripts+'/solubis_mild.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path,shell=True)
	if Command == 'CysScan' or Command == 'All':
		print 'CysScan.py'
		subprocess.call('python '+Scripts+'/CysScan.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path,shell=True)
	if Command == 'SolubisDouble':
		print 'solubis_double.py'
		subprocess.call('python '+Scripts+'/solubis_double.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path,shell=True)
	if Command == 'Supercharge':
		print 'supercharge.py'
		subprocess.call('python '+Scripts+'/supercharge.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path+' '+Charge,shell=True)
	if Command == 'DelPos':
		print 'DelPos.py'
		subprocess.call('python '+Scripts+'/DelPos.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path+' '+Charge,shell=True)
	if Command == 'Indiv':
		print 'Indiv.py'
		subprocess.call('python '+Scripts+'/Indiv.py '+Mols+' '+R_Path+' '+FoldX_Path+' '+Agadir_Path+' '+Charge,shell=True)
	subprocess.call('R < '+Scripts+'/massplots.R --no-save',shell=True)
	#subprocess.call('python '+Scripts+'/clean.py',shell=True)
	os.chdir(Start)
