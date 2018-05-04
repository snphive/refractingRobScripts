import os,sys,glob,subprocess,time
import GeneralUtilityMethods
startingDir = os.getcwd()
print 'startDir in Repair.py' + startingDir
targets = sorted(glob.glob('./PDBs/*.pdb'))
foldx = 'FoldX'
results = 'Results'

Qsub_Path = sys.argv[1]
for pdb in targets:
	name = pdb.split('/')[-1].split('.')[0]
	rawpdb = pdb.split('/')[-1]
	os.chdir('Repair')
	if not os.path.exists('RepairPDB_'+rawpdb):
		path_to_runscript = startingDir + '/Repair/'
		pdbs = rawpdb
		action = '<RepairPDB>#'
		should_print_networks = True
		GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, pdbs, action, should_print_networks)
		g = open('./job.q','w')
		g.write('#!/bin/bash\n')
		g.write('#$ -N RP_'+name+'\n')
		g.write('#$ -V\n')
		g.write('#$ -q all.q\n')
		g.write('#$ -l h_vmem=3G\n')
		g.write('#$ -l hostname=hodor1.vib\n')
		g.write('#$ -cwd\n')
		g.write('source ~/.bash_profile\n')
		g.write('/switchlab/group/tools/FoldX_2015/FoldX -runfile runscript.txt\n')
		g.close()
		subprocess.call(Qsub_Path+'qsub job.q',shell=True)
		check_qstat = subprocess.Popen(Qsub_Path+'qstat',stdout=subprocess.PIPE)
		output_qstat = check_qstat.stdout.read()
		while 'RP_' in output_qstat:
			print 'Waiting for RepairPDB to finish'
			time.sleep(10)
			check_qstat = subprocess.Popen(Qsub_Path+'qstat',stdout=subprocess.PIPE)
			output_qstat = check_qstat.stdout.read()
		subprocess.call('cp '+startingDir+'/Repair/runscript.txt runscript_repair.txt',shell=True)
		f = open(startingDir+'/Repair/runscript.txt','w')
		f.write('<TITLE>FOLDX_runscript;\n')
		f.write('<JOBSTART>#;\n')
		f.write('<PDBS>RepairPDB_'+rawpdb+';\n')
		f.write('<BATCH>#;\n')
		f.write('<COMMANDS>FOLDX_commandfile;\n')
		f.write('<SequenceDetail>#;\n')
		f.write('<AnalyseComplex>#;\n')
		f.write('<PrintNetworks>#;\n')
		f.write('<Stability>#;\n')
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
		subprocess.call(Qsub_Path+'qsub job.q',shell=True)
		check_qstat = subprocess.Popen(Qsub_Path+'qstat',stdout=subprocess.PIPE)
		output_qstat = check_qstat.stdout.read()
		while 'RP_' in output_qstat:
			print 'Waiting for AnalyseComplex to finish'
			time.sleep(5)
			check_qstat = subprocess.Popen('qstat',stdout=subprocess.PIPE)
			output_qstat = check_qstat.stdout.read()

		subprocess.call('cp '+startingDir+'/Repair/runscript.txt runscript_SequenceDetail.txt',shell=True)
	os.chdir(startingDir)
