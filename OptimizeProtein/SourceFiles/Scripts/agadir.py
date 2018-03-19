import os,glob,subprocess
fastafolder = 'Fasta'
agadir = 'Agadir'
startingDir = os.getcwd()
fastafiles = sorted(glob.glob('./'+fastafolder+'/*'))
for j in fastafiles:
	os.chdir(agadir)
	name = j.split('/')[-1].split('.')[0]
	if not os.path.exists(name):
		os.makedirs(name)
	os.chdir(name)
	print os.getcwd()
	subprocess.call('agadirwrapper ./../.'+j+' ./../Options.txt',shell=True)
	os.chdir(startingDir)

