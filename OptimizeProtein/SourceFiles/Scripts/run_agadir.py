import os
from Agadir import Agadir

cwd = os.getcwd()
print 'cwd: ' + cwd
agadir_instance = Agadir(cwd)
agadir_instance.run_agadir_with_fasta_files('./')
