import sys
from Agadir import Agadir

agadir_results_path = sys.argv[1]
agadir = Agadir(agadir_results_path)
agadir.run_agadir_with_fasta_files('./')

