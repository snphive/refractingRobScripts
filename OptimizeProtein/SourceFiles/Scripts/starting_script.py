import os
from OptProt import OptProt

option_file = open('./OptProt_Options.txt', 'r').readlines()
start_path = os.getcwd()
scripts_path = start_path + '/SourceFiles/Scripts'
opt_prot_instance = OptProt(start_path, scripts_path)
opt_prot_instance.parse_option_file(option_file)
opt_prot_instance.run_yasara_agadir_repair()
opt_prot_instance.perform_selected_computations()
