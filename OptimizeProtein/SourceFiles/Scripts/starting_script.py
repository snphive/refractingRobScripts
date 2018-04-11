import os
from OptProt import OptProt

option_file = open('./OptProt_Options.txt', 'r').readlines()
start_path = os.getcwd()
scripts_path = start_path + '/SourceFiles/Scripts'
OptProt = OptProt(start_path, scripts_path)
OptProt.parse_option_file(option_file)

if not os.path.exists('Results'):
    os.makedirs('Results')

OptProt.run_yasara_agadir_repair()
OptProt.wait_on_repair()
OptProt.perform_selected_computations()