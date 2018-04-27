from OptProt import OptProt
import yaml

option_file = open('./OptProt_Options.txt', 'r').readlines()
opt_prot_instance = ''

with open("/switchlab/group/shazib/OptimizeProteinShazibCopy/SourceFiles/Scripts/pathsAndDictionaries.yaml",
                  'r') as stream:

    try:

        paths_and_dictionaries = yaml.load(stream)
        start_path = paths_and_dictionaries['ROOT']['start_path']
        scripts_path = paths_and_dictionaries['ROOT']['scripts_path']
        r_path = paths_and_dictionaries['ROOT']['R_path']
        foldx_path = paths_and_dictionaries['ROOT']['FoldX_path']
        agadir_path = paths_and_dictionaries['ROOT']['Agadir_path']
        qsub_path = paths_and_dictionaries['ROOT']['Qsub_path']
        opt_prot_instance = OptProt(start_path, scripts_path, r_path, foldx_path, agadir_path, qsub_path)

    except yaml.YAMLError as exc:
        print(exc)

opt_prot_instance.parse_option_file(option_file)
opt_prot_instance.run_yasara_agadir_repair()
opt_prot_instance.perform_selected_computations()
