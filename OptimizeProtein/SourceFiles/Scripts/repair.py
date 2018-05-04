import os
import sys
import glob
import subprocess
import time
import GeneralUtilityMethods
import yaml

startingDir = os.getcwd()
path_to_runscript = startingDir + '/Repair/'
targets = sorted(glob.glob('./PDBs/*.pdb'))
foldx = 'FoldX'
results = 'Results'
__repair_job_prefix__ = 'RP_'
foldx_path = ''

with open("/switchlab/group/shazib/OptimizeProteinShazibCopy/SourceFiles/Scripts/pathsAndDictionaries.yaml",
          'r') as stream:
    try:

        paths_and_dictionaries = yaml.load(stream)
        foldx_path = paths_and_dictionaries['ROOT']['FoldX_path']
        qsub_path = paths_and_dictionaries['ROOT']['Qsub_path']

    except yaml.YAMLError as exc:
        print(exc)

Qsub_Path = sys.argv[1]
for pdb in targets:
    name = pdb.split('/')[-1].split('.')[0]
    rawpdb = pdb.split('/')[-1]
    os.chdir('Repair')
    if not os.path.exists('RepairPDB_'+rawpdb):
        # path_to_runscript = startingDir + '/Repair/'
        pdbs = rawpdb
        show_sequence_detail = False
        action = '<RepairPDB>#'
        print_networks = True
        calculate_stability = False
        GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, pdbs, show_sequence_detail, action,
                                                           print_networks, calculate_stability)
        grid_engine_job_name = __repair_job_prefix__ + name
        queue = 'all.q'
        max_memory = 'h_vmem=3G'
        cluster = 'hostname=hodor1.vib'
        no_python_script = ''
        GeneralUtilityMethods.GUM.build_job_q_bash(grid_engine_job_name, queue, max_memory, cluster, foldx_path,
                                                   no_python_script)
        subprocess.call(Qsub_Path+'qsub job.q',shell=True)

        message_to_print = 'RepairPDB to finish'
        GeneralUtilityMethods.GUM.wait_for_grid_engine_job_to_complete(__repair_job_prefix__, message_to_print)

        subprocess.call('cp '+startingDir+'/Repair/runscript.txt runscript_repair.txt',shell=True)

        pdbs = 'RepairPDB_' + rawpdb
        show_sequence_detail = True
        action = '<AnalyseComplex>#'
        print_networks = True
        calculate_stability = True
        GeneralUtilityMethods.GUM.build_runscript_for_pdbs(path_to_runscript, pdbs, show_sequence_detail, action,
                                                           print_networks, calculate_stability)

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
