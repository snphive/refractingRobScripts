import subprocess
import time


class GUM(object):

    @staticmethod
    def wait_for_grid_engine_job_to_complete(grid_engine_job_prefix, message_to_print):
        check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
        while grid_engine_job_prefix in output_qstat:
            print 'Waiting for ' + message_to_print
            time.sleep(10)
            check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
            output_qstat = check_qstat.stdout.read()

    # The runscript.txt is an input file for FoldX that informs the which pdbs to analyse and which programs to run
    @staticmethod
    def build_runscript_for_pdbs(path_to_runscript, pdbs, show_sequence_detail, action, print_networks,
                                 calculate_stability):
        runscript_file = open(path_to_runscript + 'runscript.txt', 'w')
        runscript_file.write('<TITLE>FOLDX_runscript;\n')
        runscript_file.write('<JOBSTART>#;\n')
        runscript_file.write('<PDBS>' + pdbs + ';\n')
        runscript_file.write('<BATCH>#;\n')
        runscript_file.write('<COMMANDS>FOLDX_commandfile;\n')
        if show_sequence_detail:
            runscript_file.write('<SequenceDetail>#;\n')
        runscript_file.write(action + ';\n')
        if print_networks:
            runscript_file.write('<PrintNetworks>#;\n')
        if calculate_stability:
            runscript_file.write('<Stability>#;\n')
        runscript_file.write('<END>#;\n')
        runscript_file.write('<OPTIONS>FOLDX_optionfile;\n')
        runscript_file.write('<Temperature>298;\n')
        runscript_file.write('<IonStrength>0.05;\n')
        runscript_file.write('<ph>7;\n')
        runscript_file.write('<moveNeighbours>true;\n')
        runscript_file.write('<VdWDesign>2;\n')
        runscript_file.write('<numberOfRuns>3;\n')
        runscript_file.write('<OutPDB>#;\n')
        runscript_file.write('<END>#;\n')
        runscript_file.write('<JOBEND>#;\n')
        runscript_file.write('<ENDFILE>#;\n')
        runscript_file.close()

    # The job.q is a script that includes all necessary information for the grid engine, in terms of which computations
    # to perform as well as how to run these computations
    @staticmethod
    def build_job_q_bash(grid_engine_job_name, queue, max_memory, cluster, using_runscript, foldx_path,
                         python_script_with_path):
        g = open('./job.q', 'w')
        g.write('#!/bin/bash\n')
        g.write('#$ -N ' + grid_engine_job_name + '\n')
        g.write('#$ -V\n')
        if queue != '':
            g.write('#$ -q ' + queue + '\n')
        if max_memory != '':
            g.write('#$ -l ' + max_memory + '\n')
        if cluster != '':
            g.write('#$ -l ' + cluster + '\n')
        g.write('#$ -cwd\n')
        g.write('source ~/.bash_profile\n')
        if using_runscript:
            g.write(foldx_path + ' -runfile runscript.txt\n')
        if python_script_with_path != '':
            g.write('python ' + python_script_with_path + '\n')
        g.close()
