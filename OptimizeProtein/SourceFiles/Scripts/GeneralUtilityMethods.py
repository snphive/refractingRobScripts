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

    @staticmethod
    def build_runscript_for_pdbs(path_to_runscript, pdbs, action, should_print_networks):
        runscript_file = open(path_to_runscript + 'runscript.txt', 'w')
        runscript_file.write('<TITLE>FOLDX_runscript;\n')
        runscript_file.write('<JOBSTART>#;\n')
        runscript_file.write('<PDBS>' + pdbs + ';\n')
        runscript_file.write('<BATCH>#;\n')
        runscript_file.write('<COMMANDS>FOLDX_commandfile;\n')
        runscript_file.write(action + ';\n')
        if should_print_networks:
            runscript_file.write('<PrintNetworks>#;\n')
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

    @staticmethod
    def build_job_q_bash(grid_engine_job_name, queue, max_memory, cluster, foldx_path, python_script_with_path):
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
        g.write(foldx_path + ' -runfile runscript.txt\n')
        if python_script_with_path != '':
            g.write('python ' + python_script_with_path + '\n')
        g.close()
