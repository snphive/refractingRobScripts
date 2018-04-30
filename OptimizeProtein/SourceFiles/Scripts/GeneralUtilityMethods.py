import subprocess
import time


class GUM(object):

    @staticmethod
    def wait_for_grid_engine_job_to_complete(prefix, message_to_print):
        check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
        output_qstat = check_qstat.stdout.read()
        while prefix in output_qstat:
            print 'Waiting for ' + message_to_print
            time.sleep(10)
            check_qstat = subprocess.Popen('qstat', stdout=subprocess.PIPE)
            output_qstat = check_qstat.stdout.read()
