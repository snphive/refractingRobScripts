import os
import glob
import subprocess


class Agadir(object):

    def __init__(self, agadir_results_path):
        self.path_Results_pdbname = agadir_results_path
        self.fasta_folder_name = 'Fasta'
        self.agadir_folder_name = 'Agadir'

    def run_agadir_with_fasta_files(self, path_to_fasta_files):
        fasta_files = sorted(glob.glob(path_to_fasta_files + self.fasta_folder_name + '/*'))

        for fasta_file in fasta_files:
            os.chdir(self.agadir_folder_name)
            pdb_name_chain = fasta_file.split('/')[-1].split('.')[0]
            if not os.path.exists(pdb_name_chain):
                os.makedirs(pdb_name_chain)
            os.chdir(pdb_name_chain)
            subprocess.call('agadirwrapper ./../.' + fasta_file + ' ./../Options.txt', shell=True)
            os.chdir(self.path_Results_pdbname)

    def get_APRs_from_tangowindow_file(self, path_to_tangowindow_file):
        tangowindow_file = open(path_to_tangowindow_file + '/PSX_tangowindow.out', 'r').readlines()
        apr_sequences = []
        for line in tangowindow_file[1:]:
            apr_sequences.append(line.split()[-4])
        return apr_sequences

    @staticmethod
    def get_APR_sequence_and_index_from_line_in_tangowindow_file(line_in_tangowindow_file):
        return [line_in_tangowindow_file.split()[-4], int(float(line_in_tangowindow_file.split()[-6]))]
