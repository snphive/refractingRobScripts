import os
import sys
import glob
import subprocess

startingDir = os.getcwd()
targets = sorted(glob.glob('./PDBs/*.pdb'))
fastafolder = 'Fasta'
# amino acid dictionaries
aa_dict_3to1 = {
				'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
				'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
				'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

for pdb in targets:
	name = pdb.split('/')[-1].split('.')[0]
	f = open(pdb).readlines()
	atomlines = []
	mols = []

	for line in f:
		if 'ATOM' == line[0:4]:
			mol = line[21]
			atomlines.append(line)
			if mol not in mols:
				mols.append(mol)

	for mol in mols:
		fastalist = []
		resnum = '0'

		for line in atomlines:

			if line[21] == mol and resnum != line[22:26].strip(' '):
				resnum = line[22:26].strip(' ')
				aa = line[17:20]

				if aa in aa_dict_3to1.keys():
					fastalist.append(aa_dict_3to1[aa])
		fasta = "".join(fastalist)
		print name + '_' + mol
		print fasta
		f = open(fastafolder + '/' + name + '_' + mol + '.fasta', 'w')
		f.write('>' + name + '_' + mol + '\n')
		f.write(fasta)
		f.close()
