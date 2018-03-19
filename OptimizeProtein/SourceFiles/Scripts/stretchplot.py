# For Whatif
import os,sys,glob,subprocess
startingDir = os.getcwd()

Mols = sys.argv[1].split('_')
R_Path = sys.argv[2]
FoldX_Path = sys.argv[3]
Agadir_Path = sys.argv[4]

# amino acid dictionaries
aa_dict_3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
aa_dict_1to3 = {'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR'}

agadfolders = []

agadirfolders_temp = sorted(glob.glob('./Agadir/*'))
for agad in agadirfolders_temp:
	temp = agad.split('/')[-1]
	if temp == 'Options.txt':
		continue
	agadfolders.append(agad)

name = agadfolders[0].split('/')[-1].split('_')[0]

Seqdet = open('./Repair/SequenceDetail_RepairPDB_'+name+'.fxout','r').readlines()
headers = Seqdet[8].split('\t')
energy_header = 'total energy'
mol_header = 'chain'
resnum_header = 'number'
AA_header = 'amino acid'

for x, head in enumerate(headers):
	if head == energy_header:
		index_energy = x
	if head == mol_header:
		index_mol = x
	if head == resnum_header:
		index_resnum = x
	if head == AA_header:
		index_AA = x

h = open('SummaryPerResSeqDet.txt','w')
for line in Seqdet[8:]:
	h.write(line)
h.close()

h = open('SummaryPerResAgad.txt', 'w')
h.write('Name\tMol\tResnum\tRes1\tTANGO\tAPRcount\n')
APRcount = 0
APR = False
for agad in agadfolders:
	print agad
	agad_lines = open(agad+'/PSX_globalresidue.out','r').readlines()[1:]
	for agad_res_line in agad_lines:
		pieces = agad_res_line.split('\t')
		nametemp = pieces[0]
		moltemp = nametemp.split('_')[-1]
		resnum = pieces[1]
		res1 = pieces[2]
		TANGOstring = pieces[3]
		TANGOscore = float(TANGOstring)
		if TANGOscore > 5:
			if APR == False:
				APR = True
				APRnew = True
				APRcount += 1
			else:
				APRnew = False
		if TANGOscore < 5 and APR == True:
			EndOfAPR = 1
			APR = False
		else:
			EndOfAPR = 0
		if APR == False:
			APRcountstring = '0'
		else:
			APRcountstring = str(APRcount)
		h.write(nametemp+'\t'+moltemp+'\t'+resnum+'\t'+res1+'\t'+TANGOstring+'\t'+APRcountstring+'\n')
h.close()

g = open('SummaryStretches.txt','w')
g.write('Name\tMol\tStretch\tTangoSum\tdGSum\n')
for mol in Mols:
	f = open('./Agadir/'+name+'_'+mol+'/PSX_tangowindow.out','r').readlines()
	allstretchlines = []
	for tangoline in f[1:]:
		allstretchlines.append(tangoline)
	for stretchline in allstretchlines:
		stretchAA = stretchline.split()[3]
		stretchlen = stretchline.split()[-1]
		stretchscore = stretchline.split()[-2]
		stretchfirst = int(stretchline.split()[1])+1
		aa_3 = []
		for aa in stretchAA:
			aa_3.append(aa_dict_1to3[aa])
		total_energy_stretch = float(0)
		for x,aa in enumerate(aa_3):
			resnum = str(stretchfirst+x)
			for seqline in Seqdet[9:]:
				seqline_pieces = seqline.split('\t')
				if seqline != "" and len(seqline_pieces)==len(headers):
					if aa == seqline_pieces[index_AA] and mol == seqline_pieces[index_mol]:
						total_energy_stretch = total_energy_stretch + float(seqline_pieces[index_energy])
		TangoSum = float(stretchscore)*float(stretchlen)
		if total_energy_stretch < -5:
			total_energy_stretch = -5
		if total_energy_stretch > 5:
			total_energy_stretch = 5
		g.write(name+'\t'+mol+'\t'+stretchAA+'\t'+str(TangoSum)+'\t'+str(total_energy_stretch)+'\n')
g.close()
