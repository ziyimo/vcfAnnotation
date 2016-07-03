gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',\
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',\
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',\
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',\
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',\
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',\
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',\
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

def main(codonChange):

	''' The format of input to the function should resemble this example: 'tcA/tcG'
	'''

	for base in codonChange[:3]:
		if base.isupper():
			position=codonChange.index(base)
			break

	refSeq=gencode[codonChange[:3].upper()]

	checklist=['','','','']
	bases='AGCT'
	for i in range(3):
		for j in range(4):
			if i != position:
				checklist[j]+=codonChange[i].upper()
			else:
				checklist[j]+=bases[j]

	#Developer option
	print (checklist)

	degeneracy=0
	for codon in checklist:
		if gencode[codon]==refSeq:
			degeneracy+=1

	return degeneracy

# developer tests

while True:
	a=input('Codon Change (e.g. tcA/tcG): ')
	print(main(a))