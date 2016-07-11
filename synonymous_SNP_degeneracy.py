import argparse
import datetime
import vcf_parser as vcfp

threeToOne={'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Asx':'B','Cys':'C','Glu':'E','Gln':'Q','Glx':'Z','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',\
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',\
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',\
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',\
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',\
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',\
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',\
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

message=''

parser=argparse.ArgumentParser(description='The program populates the INFO field of each synonymous SNP with a tag indicating its degeneracy')
parser.add_argument('vcfDir',type=str,help='Specify the vcf file directory')
parser.add_argument('-a','--annFormat',type=str,help='Run the program on ann tagged SNPs, directory of a fasta file containing CDS sequence should be specified after this flag')
parser.add_argument('-t','--targetDir',type=str,default='',help='Specify the target directory for the updated vcf file and the log file')
args=parser.parse_args()
if args.targetDir != '':
	args.targetDir+='/'

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

	degeneracy=0
	for codon in checklist:
		if gencode[codon]==refSeq:
			degeneracy+=1

	return degeneracy

if args.annFormat != None:
	CDSfasta=open(args.annFormat,'r')
	CDSseqDict={}
	for line in CDSfasta:
		if line[:1]=='>':
			CDS_ID=line[1:].strip()
			CDSseqDict[CDS_ID]=''
		else:
			CDSseqDict[CDS_ID]+=line.strip()

vcfFile=open(args.vcfDir, 'r')
vcfData=vcfp.VCF(vcfFile)
vcfFile.close()

logStats={'SynSNPrcd':0,'Anomaly':0}

for entry in vcfData.dataFields:
	degeneracy='.'
	if args.annFormat != None:
		annFields=vcfp.grab(entry[7],'ANN')
		if annFields.count('synonymous_variant')>1:
			logStats['SynSNPrcd']+=1
			message=message+'Multiple synonymous variant detected at entry #'+str(vcfData.dataFields.index(entry)+1)+'\n'
			logStats['Anomaly']+=1
		else:
			annFieldsls=annFields.split(',')
			for annField in annFieldsls:
				annList=annField.split('|')
				if 'synonymous_variant' in annList[1]:
					logStats['SynSNPrcd']+=1
					CDS=annList[12].split('/')
					CDSpos=int(CDS[0])
					refCodon=CDSseqDict[annList[6]][CDSpos-(CDSpos-1)%3-1:CDSpos-(CDSpos-1)%3+2] #Assuming CDS identifier is the same as Transcript ID
					codonChange=''
					for i in range(3):
						if i == (CDSpos-1)%3:
							codonChange+=refCodon[i].upper()
						else:
							codonChange+=refCodon[i].lower()
					degeneracy=str(main(codonChange))

	else:
		if vcfp.grab(entry[7],'SNPEFF_EFFECT')=='SYNONYMOUS_CODING':
			logStats['SynSNPrcd']+=1
			degeneracy=str(main(vcfp.grab(entry[7],'SNPEFF_CODON_CHANGE')))
	
	entry[7]=entry[7]+';DEGENERACY='+degeneracy

vcfData.addMetaInfo('##INFO=<ID=DEGENERACY,Number=1,Type=Integer,Description="Degeneracy of synonymous SNP">')

newVCF=open(args.targetDir+'updated.vcf','w')
vcfData.writeVCFFile(newVCF)
newVCF.close()

logStats['SNPrcd']=len(vcfData.dataFields)

logtxt=str(datetime.datetime.today())+'\n'
for stat in list(logStats.items()):
	logtxt=logtxt+stat[0]+'='+str(stat[1])+'\n'
logtxt=logtxt+'\n\n'+message
print (logtxt)

logfile=open(args.targetDir+'log.txt','w')
logfile.write(logtxt)
logfile.close()
