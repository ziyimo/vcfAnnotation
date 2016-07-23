import argparse
import datetime
import vcf_parser as vcfp

Allele=['A','G','C','T']
Annotation=['coding_sequence_variant', 'chromosome', 'coding_sequence_variant', 'inframe_insertion', 'disruptive_inframe_insertion', 'inframe_deletion', 'disruptive_inframe_deletion', 'downstream_gene_variant', 'exon_variant', 'exon_loss_variant', 'frameshift_variant', 'gene_variant', 'intergenic_region', 'conserved_intergenic_variant', 'intragenic_variant', 'intron_variant', 'conserved_intron_variant', 'miRNA', 'missense_variant', 'initiator_codon_variant', 'stop_retained_variant', 'rare_amino_acid_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'splice_region_variant', 'splice_region_variant', 'splice_region_variant', 'stop_lost', '5_prime_UTR_premature_start_codon_gain_variant', 'start_lost', 'stop_gained', 'synonymous_variant', 'start_retained', 'stop_retained_variant', 'transcript_variant', 'regulatory_region_variant', 'upstream_gene_variant', '3_prime_UTR_variant', '3_prime_UTR_truncation+exon_loss', '5_prime_UTR_variant', '5_prime_UTR_truncation+exon_loss_variant', 'sequence_feature+exon_loss_variant']
Impact=['HIGH','MODERATE','LOW','MODIFIER']
AA_3letter=['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala', 'Tyr', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Ser', 'Gly']
HGVS_prefix=['g.','m.','c.','n.']
nucleotide_numbering=['0','1','2','3','4','5','6','7','8','9','+','-','*']

all_clear=True
message=''
anomalyTally=0

parser=argparse.ArgumentParser(description='Placeholder for program discription')
parser.add_argument('vcfDir',type=str,help='Specify the vcf file directory')
parser.add_argument('-t','--targetDir',type=str,default='',help='Specify the target directory for the log file')
args=parser.parse_args()
if args.targetDir != '':
	args.targetDir+='/'

def intSlashint(string):
	ls=string.split('/')
	if len(ls)==2 and ls[0].isdigit() and ls[1].isdigit():
		return True
	else:
		return False

def HGVS_cSubst(string):
	'''e.g. "c.1406G>A"'''

	std=True
	if len(string)<6:
		std=False
	elif string[:2] not in HGVS_prefix:
		std=False
	elif not string[2:-3].isdigit():
		for char in string[2:-3]:
			if char not in nucleotide_numbering:
				std=False
				break
	elif string[-3] not in Allele or string[-2]!='>' or string[-1] not in Allele:
		std=False

	return std

def HGVS_p(string):
	'''e.g. "p.Gly469Glu"'''

	std=True
	if len(string)<7:
		std=False
	elif string[:2]!='p.':
		std=False
	#elif string[2]=='*':
	#	if not string[3:-3].isdigit() or string[-3:] not in AA_3letter:
	#		std=False
	#elif string[-1]=='*':
	#	if not string[5:-1].isdigit() or string[2:5] not in AA_3letter:
	#		std=False
	#elif not string[5:-3].isdigit():
	#	std=False
	#elif string[2:5] not in AA_3letter or string[-3:] not in AA_3letter:
	#	std=False

	return std

def validate(annTag,entryNo):
	'''The input to this function are the value (string) of the ANN tag in the INFO field of a vcf entry, and the number of the entry from which the ANN tag is extracted'''

	global all_clear
	global message
	global anomalyTally

	annFields=annTag.split(',')
	ErrorLs=[]

	for annField in annFields:
		annLs=annField.split('|')
		if len(annLs) < 16:
			ErrorLs.append('Incompletion')
			all_clear=False
		else:
			if annLs[0] not in Allele:
				ErrorLs.append('NonstandardUsage:Allele')
				all_clear=False
			#for effect in annLs[1].split('&'):
			#	if effect not in Annotation:
			#		print (effect)
			#		ErrorLs.append('NonstandardUsage:Annotation')
			#		all_clear=False
			#		break
			if annLs[2] not in Impact:
				ErrorLs.append('NonstandardUsage:Impact')
				all_clear=False
			if annLs[8]!='' and not intSlashint(annLs[8]):
				ErrorLs.append('NonstandardUsage:Rank/Total')
				all_clear=False
			if annLs[9]!='' and not HGVS_cSubst(annLs[9]):
				ErrorLs.append('NonstandardUsage:HGVS.c')
				all_clear=False
			if annLs[10]!='' and not HGVS_p(annLs[10]):
				ErrorLs.append('NonstandardUsage:HGVS.p')
				all_clear=False
			if annLs[11]!='' and not intSlashint(annLs[11]):
				ErrorLs.append('NonstandardUsage:cDNA_position/cDNA_len')
				all_clear=False
			if annLs[12]!='' and not intSlashint(annLs[12]):
				ErrorLs.append('NonstandardUsage:CDS_position/CDS_len')
				all_clear=False
			if annLs[13]!='' and not intSlashint(annLs[13]):
				ErrorLs.append('NonstandardUsage:Protein_position/Protein_len')
				all_clear=False

	if ErrorLs!=[]:
		message=message+'Entry#'+str(entryNo+1)+':'+','.join(ErrorLs)+'\n'
		anomalyTally+=1

vcfFile=open(args.vcfDir, 'r')
vcfData=vcfp.VCF(vcfFile)
vcfFile.close()

entryParsed=0
for entry in vcfData.dataFields:
	ann=vcfp.grab(entry[7],'ANN')
	validate(ann,vcfData.dataFields.index(entry))
	entryParsed+=1
	if anomalyTally==20:
		break

logtxt=str(datetime.datetime.today())+'\n'+'Total number of vcf file entries:'+str(len(vcfData.dataFields))+'\n'+'Number of entries parsed:'+str(entryParsed)+'\n'
if all_clear:
	logtxt+='All clear! The ANN tags in this vcf file complies to the specifications'
elif anomalyTally==20:
	logtxt=logtxt+'Significant inconsistency detected! Please check the vcf file thoroughly.\n'+message+'{Only showing the first 20 inconsistent entries}'
else:
	logtxt=logtxt+'Minor inconsistency detected!\n'+message

print (logtxt)

logfile=open(args.targetDir+'log.txt','w')
logfile.write(logtxt)
logfile.close()
