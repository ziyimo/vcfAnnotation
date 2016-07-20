import datetime
import argparse
import vcf_parser as vcfp

threeToOne={'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Asx':'B','Cys':'C','Glu':'E','Gln':'Q','Glx':'Z','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
message=''

#Implement command line parser
parser=argparse.ArgumentParser(description='Placeholder for program description') #write a short description!!
parser.add_argument('fastaDir', type=str, help='Specify the proteome fasta file directory')
parser.add_argument('pfamDir',type=str,help='Specify the pfam tsv file directory')
parser.add_argument('vcfDir',type=str,help='Specify the vcf file directory')
parser.add_argument('-t','--targetDir',type=str,default='',help='Specify the target directory for the updated vcf file and the log file')
parser.add_argument('-a','--annformat',action='store_true',help='run the program on ann formatted SNPs')
args=parser.parse_args()
if args.targetDir != '':
	args.targetDir+='/'

def main(SNPentry):

	'''The function takes the INFO field of a SNP (as a string). If the SNP is non-synonymous, the coordinate of 
	the amino acid change is extracted. The transcript ID the SNP associates with will be matched with Pfam data 
	(generated from a single tsv file). If the amino acid change falls within the coordinate of a pfam domain, 
	two tags are added to the INFO field, SNPEFF_PFAMID contains the pfam ID of all domains that the SNP can affect, 
	whereas SNPEFF_DOMAINCOOR specifies	the coordinates of the domains affected within the transcript. 
	The function returns the updated INFO field (as a string).
	In addition, a list of all unique transcript ID records in the vcf file is compiled.
	'''

	INFOField=SNPentry[7]
	INFOList=INFOField.split(';')
	INFODict={}
	PROTEIN_DOMAIN_ls=[]
	for item in INFOList:
		a=item.split('=')
		INFODict[a[0]]=a[1]
	if INFODict['SNPEFF_TRANSCRIPT_ID'] not in vcf_tscptID:
		vcf_tscptID.append(INFODict['SNPEFF_TRANSCRIPT_ID'])
	if INFODict['SNPEFF_EFFECT']=='NON_SYNONYMOUS_CODING':
		logStats['NonsynSNPrcd']+=1
		if INFODict['SNPEFF_TRANSCRIPT_ID'] not in nonsynSNP_tscptID:
			nonsynSNP_tscptID.append(INFODict['SNPEFF_TRANSCRIPT_ID'])
		pepCoor=int(INFODict['SNPEFF_AMINO_ACID_CHANGE'][1:len(INFODict['SNPEFF_AMINO_ACID_CHANGE'])-1])
		refPep=INFODict['SNPEFF_AMINO_ACID_CHANGE'][0]
		if fasta_pepDict[INFODict['SNPEFF_TRANSCRIPT_ID']][pepCoor-1]!=refPep:
			global message
			message=message+'Anomaly detected at entry #'+str(vcfData.dataFields.index(SNPentry)+1)+' of the vcf file!\n'
		for entry in pfamData:
			if entry[0]==INFODict['SNPEFF_TRANSCRIPT_ID'] and int(entry[6])<=pepCoor<=int(entry[7]):
				PROTEIN_DOMAIN_ls.append(INFODict['SNPEFF_TRANSCRIPT_ID']+'|'+INFODict['SNPEFF_AMINO_ACID_CHANGE']+'|'+entry[4]+'|'+entry[6]+'-'+entry[7])
	if PROTEIN_DOMAIN_ls==[]:
		newINFOField=INFOField+';'+'PROTEIN_DOMAIN=.'
	else:
		newINFOField=INFOField+';'+'PROTEIN_DOMAIN='+','.join(PROTEIN_DOMAIN_ls)
		logStats['matchedNonsynSNPrcd']+=1

	return newINFOField

def mainAlt(ANNField,SNPentry):

	''' The alternative version of the main function for ANN tagged vcf files. 
	'''
	DD_ls=[]
	ANNList=ANNField.split('|')
	if ANNList[6] not in vcf_tscptID:
		vcf_tscptID.append(ANNList[6])
	if 'missense_variant' in ANNList[1]:
		logStats['NonsynSNPrcd']+=1
		if ANNList[6] not in nonsynSNP_tscptID:
			nonsynSNP_tscptID.append(ANNList[6])
		pepCoor=int(ANNList[10][5:len(ANNList[10])-3])
		refPep=threeToOne[ANNList[10][2:5]]
		if ANNList[6] in fasta_pepDict and fasta_pepDict[ANNList[6]][pepCoor-1]!=refPep:
			global message
			message=message+'Anomaly detected at entry #'+str(vcfData.dataFields.index(SNPentry)+1)+' of the vcf file!\n'
		for entry in pfamData:
			if entry[0]==ANNList[6] and int(entry[6])<=pepCoor<=int(entry[7]):
				AA_Change=refPep+str(pepCoor)+threeToOne[ANNList[10][len(ANNList[10])-3:]]
				DD_ls.append(ANNList[6]+'|'+AA_Change+'|'+entry[4]+'|'+entry[6]+'-'+entry[7])

	return DD_ls

logStats={'NonsynSNPrcd':0,'matchedNonsynSNPrcd':0}

#Reading records from proteome fasta file
fasta_pep=open(args.fastaDir,'r')
fasta_pepDict={}
for line in fasta_pep:
	if line[:1]=='>':
		tscptID=line[1:].strip()
		fasta_pepDict[tscptID]=''
	else:
		fasta_pepDict[tscptID]+=line.strip()

#Putting the pfam data in a 2 dimensional list
pfamTSV=open(args.pfamDir,'r')
pfamData=[]
pfam_tscptID=[]
for line in pfamTSV:
	record=line.strip().split('\t')
	if record[0] not in pfam_tscptID:
		pfam_tscptID.append(record[0])
	pfamData.append(record)
pfamTSV.close()

#Putting the SNP data into an object -- vcfData
vcfFile=open(args.vcfDir, 'r')
vcfData=vcfp.VCF(vcfFile)
vcfFile.close()

#Performing main task on every entry
vcf_tscptID=[]
nonsynSNP_tscptID=[]
if args.annformat:
	for entry in vcfData.dataFields:
		infoList=entry[7].split(';')
		PROTEIN_DOMAIN_ls=[]
		for item in infoList:
			if item[:3]=='ANN':
				ANNFields=item[4:].split(',')
				for ANNField in ANNFields:
					PROTEIN_DOMAIN_ls+=mainAlt(ANNField,entry)
		if PROTEIN_DOMAIN_ls==[]:
			entry[7]=entry[7]+';'+'PROTEIN_DOMAIN=.'
		else:
			entry[7]=entry[7]+';'+'PROTEIN_DOMAIN='+','.join(PROTEIN_DOMAIN_ls)
			logStats['matchedNonsynSNPrcd']+=1

else:
	for entry in vcfData.dataFields:
		updatedINFO=main(entry)
		entry[7]=updatedINFO

#Adding meta information lines
vcfData.addMetaInfo('##INFO=<ID=PROTEIN_DOMAIN,Number=.,Type=String,Description=\"Predicted protein domain disruption.Format:\'Transcript_ID|Amino_Acid_Change|Domain_ID|NTerm-CTerm\'\">')

#Creating new vcf file
newVCF=open(args.targetDir+'updated.vcf','w')
vcfData.writeVCFFile(newVCF)
newVCF.close()

#Prepare log file stats
logStats['SNPrcd']=len(vcfData.dataFields)
logStats['PFAMrcd']=len(pfamData)
logStats['tscptID_vcf']=len(vcf_tscptID)
logStats['tscptID_NonsynSNP']=len(nonsynSNP_tscptID)
logStats['tscptID_pfam']=len(pfam_tscptID)
logStats['%matched']=str(logStats['matchedNonsynSNPrcd']/logStats['NonsynSNPrcd']*100)+'%'

tscptID_vcf_fasta=0
for tscptID in vcf_tscptID:
	if tscptID in fasta_pepDict:
		tscptID_vcf_fasta+=1
tscptID_pfam_fasta=0
for tscptID in pfam_tscptID:
	if tscptID in fasta_pepDict:
		tscptID_pfam_fasta+=1
logStats['vcf/fasta_pep']=str(tscptID_vcf_fasta/len(vcf_tscptID)*100)+'%'
logStats['pfam/fasta_pep']=str(tscptID_pfam_fasta/len(pfam_tscptID)*100)+'%'

logtxt=str(datetime.datetime.today())+'\n'
for stat in list(logStats.items()):
	logtxt=logtxt+stat[0]+'='+str(stat[1])+'\n'
logtxt=logtxt+'\n\n'+message
print (logtxt)

logfile=open(args.targetDir+'log.txt','w')
logfile.write(logtxt)
logfile.close()
