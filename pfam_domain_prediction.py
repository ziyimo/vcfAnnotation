import datetime
import argparse
import synonymous_SNP_degeneracy as sSd

threeToOne={'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Asx':'B','Cys':'C','Glu':'E','Gln':'Q','Glx':'Z','Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V'}
message=''

#Implement command line parser
parser=argparse.ArgumentParser(description='Placeholder for program description') #write a short description!!
parser.add_argument('fastaDir', type=str, help='Specify the proteome fasta file directory')
parser.add_argument('pfamDir',type=str,help='Specify the pfam tsv file directory')
parser.add_argument('vcfDir',type=str,help='Specify the vcf file directory')
parser.add_argument('-t','--targetDir',type=str,default='',help='Specify the target directory for the updated vcf file and the log file')
parser.add_argument('-a','--annformat',action='store_true',help='run the program on ann formatted SNPs')
parser.add_argument('-s','--synSNPdegeneracy',action='store_true',help='annotate the synonymous SNPs with degeneracy')
args=parser.parse_args()
if args.targetDir != '':
	args.targetDir+='/'

class VCF:
	# parsing the vcf file when instantiating a VCF object
	def __init__(self,vcfFile):
		self.metaInfo=[]
		self.headerLine=[]
		self.dataFields=[]
		for line in vcfFile:
			if line[:2]=='##':
				self.metaInfo.append(line.strip())
			elif line[:1]=='#':
				self.headerLine=line.strip().split('\t')
			else:
				self.dataFields.append(line.strip().split('\t'))

	# write the information in a VCF object in a text file
	def writeVCFFile(self,newFile):
		newFile.write('\n'.join(self.metaInfo)+'\n')
		newFile.write('\t'.join(self.headerLine)+'\n')
		for SNPinfo in self.dataFields:
			newFile.write('\t'.join(SNPinfo)+'\n')


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
				if 'SNPEFF_PFAMID' not in INFODict:
					INFODict['SNPEFF_PFAMID']=entry[4]
					INFODict['SNPEFF_DOMAINCOOR']=entry[6]+'-'+entry[7]
					logStats['matchedNonsynSNPrcd']+=1
				else:
					INFODict['SNPEFF_PFAMID']= INFODict['SNPEFF_PFAMID']+','+entry[4]
					INFODict['SNPEFF_DOMAINCOOR']= INFODict['SNPEFF_DOMAINCOOR']+','+entry[6]+'-'+entry[7]
	elif args.synSNPdegeneracy and INFODict['SNPEFF_EFFECT']=='SYNONYMOUS_CODING':
		logStats['synSNPrcd']+=1
		INFODict['SNPEFF_DEGENERACY']=str(sSd.main(INFODict['SNPEFF_CODON_CHANGE']))

	if args.synSNPdegeneracy and 'SNPEFF_DEGENERACY' not in INFODict:
		INFODict['SNPEFF_DEGENERACY']='.'
	if 'SNPEFF_PFAMID' not in INFODict:
		INFODict['SNPEFF_PFAMID']='.'
		INFODict['SNPEFF_DOMAINCOOR']='.'
	
	newINFOList=[]
	for item in list(INFODict.items()):
		newINFOList.append(item[0]+'='+item[1])
	newINFOList.sort()

	return ';'.join(newINFOList)

def mainAlt(ANNField,SNPentry):

	''' The alternative version of the main function for ANN tagged vcf files. 
	'''

	pfamIDTag='.'
	coorTag='.'
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
				if pfamIDTag == '.':
					pfamIDTag=entry[4]
					coorTag=entry[6]+'-'+entry[7]
					logStats['matchedNonsynSNPrcd']+=1
				else:
					pfamIDTag= pfamIDTag+','+entry[4]
					coorTag= coorTag+','+entry[6]+'-'+entry[7]

	ANNList.insert(14,pfamIDTag)
	ANNList.insert(15,coorTag)

	return '|'.join(ANNList)

logStats={'NonsynSNPrcd':0,'matchedNonsynSNPrcd':0}
if args.synSNPdegeneracy:
	logStats['synSNPrcd']=0

#Readind records from proteome fasta file
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
vcfData=VCF(vcfFile)
vcfFile.close()

#Performing main task on every entry
vcf_tscptID=[]
nonsynSNP_tscptID=[]
if args.annformat:
	for entry in vcfData.dataFields:
		infoList=entry[7].split(';')
		for item in infoList:
			if item[:3]=='ANN':
				ANNFields=item[4:].split(',')
				for ANNField in ANNFields:
					updatedANNField=mainAlt(ANNField,entry)
					ANNFields[ANNFields.index(ANNField)]=updatedANNField
				infoList[infoList.index(item)] = 'ANN='+','.join(ANNFields)
		entry[7]=';'.join(infoList)

else:
	for entry in vcfData.dataFields:
		updatedINFO=main(entry)
		entry[7]=updatedINFO

#Adding meta information lines

if args.annformat:
	for line in vcfData.metaInfo:
		if line[:14]=='##INFO=<ID=ANN':
			vcfData.metaInfo[vcfData.metaInfo.index(line)]='##INFO=<ID=ANN,Number=.,Type=String,Description=\"Predicted effects for this variant.Format: \'Allele|Annotation|Putative_impact|Gene_Name|Gene_ID|Feature_type|Feature_ID|Transcript_biotype|Rank/total|HGVS.c|HGVS.p|cDNA_position/cDNA_len|CDS_position/CDS_len|Protein_position/Protein_len|Pfam_ID|Protein_domain_coordinate[|Distance_to_feature|Errors_Warnings_InfoMessages]\'\">'
			break
else:
	for line in vcfData.metaInfo:
		if line[:6]!='##INFO' and vcfData.metaInfo[vcfData.metaInfo.index(line)-1][:6]=='##INFO':
			vcfData.metaInfo.insert(vcfData.metaInfo.index(line),'##INFO=<ID=SNPEFF_PFAMID,Number=.,Type=String,Description="Pfam ID of protein domain affected">')
			vcfData.metaInfo.insert(vcfData.metaInfo.index(line),'##INFO=<ID=SNPEFF_DOMAINCOOR,Number=.,Type=String,Description="Coordinates of protein domain affected">')
			if args.synSNPdegeneracy:
				vcfData.metaInfo.insert(vcfData.metaInfo.index(line),'##INFO=<ID=SNPEFF_DEGENERACY,Number=1,Type=Integer,Description="Degeneracy of synonymous SNP">')
			break

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
