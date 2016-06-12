import datetime

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



def main(INFOField):

	'''The funtion takes the INFO field of a SNP (as a string). If the SNPEFF_IMPACT value is MODIFIER, the transcript ID
	the SNP associates with will be matched with Pfam data (generated from a single tsv file). Two tags are added to the 
	INFO field, SNPEFF_PFAMID contains the pfam ID of all domains that the SNP can affect, whereas SNPEFF_DOMAINCOOR specifies
	the coordinates of the domains affected within the transcript. The function returns the updated INFO field (as a string).
	In addition, a list of all unique transcript ID records in the vcf file is compiled.
	'''

	INFOList=INFOField.split(';')
	INFODict={}
	for item in INFOList:
		a=item.split('=')
		INFODict[a[0]]=a[1]
	if INFODict['SNPEFF_TRANSCRIPT_ID'] not in vcf_tscptID:
		vcf_tscptID.append(INFODict['SNPEFF_TRANSCRIPT_ID'])
	if INFODict['SNPEFF_IMPACT']=='MODIFIER':
		logStats['ModifierSNPrcd']+=1
		if INFODict['SNPEFF_TRANSCRIPT_ID'] not in modifierSNP_tscptID:
			modifierSNP_tscptID.append(INFODict['SNPEFF_TRANSCRIPT_ID'])
		for entry in pfamData:
			if entry[0]==INFODict['SNPEFF_TRANSCRIPT_ID']:
				if 'SNPEFF_PFAMID' not in INFODict:
					INFODict['SNPEFF_PFAMID']=entry[4]
					INFODict['SNPEFF_DOMAINCOOR']=entry[6]+'-'+entry[7]
					logStats['matchedModifierSNPrcd']+=1
				else:
					INFODict['SNPEFF_PFAMID']= INFODict['SNPEFF_PFAMID']+','+entry[4]
					INFODict['SNPEFF_DOMAINCOOR']= INFODict['SNPEFF_DOMAINCOOR']+','+entry[6]+'-'+entry[7]
	if 'SNPEFF_PFAMID' not in INFODict:
		INFODict['SNPEFF_PFAMID']='.'
		INFODict['SNPEFF_DOMAINCOOR']='.'
	
	newINFOList=[]
	for item in list(INFODict.items()):
		newINFOList.append(item[0]+'='+item[1])
	newINFOList.sort()

	return ';'.join(newINFOList)

logStats={'ModifierSNPrcd':0,'matchedModifierSNPrcd':0}

#Readind transcriptID records from proteome fasta file
fasta_pepDirectory=input('Specify the proteome fasta file directory: ')
fasta_pep=open(fasta_pepDirectory,'r')
fasta_pep_tscptID=[]
for line in fasta_pep:
	if line[:1]=='>':
		fasta_pep_tscptID.append(line[1:].strip())

#Putting the pfam data in a 2 dimensional list
pfamTSVDirectory=input('Specify the Pfam tsv file directory: ')
pfamTSV=open(pfamTSVDirectory,'r')
pfamData=[]
pfam_tscptID=[]
for line in pfamTSV:
	record=line.strip().split('\t')
	if record[0] not in pfam_tscptID:
		pfam_tscptID.append(record[0])
	pfamData.append(record)
pfamTSV.close()

#Putting the SNP data into an object -- vcfData
vcfFileDirectory=input('Specify vcf file directory: ')
vcfFile=open(vcfFileDirectory, 'r')
vcfData=VCF(vcfFile)
vcfFile.close()

#Performing main task on every entry
vcf_tscptID=[]
modifierSNP_tscptID=[]
for entry in vcfData.dataFields:
	updatedINFO=main(entry[7])
	entry[7]=updatedINFO

#Adding meta information lines
for line in vcfData.metaInfo:
	if line[:6]!='##INFO' and vcfData.metaInfo[vcfData.metaInfo.index(line)-1][:6]=='##INFO':
		vcfData.metaInfo.insert(vcfData.metaInfo.index(line),'##INFO=<ID=SNPEFF_PFAMID,Number=.,Type=String,Description="Pfam ID of protein domain affected">')
		vcfData.metaInfo.insert(vcfData.metaInfo.index(line),'##INFO=<ID=SNPEFF_DOMAINCOOR,Number=.,Type=String,Description="Coordinates of protein domain affected">')
		break

#Creating new vcf file
targetDirectory=input('Specify target directory: ')
newVCF=open(targetDirectory+'/updated.vcf','w')
vcfData.writeVCFFile(newVCF)
newVCF.close()

#Prepare log file stats
logStats['SNPrcd']=len(vcfData.dataFields)
logStats['PFAMrcd']=len(pfamData)
logStats['tscptID_vcf']=len(vcf_tscptID)
logStats['tscptID_ModifierSNP']=len(modifierSNP_tscptID)
logStats['tscptID_pfam']=len(pfam_tscptID)
logStats['%matched']=str(logStats['matchedModifierSNPrcd']/logStats['ModifierSNPrcd']*100)+'%'

tscptID_vcf_fasta=0
for tscptID in vcf_tscptID:
	if tscptID in fasta_pep_tscptID:
		tscptID_vcf_fasta+=1
tscptID_pfam_fasta=0
for tscptID in pfam_tscptID:
	if tscptID in fasta_pep_tscptID:
		tscptID_pfam_fasta+=1
logStats['vcf/fasta_pep']=str(tscptID_vcf_fasta/len(vcf_tscptID)*100)+'%'
logStats['pfam/fasta_pep']=str(tscptID_pfam_fasta/len(pfam_tscptID)*100)+'%'

logtxt=str(datetime.datetime.today())+'\n'
for stat in list(logStats.items()):
	logtxt=logtxt+stat[0]+'='+str(stat[1])+'\n'
print (logtxt)

logfile=open(targetDirectory+'/log.txt','w')
logfile.write(logtxt)
logfile.close()
