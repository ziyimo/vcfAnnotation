import argparse
import datetime
import vcf_parser as vcfp

parser=argparse.ArgumentParser(description='Placeholder for program discription')
parser.add_argument('vcfDir',type=str,help='Specify the vcf file directory')
parser.add_argument('CDSfastaDir',type=str,help='Specify the directory of a fasta file containing CDS sequences')
parser.add_argument('-t','--targetDir',type=str,default='',help='Specify the target directory for the updated vcf file and the log file')
args=parser.parse_args()
if args.targetDir != '':
	args.targetDir+='/'

CDSfasta=open(args.CDSfastaDir,'r')
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

logStats={'Codon_Change_Annotated':0}



for entry in vcfData.dataFields:
	infoList=entry[7].split(';')
	for item in infoList:
		if item[:3]=='ANN':
			ANNFields=item[4:].split(',')
			for ANNField in ANNFields:
				annList=annField.split('|')	
				if annList[12]!='':
					CDS=annList[12].split('/')
					CDSpos=int(CDS[0])
					refCodon=CDSseqDict[annList[6]][CDSpos-(CDSpos-1)%3-1:CDSpos-(CDSpos-1)%3+2] #Assuming CDS identifier is the same as Transcript ID
					codonChange=''
					for i in range(3):
						if i == (CDSpos-1)%3:
							codonChange+=refCodon[i].upper()
						else:
							codonChange+=refCodon[i].lower()
					codonChange+='/'
					for i in range(3):
						if i == (CDSpos-1)%3:
							codonChange+=annList[0]
						else:
							codonChange+=refCodon[i].lower()
					annList.insert(15,codonChange)
					logStats['Codon_Change_Annotated']+=1
				else:
					annList.insert(15,'')

				updatedANNField='|'.join(annList)
				ANNFields[ANNFields.index(ANNField)]=updatedANNField
			infoList[infoList.index(item)] = 'ANN='+','.join(ANNFields)
	entry[7]=';'.join(infoList)

vcfData.updateMetaInfo('ANN','##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: \'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | Codon_Change | ERRORS / WARNINGS / INFO\' ">')

newVCF=open(args.targetDir+'updated.vcf','w')
vcfData.writeVCFFile(newVCF)
newVCF.close()

logStats['SNPrcd']=len(vcfData.dataFields)

logtxt=str(datetime.datetime.today())+'\n'
for stat in list(logStats.items()):
	logtxt=logtxt+stat[0]+'='+str(stat[1])+'\n'

print (logtxt)

logfile=open(args.targetDir+'log.txt','w')
logfile.write(logtxt)
logfile.close()
