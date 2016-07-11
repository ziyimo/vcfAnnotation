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

	def addMetaInfo(self,metaInfoLine):
		for i in range(1,len(self.metaInfo)):
			if self.metaInfo[i][:6]!='##INFO' and self.metaInfo[i-1][:6]=='##INFO':
				self.metaInfo.insert(i,metaInfoLine)
				break

	# write the information in a VCF object in a text file
	def writeVCFFile(self,newFile):
		newFile.write('\n'.join(self.metaInfo)+'\n')
		newFile.write('\t'.join(self.headerLine)+'\n')
		for SNPinfo in self.dataFields:
			newFile.write('\t'.join(SNPinfo)+'\n')

def grab(infoField,infoTag):
	'''takes a complete info record and the name of an INFO tag and returns the value of that tag'''
	infoList=infoField.split(';')
	for item in infoList:
		if item[:item.find('=')]==infoTag:
			return item[item.find('=')+1:]