# vcf Annotation Toolkit v1.0

## Table of Contents
1. Overview
2. vcf Parser Module
3. Protein Domain Disruption Annotator
4. Synonymous SNP Degeneracy Annotator
5. ANN Tag Codon Change Annotator
6. ANN Tag Syntax Validator

## 1\.Overview
This vcf Annotation Toolkit includes 4 command line tools that read and write information in the `INFO` tags of [vcf](https://samtools.github.io/hts-specs/VCFv4.2.pdf) files. The meta\-information lines \(a\.k\.a\. header\) of a vcf file are updated accordingly as the file is annotated\.

The toolkit requires variant annotations by either [GATK](https://software.broadinstitute.org/gatk/) or [SnpEff](http://snpeff.sourceforge.net/index.html) \(**ANN** annotation format\)\.

## 2\.vcf Parser Module
The [vcf_parser](vcf_parser.py) module aids the reading and writing of meta\-information lines as well as `INFO` tags of vcf files\. It is implemented in all 4 scripts of this toolkit\.

### Object
```python
VCF(fileObject)
```
Reads the vcf file `fileObject` points to and returns a `VCF` object containing the parsed information\.

### Attributes
```python
VCF.metaInfo
```
A *list* of all meta\-information lines

```python
VCF.headerLine
```
A *list* of fields present in the data lines

```python
VCF.dataFields
```
A 2\-dimensional *list* of all fields of all entries

### Functions
```python
VCF.addMetaInfo(metaInfoLine)
```
Adds a new meta\-information line `metaInfoLine` to `VCF.metaInfo`

```python
VCF.updateMetaInfo(infoTag,metaInfoLine)
```
Replaces the meta\-information line of `infoTag` in `VCF.metaInfo` with `metaInfoLine`

```python
VCF.writeVCFFile(fileObject)
```
Writes the information in `VCF` object to a vcf file `fileObject` points to

```python
grab(infoField,infoTag)
```
Searches `infoField` and returns the value of `infoTag`

## 3\.Protein Domain Disruption Annotator
This tool [pfam_domain_prediction](pfam_domain_prediction.py) examines any non\-synonymous SNP record. The coordinate of the amino acid change is extracted. The transcript ID the SNP associates with will be matched with Pfam data. If the amino acid change falls within the coordinate of any protein domain, its information will be registered in `PROTEIN_DOMAIN` tag\.

### Usage
```
usage: pfam_domain_prediction.py [-h] [-t TARGETDIR] [-a]
                                 fastaDir pfamDir vcfDir

positional arguments:
  fastaDir              Specify the proteome fasta file directory
  pfamDir               Specify the pfam tsv file directory
  vcfDir                Specify the vcf file directory

optional arguments:
  -h, --help            show this help message and exit
  -t TARGETDIR, --targetDir TARGETDIR
                        Specify the target directory for the updated vcf file
                        and the log file
  -a, --annformat       run the program on ann formatted SNPs
```

### Format of Additional Input Files
* Proteome Fasta
```
>Transcript_ID (Must correspond to the identifier used in vcf variant annotation)
SEQUENCESEQUENCE...
```
  Example:
```
>LOC_Os01g01010.1
MSSAAGQDNGDTAGDYIKWMCGAGGRAGGAMANLQRGVGSLVRDIGDPCLNPSPVKGSKM
LKPEKWHTCFDNDGKVIGFR...
>LOC_Os01g01010.2
MSSAAGQDNGDTAGDYIKWM...
```
* Pfam Data \([Tab Delimited](https://github.com/ebi-pf-team/interproscan/wiki/InterProScan5OutputFormats)\)
```
Transcript_ID	Sequence_MD5_digest	Length	Analysis	Domain_Identifier	Domain_Name	N_Terminal	C_Terminal	Score	Status	Date
```
  Example:
```
LOC_Os01g01030.1	2166a1eddf92e1e6f938f3fb37b16067	593	Pfam	PF07732	Multicopper oxidase	45	152	5.2E-36	T	26-04-2016
LOC_Os01g01010.2	020b91988c997c9ec66e1f078bf4af57	616	Pfam	PF00566	Rab-GTPase-TBC domain	390	536	1.6E-32	T	26-04-2016
...
```

### Annotation
A new `INFO` tag of the following format is added
```
PROTEIN_DOMAIN=Transcript_ID|Amino_Acid_Change|Domain_ID|NTerm-CTerm
```
Example:
```
PROTEIN_DOMAIN=LOC_Os01g01010.2|A6W|PF07731|2-8
```
The tag can have multiple subfields if 1\) there are multiple domains affected by one SNP, 2\) the ANN tag has multiple non-synonymous annotations on one SNP.

A new meta\-information line documenting the new tag is added to the output vcf file:
```
##INFO=<ID=PROTEIN_DOMAIN,Number=.,Type=String,Description="Predicted protein domain disruption.Format:'Transcript_ID|Amino_Acid_Change|Domain_ID|NTerm-CTerm'">
```

### Log File

The log file contains the following items:

* The date and time when the program finished running
* `SNPrcd`: The total number of SNP records in the vcf file evaluated
* `NonsynSNPrcd`: The total number of non-synonymous SNP records
* `PFAMrcd`: The total number of pfam records evaluated
* `tscptID_vcf`: The total number of unique transcript IDs contained in the vcf file
* `tscptID_ModifierSNP`: The total number of unique transcript IDs that correspond to modifier SNPs
* `tscptID_pfam`: The total number of unique transcript IDs in the pfam record
* `vcf/fasta_pep`: Percentage of transcript IDs in the vcf file that matches the protein fasta record
* `pfam/fasta_pep`: Percentage of transcript IDs in the pfam record that matches the protein fasta record
* `matchedmodifierSNPrcd`: The number of modifier SNPs that are identified by the program to have one or more predicted pfam domain
* `%matched`: The percentage of modifier SNPs that are identified by the program to have one or more predicted pfam domain
* If the coordinate and amino acid code in the vcf entry do not match the sequence in the fasta file, a message will appear in the log file reporting the entry number of such anomaly.

## 4\.Synonymous SNP Degeneracy Annotator
This tool [synonymous_SNP_degeneracy](synonymous_SNP_degeneracy.py) calculates the degeneracy of each synonymous SNP and writes the result to a new tag `DEGENERACY`\.

### Usage
```
usage: synonymous_SNP_degeneracy.py [-h] [-a ANNFORMAT] [-t TARGETDIR] vcfDir

positional arguments:
  vcfDir                Specify the vcf file directory

optional arguments:
  -h, --help            show this help message and exit
  -a ANNFORMAT, --annFormat ANNFORMAT
                        Run the program on ann tagged SNPs, directory of a
                        fasta file containing CDS sequence should be specified
                        after this flag
  -t TARGETDIR, --targetDir TARGETDIR
                        Specify the target directory for the updated vcf file
                        and the log file
```

### Format of Additional Input File
* CDS fasta \(only required when the variant annotation is in **ANN** format\)
```
>CDS_Identifier (Must correspond to the identifier used in vcf variant annotation)
SEQUENCESEQUENCE...
```
  Example:
```
>LOC_Os01g01010.1
ATGATCTGGGGAGAAAGGCGGTAG...
>LOC_Os01g01010.2
ATGCACCCTCTTCTGTGCAGGAAA...
```

### Annotation
A new `DEGENERACY` tag is added

Example:
```
DEGENERACY=2
```
A new meta\-information line documenting the new tag is added to the output vcf file:
```
##INFO=<ID=DEGENERACY,Number=1,Type=Integer,Description="Degeneracy of synonymous SNP">'
```

### Log File
The log file contains the following items:

* The date and time when the program finished running
* `SNPrcd`: The total number of SNP records in the vcf file evaluated
* `SynSNPrcd`: The total number of synonymous SNP records
* `Anomaly`: The total number of cases where multiple synonymous variants occurs at one SNP
* Message containing the entry number of any anomaly.

## 5\.ANN Tag Codon Change Annotator
This tool [codon_change_ann](codon_change_ann.py) adds codon change information extracted from CDS sequence to the `ANN` tag. Since synonymous and non-synonymous variants are not the only variants that causes codon change, instead of filtering by variant type, the program performs the codon change annotation on ANN tags as long as the "CDS_pos/CDS_length" information is present \(such as start/stop codon gains and losses\)\.

### Usage
```
usage: codon_change_ann.py [-h] [-t TARGETDIR] vcfDir CDSfastaDir


positional arguments:
  vcfDir                Specify the vcf file directory
  CDSfastaDir           Specify the directory of a fasta file containing CDS
                        sequences

optional arguments:
  -h, --help            show this help message and exit
  -t TARGETDIR, --targetDir TARGETDIR
                        Specify the target directory for the updated vcf file
                        and the log file
```

### Format of Additional Input File
* CDS fasta
```
>CDS_Identifier (Must correspond to the identifier used in vcf variant annotation)
SEQUENCESEQUENCE...
```
  Example:
```
>LOC_Os01g01010.1
ATGATCTGGGGAGAAAGGCGGTAG...
>LOC_Os01g01010.2
ATGCACCCTCTTCTGTGCAGGAAA...
```

### Annotation
The codon change information extracted from CDS sequences is inserted to the second to last subfield of the ANN tag.
Example:
```
ANN=T|synonymous_variant|LOW|id3|GENE_id3|transcript|test0001|protein_coding|1/1|c.6C>T|p.Ile2Ile|9/2689|6/24|2/7||**atC/atT**|WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS
```
The meta\-information line of `ANN` tag is updated as follows:
```
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | **Codon_Change** | ERRORS / WARNINGS / INFO' ">
```

### Log File
The log file contains the following items:

* The date and time when the program finished running
* `SNPrcd`: The total number of SNP records in the vcf file evaluated
* `Codon_Change_Annotated`: The total number of SNP records with codon change information annotated

## 6\.ANN Tag Syntax Validator
This tool [ann_field_validator](ann_field_validator.py) checks whether the format of `ANN` tags is consistent with the specifications.

### Usage
```
usage: ann_field_validator.py [-h] [-t TARGETDIR] vcfDir

positional arguments:
  vcfDir                Specify the vcf file directory

optional arguments:
  -h, --help            show this help message and exit
  -t TARGETDIR, --targetDir TARGETDIR
                        Specify the target directory for the log file
```

### Output
* When the `ANN` tags in all entries of the vcf file pass the check, the following message shows in the log file:
> All clear\! The ANN tags in this vcf file complies to the specifications
* The program stops whenever 20 syntactic anomalies are detected and gives the following message:
> Significant inconsistency detected\! Please check the vcf file thoroughly\.
* If fewer than 20 anomalies are detected by the time all entries of the file are checked, the following message is given:
> Minor inconsistency detected!

The log file reports the entry no\. of the problematic `ANN` tags followed by specific error messages \(see table below\)\.

Implementation | ANN sub\-field | Check | Error Message
-------------- | -------------- | ----- | -------------
[x] | n/a | All 16 mandatory fields are present | Incompletion
[x] | Allele | One of A,G,C,T | NonstandardUsage:Allele
[ ] | Annotation | In Sequence Ontology Terms | NonstandardUsage:Annotation
[x] | Impact | One of 'HIGH','MODERATE','LOW','MODIFIER' | NonstandardUsage:Impact
[x] | Rank/Total | Integer/Integer | NonstandardUsage:Rank/Total
[x] | HGVS.c | Complies with HGVS format of DNA substitution variant | NonstandardUsage:HGVS.c
Partial | HGVS.p | Complies with HGVS format of protein variant | NonstandardUsage:HGVS.p
[x] | cDNA_position/cDNA_len | Integer/Integer | NonstandardUsage:cDNA_position/cDNA_len
[x] | CDS_position/CDS_len | Integer/Integer | NonstandardUsage:CDS_position/CDS_len
[x] | Protein_position/Protein_len | Integer/Integer | NonstandardUsage:Protein_position/Protein_len

The contents of other customizable sub\-fields are not checked.

### Log File
The log file contains the following items:

* The date and time when the program finished running
* Total number of vcf file entries
* Number of entries parsed
* Message and Error Messages