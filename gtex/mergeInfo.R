## Merge GTEX info
library('devtools')

## Load data
### sample level data
pd <-  read.delim("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/44735/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v5.p1.c1.GRU/PhenotypeFiles/phs000424.v5.pht002743.v5.p1.c1.GTEx_Sample_Attributes.GRU.txt",
	as.is = TRUE, skip = 10)
#### subject level data
pheno <-  read.delim("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/44735/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v5.p1.c1.GRU/PhenotypeFiles/phs000424.v5.pht002742.v5.p1.c1.GTEx_Subject_Phenotypes.GRU.txt",
	as.is = TRUE, skip = 10)

#### linking sample or subject to Bioproject ID
matching <-  read.delim("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/44735/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v5.p1.c1.GRU/PhenotypeFiles/phs000424.v5.pht002741.v5.p1.GTEx_Sample.MULTI.txt",
	as.is = TRUE, skip = 10)

### linking bioproject ID to SRA accession
manifest <-  read.csv("/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/manifest_33618_06-15-2015_s.csv",
	as.is = TRUE)



## Start merging matching and pheno data
colnames(matching)[colnames(matching) %in% colnames(pheno)]
dim(matching)
dim(pheno)
merged <- merge(matching, pheno, all = TRUE)
dim(merged)

## Next merge in the sample data
dim(pd)
colnames(merged)[colnames(merged) %in% colnames(pd)]
merged <- merge(merged, pd, all = TRUE)
dim(merged)

## Next find how to merge the manifest data
colnames(merged)[colnames(merged) %in% colnames(manifest)]

table(manifest$gap_sample_id %in% merged$dbGaP_Sample_ID)
colnames(manifest)[colnames(manifest) == 'gap_sample_id'] <- 'dbGaP_Sample_ID'

table(gsub(' ', '', manifest$bio_sample_id) %in% merged$BioSample.Accession)
colnames(manifest)[colnames(manifest) == 'bio_sample_id'] <- 'BioSample.Accession'
manifest$BioSample.Accession <- gsub(' ', '', manifest$BioSample.Accession)

table(gsub(' ', '', manifest$submitted_sample_id) %in% merged$SAMPID)
colnames(manifest)[colnames(manifest) == 'submitted_sample_id'] <- 'SAMPID'
manifest$SAMPID <- gsub(' ', '', manifest$SAMPID)

## Proceed to merge the manifest data
colnames(merged)[colnames(merged) %in% colnames(manifest)]
dim(manifest)
merged <- merge(merged, manifest, all = TRUE)
dim(merged)

## Note that not all have files
table(is.na(merged$file))

## Save info
gtexPd <- merged
save(gtexPd, file = '/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/gtexPd.Rdata')

## Reproducibility info
proc.time()
Sys.time()
options(width = 120)
session_info()
