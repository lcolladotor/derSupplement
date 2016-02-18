## Select samples to study
## Usage:
# mkdir -p logs
# Rscript select_samples.R > logs/select_samples_log.txt 2>&1
library('devtools')

load('/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/gtexPd.Rdata')
dim(gtexPd)
gtexPd$sra_accession <- gsub(' ', '', gtexPd$sra_accession)

## Subset to valid options: they have a SRR, are from the tissues of interest, 
## and have a RIN > 7
gtex <- gtexPd[gtexPd$SMTS %in% c('Heart', 'Liver', 'Testis') & !is.na(gtexPd$sra_accession) & gtexPd$SMRIN > 7, ]
dim(gtex)

valid_subj <- tapply(gtex$SMTSD, as.factor(gtex$dbGaP_Subject_ID), function(x) {
    all(c('Heart - Left Ventricle', 'Liver', 'Testis') %in% x)
})
valid_subj <- as.integer(names(valid_subj)[valid_subj])

## Number of valid subjects
length(valid_subj)

if(length(valid_subj) >= 12) {
    set.seed(20160218)
    selected_subj <- sample(valid_subj, 12)
} else {
    selected_subj <- valid_subj
}

#tapply(gtex$SMTS[gtex$dbGaP_Subject_ID %in% selected_subj], gtex$dbGaP_Subject_ID[gtex$dbGaP_Subject_ID %in% selected_subj], identity)

selected_samples <- sapply(selected_subj, function(subj) {
    ids <- gtex$dbGaP_Subject_ID
    tissue <- gtex$SMTS
    heart <- which(ids == subj & tissue == 'Heart')
    liver <- which(ids == subj & tissue == 'Liver')
    testis <- which(ids == subj & tissue == 'Testis')
    
    c(heart[1], liver[1], testis[1])
})

selected_srr <- gtex$sra_accession[as.vector(selected_samples)]
selected_srr

## Some checks
test <- gtex[selected_samples[, 1], c('SMTS', 'dbGaP_Subject_ID')]
stopifnot(all(test$SMTS %in% c('Heart', 'Liver', 'Testis')))
stopifnot(length(unique(test$dbGaP_Subject_ID)) == 1)
stopifnot(identical(
    apply(selected_samples, 2, function(x) { gtex$SMTS[x] }),
    matrix(rep(c('Heart', 'Liver', 'Testis'), length(selected_subj)), ncol = length(selected_subj))
))
stopifnot(table(gtex$SMTS[as.vector(selected_samples)]) - length(selected_subj) == rep(0, 3))

## Subset pheno data
pd1 <- gtex[gtex$sra_accession %in% selected_srr, ]
dim(pd1)

## Save results
save(pd1, file="gtex_pheno.rda")

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
