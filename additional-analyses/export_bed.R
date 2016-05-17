## Export ER-level DERs for BrainSpan to a BED file
# qrsh -l mem_free=130G,h_vmem=150G
# module load R/3.3
# mkdir -p logs
# Rscript export_bed.R > logs/export_bed_log.txt 2>&1

library('GenomicRanges')
library('rtracklayer')
library('limma')

## Load pheno data
load('/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda')
bad_samples <- which(rownames(pdSpan) %in% c('216', '218', '219'))
if(nrow(pdSpan) == 487) pdSpan <- pdSpan[-bad_samples, ]
stopifnot(nrow(pdSpan) == 484)

pdSpan$bigwig <- paste0('http://download.alleninstitute.org/brainspan/MRF_BigWig_Gencode_v10/bigwig/', gsub('/nexsan2/disk3/ajaffe/BrainSpan/RNAseq/bigwig/', '', pdSpan$wig))
pd_sub <- pdSpan[, c('gender', 'lab', 'Age', 'structure_acronym', 'structure_name', 'bigwig')]

## Save pheno data subset
write.table(pd_sub, file = 'BrainSpan_bigwig_files.tsv', row.names = FALSE, sep = '\t', quote = FALSE)

## Load models
load('/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/derAnalysis/run5-v1.5.30/models.Rdata')
stopifnot(nrow(models$mod) == nrow(models$mod0))
if(unique(sapply(models, nrow)) == 487) {
    models$mod <- models$mod[-bad_samples, ]
    models$mod0 <- matrix(models$mod0[-bad_samples, ], ncol = 1)
}
stopifnot(nrow(models$mod) == 484)

## Load regions
load('/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.25.Rdata')
regList <- lapply(regionMat, function(x) x$regions)
fullRegionGR <- unlist(GRangesList(regList))

# coverage matrix
fullRegionMat = do.call('rbind',
	lapply(regionMat, function(x) x$coverageMatrix))
## Drop bad samples
if(ncol(fullRegionMat) == 487) fullRegionMat <- fullRegionMat[, -bad_samples]

## drop regions shorter than 6 bp
keepIndex <- which(width(fullRegionGR) >= 6)
fullRegionGR <- fullRegionGR[keepIndex]
fullRegionMat <- fullRegionMat[keepIndex,]

## log transform
y <- log2(fullRegionMat + 1)
rownames(y) <- NULL

## DE analysis
fit <- lmFit(y, models$mod)
fit0 <- lmFit(y, models$mod0)
ff <- getF(fit,fit0, y)

print('Number and percent of ER-level DERs that are significant')
sum(p.adjust(ff$f_pval, 'bonf') < 0.05)
round(mean(p.adjust(ff$f_pval, 'bonf') < 0.05) * 100, 2)
sigIndex <- which(p.adjust(ff$f_pval, 'bonf') < 0.05)

## DERs from ER-level
ders <- fullRegionGR[sigIndex]

## Save DERs
export(ders, 'BrainSpan_ER_DERs_cutoff_0.25.bed', format = 'BED')


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
