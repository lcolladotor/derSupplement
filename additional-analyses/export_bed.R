## Export ER-level DERs for BrainSpan to a BED file
# qrsh -l mem_free=130G,h_vmem=150G
# module load R/3.3
# mkdir -p logs
# Rscript analyze_brainspan.R > logs/analyze_brainspan_log.txt 2>&1

library('GenomicRanges')
library('rtracklayer')

## Load pheno data
load('/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda')
bad_samples <- which(rownames(pdSpan) %in% c('216', '218', '219'))
if(nrow(pdSpan) == 487) pdSpan <- pdSpan[-bad_samples, ]
stopifnot(nrow(pdSpan) == 484)

pdSpan$bigwig <- paste0('http://download.alleninstitute.org/brainspan/MRF_BigWig_Gencode_v10/bigwig/', gsub('/nexsan2/disk3/ajaffe/BrainSpan/RNAseq/bigwig/', '', pdSpan$wig))
pd_sub <- pdSpan[, c('gender', 'lab', 'Age', 'structure_acronym', 'structure_name', 'bigwig')]

## Save pheno data subset
write.table(pd_sub, file = 'BrainSpan_bigwig_files.tsv', row.names = FALSE, sep = '\t', quote = FALSE)


load('/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.25.Rdata')
regList <- lapply(regionMat, function(x) x$regions)
fullRegionGR <- unlist(GRangesList(regList))
export(fullRegionGR, 'BrainSpan_ER_DERs_cutoff_0.25.bed', format = 'BED')


## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
