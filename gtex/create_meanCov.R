## Usage:
# mkdir -p logs
# Needs about 12-20G of RAM
# Rscript create_meanCov.R > logs/create_meanCov.txt 2>&1
library('derfinder')
library('rtracklayer')
library('GenomicRanges')
library('devtools')

## Options
targetSize <- 40e6

## Load GTEx pheno table
## This table is created by https://github.com/nellore/runs/blob/master/gtex/DER_analysis/pheno/format_pheno.R
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/pheno/pheno_missing_less_10.Rdata')

## Load pheno info from samples of interest
## This table is created by the select_samples.R script
load('/dcl01/lieber/ajaffe/derRuns/derSupplement/gtex/gtex_pheno.rda')

## Subset to use only samples of interest
pheno <- pheno[match(pd1$sra_accession, pheno$Run), ]


## Identify sample files
sampleFiles <- pheno$BigWigPath
names(sampleFiles) <- gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', sampleFiles)
stopifnot(length(sampleFiles) == nrow(pd1))

## Find count files
counts_files <- file.path(dir('/dcl01/leek/data/gtex', pattern = 'batch', full.names = TRUE), 'cross_sample_results', 'counts.tsv.gz')
names(counts_files) <- dir('/dcl01/leek/data/gtex', pattern = 'batch')

## Read in counts info
counts <- lapply(counts_files, function(i) {
    read.table(i, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
})
counts <- do.call(rbind, counts)
counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads, ','), '[[', 1))


## Match files to counts
map <- match(gsub('/dcl01/leek/data/gtex/batch_[0-9]*/coverage_bigwigs/|.bw', '', sampleFiles), counts$X)
counts <- counts[map, ]

## Check distribution
summary(counts$totalMapped)
summary(counts$totalMapped) / 1e6

## Create updated pheno table
pd2 <- pd1
pd2$sampleFile <- sampleFiles
pd2$totalMapped <- counts$totalMapped
pd2$SumCoverage <- pheno$SumCoverage
pd2$avgLength <- pheno$avgLength
pd2$LibraryLayout <- pheno$LibraryLayout
save(pd2, file = 'gtex_pheno_with_mapped.Rdata')

## Get chr length
chrs <- paste0('chr', c(1:22, 'X', 'Y'))
#chrs <- paste0('chr', c('X', 'Y')) # For testing
chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
chrlens <- chrInfo$length[chrInfo$chr %in% chrs]

## Load coverage
bws <- BigWigFileList(pd2$sampleFile)
names(bws) <- pd2$sra_accession
fullCov <- fullCoverage(files = bws, chrs = chrs, totalMapped = pd2$totalMapped, targetSize = targetSize, chrlens = chrlens)

## Calculate mean
meanCov <- lapply(fullCov, function(cov) {
    DataFrame('normalizedMean' = Reduce('+', cov) / length(cov))
})

## Write BigWig file with normalized mean
createBw(meanCov)

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
