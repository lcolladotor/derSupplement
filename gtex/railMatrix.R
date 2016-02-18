library('derfinder')
library('rtracklayer')
library('BiocParallel')
library('devtools')

## Options
cutoff <- 5L
maxClusterGap <- 3000L
targetSize <- 40e6

message(Sys.time())
timeinfo <- NULL
timeinfo <- c(timeinfo, list(Sys.time()));

## Get chr length
chrs <- paste0('chr', c(1:22, 'X', 'Y'))
#chrs <- paste0('chr', c('X', 'Y')) # For testing
chrInfo <- read.table('/dcl01/leek/data/gtex_work/runs/gtex/hg38.sizes', header = FALSE, stringsAsFactors = FALSE, col.names = c('chr', 'length'))
chrlens <- chrInfo$length[chrInfo$chr %in% chrs]

## Load sample info
load('/dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/gtex_pheno_with_mapped.Rdata')

## Define summary and sample files
summaryFiles <- rep('/dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/normalizedMean.bw', length(chrs))
sampleFiles <- pd2$sampleFile
names(sampleFiles) <- pd2$sra_accession

## Read length
L <- pd2$avgLength / ifelse(pd2$LibraryLayout == 'SINGLE', 1, 2)
## What is the observed read length
table(L) # They are all L = 76

## Calculate regionMatrix
regionMat <- railMatrix(chrs, summaryFiles, sampleFiles, L = L, cutoff = cutoff, targetSize = targetSize, totalMapped = pd2$totalMapped, file.cores = 8L, chunksize = 10000, verbose.load = FALSE, chrlens = chrlens, maxClusterGap = maxClusterGap)

timeinfo <- c(timeinfo, list(Sys.time()))

## Save results
save(regionMat, file=paste0('regionMat-cut', cutoff, '.Rdata'))
timeinfo <- c(timeinfo, list(Sys.time()))

## Save time information
save(timeinfo, file=paste0('timeinfo-cut', cutoff, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
