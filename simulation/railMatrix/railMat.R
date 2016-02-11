library('derfinder')
library('devtools')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'replicate', 'r', 1, 'integer', 'Replicate number. Either 1, 2 or 3.',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) opt <- list(replicate = 1)

## Locate files
summaryFile <- file.path('..', 'rail', paste0('rail-rna_out-R', opt$replicate), 'coverage_bigwigs', 'mean.chr17.bw')
sampleFiles <- rawFiles(datadir = file.path('..', 'rail', paste0('rail-rna_out-R', opt$replicate), 'coverage_bigwigs'), samplepatt = paste0('R', opt$replicate, '.bw$'), fileterm = NULL)
names(sampleFiles) <- gsub('.bw|-', '', names(sampleFiles))

## Find total number of mapped reads
counts <- read.table(file.path('..', 'rail', paste0('rail-rna_out-R', opt$replicate), 'cross_sample_results', 'counts.tsv.gz'), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
counts$sample <- gsub('-', '', counts$X)
counts$totalMapped <- as.numeric(sapply(strsplit(counts[["total.mapped.reads"]], ','), '[[', 2))
totalMapped <- counts$totalMapped[match(names(sampleFiles), counts$sample)]
names(totalMapped) <- names(sampleFiles)

## Create region matrix
regionMat <- railMatrix(chrs = 'chr17', summaryFiles = summaryFile, sampleFiles = sampleFiles, L = 100L, maxClusterGap = 3000L, returnBP = FALSE, cutoff = 2.5, totalMapped = totalMapped, targetSize = 40e6, chunksize = 15e3)
print(object.size(regionMat), units = 'Mb')

## Basic summary info
summary(width(regionMat$chr17$regions))
length(regionMat$chr17$regions)

## Save results
save(regionMat, file = paste0('regionMat-R', opt$replicate, '.Rdata'))

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
