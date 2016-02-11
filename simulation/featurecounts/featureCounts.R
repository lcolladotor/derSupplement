library('Rsubread')
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

## Identify BAM files
files <- dir('/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/hisat', pattern = paste0('[0-9]R', opt$replicate, '\\.bam$'), full.names = TRUE)
names(files) <- gsub('\\.bam', '', dir('/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/hisat', pattern = paste0('[0-9]R', opt$replicate, '\\.bam$')))


message(paste(Sys.time(), 'running featureCounts'))
featCounts <- featureCounts(files = files, annot.ext = '/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/gtf/chr17.gtf', isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

## Fix sample names
colnames(featCounts$counts) <- gsub('.bam', '', gsub('X.dcs01.ajaffe.Brain.derRuns.derSupplement.simulation.hisat.', '', colnames(featCounts$counts)))

colnames(featCounts$stat) <- gsub('.bam', '', gsub('X.dcs01.ajaffe.Brain.derRuns.derSupplement.simulation.hisat.', '', colnames(featCounts$stat)))

## Print stats
stat <- t(featCounts$stat[, -1])
colnames(stat) <- featCounts$stat[, 1]
rownames(stat) <- colnames(featCounts$stat)[-1]

## Summary first, then all the info
summary(stat)
options(width = 120)
stat

message(paste(Sys.time(), 'saving featureCounts output'))
save(featCounts, stat, file = paste0('featureCounts-R', opt$replicate, '.Rdata'))

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
