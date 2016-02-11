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
files <- rawFiles(datadir = file.path('..', 'hisat'), samplepatt = paste0('R', opt$replicate, '.bam$'), fileterm = NULL)
names(files) <- gsub('.bam', '', names(files))

## Find number of total mapped reads
load(file.path('..', 'generateReads', 'simulation_info.Rdata'))
percentMapped <- as.numeric(gsub('% overall alignment rate', '', 
    sapply(names(files), function(x) { system(paste0("grep 'overall alignment rate' ../hisat/logs/", x, '.hisat.e*'), intern = TRUE) } )
))

## Note that reads are paired-end
totalMapped <- round(colSums(readmat[, names(files)]) * percentMapped / 100, 0)

## Load coverage info, normalize to 40 million paired-end reads libraries
## (same as 80 million single-end)
fullCov <- fullCoverage(files, chrs = 'chr17', totalMapped = totalMapped, targetSize = 40e6, mc.cores = 4L)
print(object.size(fullCov), units = 'Mb')


## Save coverage info
save(fullCov, file = paste0('fullCov-R', opt$replicate, '.Rdata'))


## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
