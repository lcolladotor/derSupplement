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

## Load coverage info
system.time(load(paste0('fullCov-R', opt$replicate, '.Rdata')))

## Create region matrix
regionMat <- regionMatrix(fullCov, L = 100L, maxClusterGap = 3000L, returnBP = FALSE, cutoff = 2.5)
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
