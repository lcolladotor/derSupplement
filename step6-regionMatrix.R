library('derfinder')
library('getopt')
library('devtools')


## Specify parameters
spec <- matrix(c(
    'maindir', 'm', 1, 'character', 'Main directory',
	'cutoff', 't', 1, 'numeric', 'Cutoff to use',
	'chr', 'c', 1, 'character', 'Chromosome under analysis',
	'readLen', 'r', 1, 'integer', 'Read length',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Input options used
maindir <- opt$maindir
cutoff <- opt$cutoff
readLen <- opt$readLen
chr <- opt$chr

message(Sys.time())
timeinfo <- NULL
timeinfo <- c(timeinfo, list(Sys.time()));

## Load data
load(file.path(maindir, 'CoverageInfo', paste0(chr, 'CovInfo.Rdata')))
fullCov <- list(get(paste0(chr, 'CovInfo')))
names(fullCov) <- chr
timeinfo <- c(timeinfo, list(Sys.time()))
proc.time()
message(Sys.time())

## run regionMatrix
regionMat <- regionMatrix(fullCov, maxClusterGap = 3000L, L = readLen, cutoff = cutoff, returnBP = FALSE)
timeinfo <- c(timeinfo, list(Sys.time()))

## Save results
save(regionMat, file=paste0('regionMat-cut', cutoff, '-', chr, '.Rdata'))
timeinfo <- c(timeinfo, list(Sys.time()))

## Save time information
save(timeinfo, file=paste0('timeinfo-', chr, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
