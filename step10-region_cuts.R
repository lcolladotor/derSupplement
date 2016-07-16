library('derfinder')
library('getopt')
library('devtools')
library('GenomicRanges')


## Specify parameters
spec <- matrix(c(
    'maindir', 'm', 1, 'character', 'Main directory',
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
readLen <- opt$readLen
chr <- opt$chr


## Load data
message(paste(Sys.time(), 'loading', paste0(chr, 'covInfo.Rdata')))
load(file.path(maindir, 'CoverageInfo', paste0(chr, 'CovInfo.Rdata')))
chrCov <- get(paste0(chr, 'CovInfo'))

## Calculate mean
message(paste(Sys.time(), 'calculating mean'))
meanCov <- Reduce('+', chrCov$coverage) / length(chrCov$coverage)

## Loop over several cutoffs
region_cuts <- lapply(seq(0.025, 0.5, by = 0.025), function(cutoff) {
    message(paste(Sys.time(), 'finding regions with cutoff', cutoff))
    
    ## Find regions
    regs <- findRegions(position = meanCov > cutoff,
        fstats = meanCov, chr = chr, cutoff = cutoff,
        maxClusterGap = 3000L, L = readLen)    
    return(regs)
})
names(region_cuts) <- seq(0.025, 0.5, by = 0.025)

## Save results
save(region_cuts, file = paste0('region_cuts-', chr, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
