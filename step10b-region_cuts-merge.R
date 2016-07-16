library('derfinder')
library('getopt')
library('devtools')


## Specify parameters
spec <- matrix(c(
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Load data
chrs <- paste0('chr', c(1:22, 'X', 'Y'))
region_cuts_raw <- lapply(chrs, function(chr) {
    message(paste(Sys.time(), 'loading data for', chr))
    load(paste0('region_cuts-', chr, '.Rdata'))
    res <- region_cuts
    return(res)
})

## Merge
message(paste(Sys.time(), 'merging data'))
region_cuts <- lapply(seq(0.025, 0.5, by = 0.025), function(cutoff) {
    regions <- lapply(region_cuts_raw, '[[', as.character(cutoff))
    res <- do.call(c, regions)
    return(res)
})
names(region_cuts) <- seq(0.025, 0.5, by = 0.025)

## Save
save(region_cuts, file = 'region_cuts.Rdata')

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
