library('ballgown')
library('devtools')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'replicate', 'r', 1, 'integer', 'Replicate number. Either 1, 2 or 3.',
    'complete', 'c', 1, 'character', "'yes' or 'no' for whether the GTF was complete or not"
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
if(FALSE) opt <- list(replicate = 1, complete = 'yes')

## Locate files
sampleDirs <- dir(file.path('..', ifelse(opt$complete == 'yes', 'stringtie', 'stringtie-inc')), pattern = paste0('sample.*R', opt$replicate))
samplePaths <- dir(file.path('..', ifelse(opt$complete == 'yes', 'stringtie', 'stringtie-inc')), pattern = paste0('sample.*R', opt$replicate), full.names = TRUE)[grepl('assembly', sampleDirs)]
sampleDirs <- sampleDirs[grepl('assembly', sampleDirs)]
pData <- data.frame(sampleDir = sampleDirs, sample = as.integer(gsub('sample|G.*', '', sampleDirs)), group = as.integer(gsub('sample.*G|R.*', '', sampleDirs)), replicate = opt$replicate)

## Load data
bg <- ballgown(samples = samplePaths, pData = pData)
print(object.size(bg), units = 'Mb')

## Save ballgown object
save(bg, file = paste0('bg-R', opt$replicate, '-', ifelse(opt$complete == 'yes', 'comp', 'inc'), '.Rdata'))

## Perform statistical tests
stat_results <- stattest(gown = bg, feature = 'trans', getFC = TRUE, meas = 'FPKM', covariate = 'group')


## Save results
save(stat_results, file = paste0('stat-R', opt$replicate, '-', ifelse(opt$complete == 'yes', 'comp', 'inc'), '.Rdata'))


## Reproducibility info
Sys.time()
proc.time()
session_info()
