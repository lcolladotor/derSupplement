## Helper script for comparing ER-level DERs vs single base-level DERs
library('getopt')
library('rmarkdown')
library('knitr')


## Specify parameters
spec <- matrix(c(
	'run', 'r', 1, "character", "Name of the run, for example 'run1-v0.0.42'",
    'wdir', 'w', 1, 'character', 'Path to working directory',
    'maindir', 'm', 1, 'character', 'Path to main directory',
    'rootdir', 'r', 1, 'character', 'Path to root directory',
    'cutoff', 'c', 1, 'numeric', 'Cutoff used in the regionMatrix analysis',
	'help' , 'h', 0, "logical", "Display help"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

analysisPath <- opt$wdir
load(paste0(opt$maindir, '/regionMatrix/regionMat-cut', opt$cutoff, '.Rdata'))
proc.time()
load(file.path(opt$maindir, 'derAnalysis', opt$run, 'fullRegions.Rdata'))
proc.time()
opts_chunk$set(dev = 'CairoPNG')
render(file.path(opt$rootdir, 'step7-regMatVsDERs.Rmd'),
    output_file= file.path(opt$wdir, 'step7-regMatVsDERs.html'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
