## Get summary information for paper

## Load libraries

## Available from http://cran.r-project.org/web/packages/getopt/index.html
# install.packages("getopt")
library("getopt")

## Available from http://cran.at.r-project.org/web/packages/rmarkdown/index.html
# install.packages("rmarkdown")
library('rmarkdown')

# install.packages("knitr")
library('knitr')

## Available from http://www.bioconductor.org/packages/release/bioc/html/GenomicRanges.html
# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicRanges")
suppressMessages(library("GenomicRanges"))

## Available from http://www.bioconductor.org/packages/release/bioc/html/ggbio.html
# source("http://bioconductor.org/biocLite.R")
# biocLite("ggbio")
suppressMessages(library("ggbio"))

## Available from
# source("http://bioconductor.org/biocLite.R")
# biocLite("derfinder")
suppressMessages(library("derfinder"))

## Available from
# source("http://bioconductor.org/biocLite.R")
# biocLite("derfinderPlot")
suppressMessages(library("derfinderPlot"))

## Available from http://www.bioconductor.org/packages/release/bioc/html/TxDb.Hsapiens.UCSC.hg19.knownGene.html
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))

## Specify parameters
spec <- matrix(c(
	'short', 's', 1, "character", "Short name of project, for example 'Hippo'",
	'run', 'r', 1, "character", "Name of the run, for example 'run1-v0.0.42'",
	'example', 'p', 1, "character", "Ids of the example cluster plots to make, for example 'c(3, 4, 8)'",
	'verbose' , 'v', 2, "logical", "Print status updates",
	'help' , 'h', 0, "logical", "Display help"
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## Testing the script
test <- FALSE
if(test) {
	## Speficy it using an interactive R session and testing
	test <- TRUE
}

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Default value for verbose = TRUE
if (is.null(opt$verbose)) opt$verbose <- TRUE
	
## Save start time for getting the total processing time
startTime <- Sys.time()

## Paths
rootdir <- '/dcl01/lieber/ajaffe/derRuns/derSupplement'
resdir <- file.path(rootdir, opt$short, 'summaryInfo', opt$run)

## results path
dir.create(resdir, recursive=TRUE)

opts_chunk$set(dev = 'CairoPNG')
render(file.path(rootdir, 'step9-summaryInfo.Rmd'), output_file=file.path(resdir, 'summaryInfo.html'), clean = FALSE)


## Done
if(opt$verbose) {
	print(proc.time())
	print(sessionInfo(), locale=FALSE)
}
