## Load the data without a filter, save it, then filter it for derfinder processing steps

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('BiocParallel')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'datadir', 'd', 1, 'character', 'Data directory, matched with rawFiles(datadir)',
	'pattern', 'p', 1, 'character', 'Sample pattern',
	'cutoff', 'c', 1, 'numeric', 'Filtering cutoff used',
	'mcores', 'm', 1, 'integer', 'Number of cores',
    'fileStyle', 'f', 2, 'character', 'FileStyle used for naming the chromosomes',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)

## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Default value for fileStyle
if (is.null(opt$fileStyle)) opt$fileStyle <- 'UCSC'

## Normalize filtered coverage
totalMapped <- NULL
targetSize <- 80e6


## Identify the data directories 
load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")
bad_samples <- which(rownames(pdSpan) %in% c('216', '218', '219'))
pdSpan[bad_samples, ]
if(nrow(pdSpan) == 487) pdSpan <- pdSpan[-bad_samples, ]
stopifnot(nrow(pdSpan) == 484)

files <- pdSpan$wig
names(files) <- pdSpan$lab



## Load the coverage information without filtering
chrnums <- c(1:22, 'X', 'Y')

fullCov <- fullCoverage(files = files, chrs = chrnums, mc.cores = opt$mcores, fileStyle = opt$fileStyle, outputs = 'auto')

message(paste(Sys.time(), 'Saving the full (unfiltered) coverage data'))
save(fullCov, file='fullCov.Rdata')

rm(fullCov)

fullCov_files <- as.list(dir(pattern = 'chr'))
names(fullCov_files) <- sapply(fullCov_files, function(x) gsub('CovInfo.Rdata', '', x))

## Filter the data and save it by chr
myFilt <- function(chr, rawData_file, cutoff, totalMapped = NULL, targetSize = 80e6) {
    library('derfinder')
    
    ## Load raw data
    message(paste(Sys.time(), 'Loading raw file', rawData_file, 'for', chr))
    load(rawData_file)
    rawData <- get(paste0(chr, 'CovInfo'))
    
	## Filter the data
    message(paste(Sys.time(), 'Filtering chromosome', chr))
	res <- filterData(data = rawData$coverage, cutoff = cutoff, index = NULL,
        totalMapped = totalMapped, targetSize = targetSize)
	
	## Save it in a unified name format
	varname <- paste0(chr, 'CovInfo')
	assign(varname, res)
	output <- paste0(varname, '-filtered.Rdata')
	
	## Save the filtered data
	save(list = varname, file = output, compress='gzip')
	
	## Finish
	return(invisible(NULL))
}

message(paste(Sys.time(), 'Filtering and saving the data with cutoff', opt$cutoff))
filteredCov <- bpmapply(myFilt, names(fullCov_files), fullCov_files, BPPARAM = SnowParam(workers = opt$mcores), MoreArgs = list(cutoff = opt$cutoff, totalMapped = totalMapped, targetSize = targetSize))

## Done!
proc.time()
options(width = 120)
session_info()
