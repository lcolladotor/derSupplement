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
if(opt$datadir == '/dcs01/ajaffe/Snyder/RNAseq/TopHat') {
    files <- rawFiles(datadir=opt$datadir, samplepatt=opt$pattern, fileterm="accepted_hits.bam")
    names(files) <- gsub('_out', '', names(files))
} else if(opt$datadir == '/dcs01/ajaffe/Hippo/TopHat') {
    files <- rawFiles(datadir=opt$datadir, samplepatt=opt$pattern)
    
    ## In some cases, you might want to modify the names of the files object
    ## These names specify the column names used in the DataFrame objects.
    ## For example, they could end with _out
    names(files) <- gsub('_out', '', names(files))
} else if (opt$datadir == '/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/tophat') {
    files <- rawFiles(datadir = opt$datadir, samplepatt = 'SRR')
    load('/dcl01/lieber/ajaffe/PublicData/SRA_GTEX/gtexPd.Rdata')
    
    if(file.exists('/dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/mappedInfo.Rdata')) {
        load('/dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/mappedInfo.Rdata')
    } else {
        ## Taken from /home/epi/ajaffe/Lieber/lieber_functions_aj.R
        getTotalMapped = function(bamFile, mc.cores=1, returnM = TRUE) {
        	thecall = paste("samtools idxstats", bamFile)
        	tmp = parallel::mclapply(thecall, function(x) {
        		cat(".")
        		xx = system(x,intern=TRUE)
        		xx = do.call("rbind", strsplit(xx, "\t"))
        		d = data.frame(chr=xx[,1], L=xx[,2], mapped = xx[,3],
        			stringsAsFactors=FALSE)
        		d
        	},mc.cores=mc.cores)

        	out = list(totalMapped = sapply(tmp, function(x) sum(as.numeric(x$mapped[x$chr %in% paste0("chr", c(1:22,"X","Y"))]))),
        		mitoMapped = sapply(tmp, function(x) as.numeric(x$mapped[x$chr=="chrM"])))
        	return(out)
        }

        libSize <- getTotalMapped(files, mc.cores = opt$mcores)
        mappedInfo <- data.frame(sample = names(files), file = files,
            totalMapped = libSize$totalMapped, mitoMapped = libSize$mitoMapped)
        rownames(mappedInfo) <- NULL
        ## Display the info
        print(mappedInfo)
        save(mappedInfo, file = '/dcs01/ajaffe/Brain/derRuns/derSupplement/gtex/mappedInfo.Rdata')
    }
    totalMapped <- mappedInfo$totalMapped
} else {
    load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")
    files <- pdSpan$wig
    names(files) <- pdSpan$lab
}


## Load the coverage information without filtering
chrnums <- c(1:22, 'X', 'Y')

fullCov <- fullCoverage(files = files, chrs = chrnums, mc.cores = opt$mcores, fileStyle = opt$fileStyle)

message(paste(Sys.time(), 'Saving the full (unfiltered) coverage data'))
save(fullCov, file='fullCov.Rdata')

## Filter the data and save it by chr
myFilt <- function(chr, rawData, cutoff, totalMapped = NULL, targetSize = 80e6) {
    library('derfinder')
    message(paste(Sys.time(), 'Filtering chromosome', chr))
    
	## Filter the data
	res <- filterData(data = rawData, cutoff = cutoff, index = NULL,
        totalMapped = totalMapped, targetSize = targetSize)
	
	## Save it in a unified name format
	varname <- paste0(chr, 'CovInfo')
	assign(varname, res)
	output <- paste0(varname, '.Rdata')
	
	## Save the filtered data
	save(list = varname, file = output, compress='gzip')
	
	## Finish
	return(invisible(NULL))
}

message(paste(Sys.time(), 'Filtering and saving the data with cutoff', opt$cutoff))
filteredCov <- bpmapply(myFilt, names(fullCov), fullCov, BPPARAM = SnowParam(opt$mcores, outfile = Sys.getenv('SGE_STDERR_PATH')), MoreArgs = list(cutoff = opt$cutoff, totalMapped = totalMapped, targetSize = targetSize))

## Done!
proc.time()
options(width = 120)
session_info()
