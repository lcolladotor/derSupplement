## Get exon coverage table for UCSC or Ensembl annotation

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')

## Specify parameters
spec <- matrix(c(
    'experiment', 'e', 1, 'character', 'Experiment. Only brainspan',
    'annotation', 'a', 1, 'character', 'Annotation to use. Either ensembl or ucsc',
    'readlen', 'r', 1, 'integer', 'Read length',
    'mc.cores', 'c', 1, 'integer', 'Number of cores to use',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Check experiment input
stopifnot(opt$experiment %in% c('brainspan'))



## Setup
message(paste(Sys.time(), 'loading coverage data'))
load("../CoverageInfo/fullCov.Rdata")

message(paste(Sys.time(), 'loading annotation data'))
if(opt$annotation == 'ensembl') {
    load('/users/ajaffe/Lieber/Projects/RNAseq/derannotator/rdas/GenomicState.Hsapiens.ensembl.GRCh37.p11.rda')
    anno <- GenomicState.Hsapiens.ensembl.GRCh37.p11$fullGenome

} else if (opt$annotation == 'ucsc') {
    load('/users/ajaffe/Lieber/Projects/RNAseq/derannotator/rdas/GenomicState.Hsapiens.UCSC.hg19.knownGene.rda')
    anno <- GenomicState.Hsapiens.UCSC.hg19.knownGene$fullGenome
}


## get table
message(paste(Sys.time(), 'running coverageToExon'))
covToEx <- coverageToExon(fullCov, anno, L=opt$readlen, strandCores = 1, mc.cores = opt$mc.cores)

## Save results
message(paste(Sys.time(), paste0("saving covToEx-", opt$annotation, ".Rdata")))
save(covToEx, file=paste0("covToEx-", opt$annotation, ".Rdata"))

## Done
proc.time()
sessionInfo()
