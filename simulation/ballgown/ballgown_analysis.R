library('ballgown')
library('GenomicRanges')
library('devtools')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'replicate', 'r', 1, 'integer', 'Replicate number. Either 1, 2 or 3.',
    'complete', 'c', 1, 'character', "'yes' or 'no' for whether the GTF was complete or not",
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

## Perform statistical tests: transcript-level
stat_results <- stattest(gown = bg, feature = 'trans', meas = 'cov', covariate = 'group')
stat_results$pbonf <- p.adjust(stat_results$pval, 'bonferroni')

## Create nice GRangesList: transcript-level
bgres <- structure(bg)$trans
mcols(bgres) <- cbind(mcols(bgres), stat_results)

## Save results
save(stat_results, file = paste0('stat-trans-R', opt$replicate, '-', ifelse(opt$complete == 'yes', 'comp', 'inc'), '.Rdata'))
save(bgres, file = paste0('bgres-trans-R', opt$replicate, '-', ifelse(opt$complete == 'yes', 'comp', 'inc'), '.Rdata'))

## Clean up
rm(stat_results, bgres)

## Perform statistical tests: exon-level
exon_cov  <- eexpr(bg, 'cov')
stat_results <- stattest(gowntable = exon_cov, feature = 'exon', pData = pData(bg), mod = model.matrix(~ pData(bg)$group), mod0 = model.matrix( ~ 0 + rep(1, nrow(pData(bg)))))
stat_results$pbonf <- p.adjust(stat_results$pval, 'bonferroni')

## Create nice GRanges (exon-level)
bgres <- structure(bg)$exon
mcols(bgres) <- cbind(mcols(bgres), stat_results)

## Save results
save(stat_results, file = paste0('stat-exon-R', opt$replicate, '-', ifelse(opt$complete == 'yes', 'comp', 'inc'), '.Rdata'))
save(bgres, file = paste0('bgres-exon-R', opt$replicate, '-', ifelse(opt$complete == 'yes', 'comp', 'inc'), '.Rdata'))


## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
