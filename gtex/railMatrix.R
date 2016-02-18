library('derfinder')
library('rtracklayer')
library('BiocParallel')
library('devtools')

## Options
cutoff <- 5L
L <- 76L
maxClusterGap <- 3000L
targetSize <- 40e6

message(Sys.time())
timeinfo <- NULL
timeinfo <- c(timeinfo, list(Sys.time()));

## Identify the summary files
sfiles <- rawFiles(datadir = '/dcl01/leek/data/cwilks/gtex36/coverage_bigwigs', samplepatt = 'mean', fileterm = NULL)
sfiles <- sfiles[!grepl('unique', sfiles)]
names(sfiles) <- gsub('.bw', '', names(sfiles))
sfiles <- BigWigFileList(sfiles)

## Reorder to match chr order
chrs <- paste0('chr', c(1:22, 'X', 'Y'))
names(chrs) <- chrs
sfiles <- sfiles[match(chrs, gsub('mean\\.', '', names(sfiles)))]

## Identify the data files
bws <- rawFiles(datadir = '/dcl01/leek/data/cwilks/gtex36/coverage_bigwigs', samplepatt = 'bw', fileterm = NULL)
bws <- bws[!grepl('mean|median|unique', bws)]
names(bws) <- gsub('.bw', '', names(bws))
bws <- BigWigFileList(bws)

## Find counts
counts <- read.table('/dcl01/leek/data/cwilks/gtex36/cross_sample_results/counts.tsv.gz', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
counts$totalMapped <- as.integer(sapply(strsplit(counts$total.mapped.reads, ','), '[[', 1))
mapped <- counts$totalMapped[match(names(bws), counts$X)]
names(mapped) <- names(bws)

## Check distribution
summary(mapped)
summary(mapped) / 1e6

## Load mean coverage, find regions and then get the coverage matrix
railMatrix <- function(chr, summaryFile, sampleFiles, L, maxClusterGap, cutoff, totalMapped, targetSize) {
    meanCov <- loadCoverage(files = summaryFile, chr = chr)
    regs <- findRegions(position = Rle(TRUE, length(meanCov$coverage[[1]])), fstats = meanCov$coverage[[1]], chr = chr, maxClusterGap = maxClusterGap, cutoff = cutoff)
    
    ## Format appropriately
    names(regs) <- seq_len(length(regs))
    
    ## Set the length
    seqlengths(regs) <- length(meanCov$coverage[[1]])
    
    ## Get the region coverage matrix
    regionCov <- getRegionCoverage(regions = regs, files = sampleFiles, which = regs, totalMapped = totalMapped, targetSize = targetSize, mc.cores = 1L, protectWhich = 0)
    covMat <- lapply(regionCov, colSums)
    covMat <- do.call(rbind, covMat) / L
    
    ## Finish
    res <- list(regions = regs, coverageMatrix = covMat)
    return(res)
}

## Run it all in parallel
regionMat <- bpmapply(railMatrix, chrs, sfiles, SIMPLIFY = FALSE, MoreArgs = list(sampleFiles = bws, L = L, maxClusterGap = maxClusterGap, cutoff = cutoff, totalMapped = mapped, targetSize = targetSize), BPPARAM = MulticoreParam(workers = 10, outfile = Sys.getenv('SGE_STDERR_PATH')))

## Test
#regionMat <- bpmapply(railMatrix, chrs[23:24], sfiles[23:24], SIMPLIFY = FALSE, MoreArgs = list(sampleFiles = bws, L = L, maxClusterGap = maxClusterGap, cutoff = cutoff, totalMapped = mapped, targetSize = targetSize), BPPARAM = SerialParam())

timeinfo <- c(timeinfo, list(Sys.time()))

## Save results
save(regionMat, file=paste0('regionMat-cut', cutoff, '.Rdata'))
timeinfo <- c(timeinfo, list(Sys.time()))

## Save time information
save(timeinfo, file=paste0('timeinfo-cut', cutoff, '.Rdata'))

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
