## Setup
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('derfinder')
library('Rsamtools')
library('GenomicAlignments')
library('parallel')
options(mc.cores=24)

## Exons by gene
ex <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = 'gene')

## Make bamFileList
files <- rawFiles(datadir='/dcs01/ajaffe/Hippo/TopHat', samplepatt="out$",
    fileterm="accepted_hits.bam")
names(files) <- gsub('_out', '', names(files))
bai <- paste0(files, ".bai")
bList <- BamFileList(files, bai)

## Compute the overlaps
message(paste(Sys.time(), "summarizeOverlaps: Running summarizeOverlaps()"))
summOv <- summarizeOverlaps(ex, bList, mode="Union",
    singleEnd=TRUE, ignore.strand=TRUE)

## Finish
message(paste(Sys.time(), "summarizeOverlaps: Saving summOverlaps"))
save(summOv, file="summOv.Rdata")

proc.time()
options(width = 120)
devtools::session_info()
