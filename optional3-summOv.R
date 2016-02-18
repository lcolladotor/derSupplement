## Following http://bioconductor.org/packages/release/data/experiment/vignettes/parathyroidSE/inst/doc/parathyroidSE.pdf

## Setup
library("Rsamtools")
library("derfinder")
library("GenomicRanges")

## Load pre-computed exonic parts from Ensembl annotation
# load("/dcs01/ajaffe/Brain/derRuns/derMisc/ensemblExons/exonicParts.Rdata")

## For comparability, use the exons from genomic state 
load("/home/epi/ajaffe/Lieber/Projects/RNAseq/derannotator/rdas/GenomicState.Hsapiens.ensembl.GRCh37.p11.rda")
exonicParts <- GenomicState.Hsapiens.ensembl.GRCh37.p11$fullGenome[ GenomicState.Hsapiens.ensembl.GRCh37.p11$fullGenome$theRegion == "exon"]

## Make bamFileList
files <- rawFiles(datadir=datadir, samplepatt="out$",
    fileterm="accepted_hits.bam")
names(files) <- gsub('_out', '', names(files))
bai <- paste0(files, ".bai")
bList <- BamFileList(files, bai)

## Compute the overlaps
message(paste(Sys.time(), "summarizeOverlaps: Running summarizeOverlaps()"))
summOverlaps <- summarizeOverlaps(exonicParts, bList, mode="Union",
    singleEnd=TRUE, ignore.strand=TRUE, inter.feature=FALSE, mc.cores=cores)

## Finish
message(paste(Sys.time(), "summarizeOverlaps: Saving summOverlaps"))
save(summOverlaps, file="summOverlaps.Rdata")

proc.time()
sessionInfo()
