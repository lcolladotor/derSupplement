# Setup
library("GenomicRanges")
library("rtracklayer")

load('/home/epi/ajaffe/Lieber/Projects/RNAseq/derannotator/rdas/GenomicState.Hsapiens.ensembl.GRCh37.p11.rda')
load('/home/epi/ajaffe/Lieber/Projects/RNAseq/derannotator/rdas/GenomicState.Hsapiens.UCSC.hg19.knownGene.rda')

makeGFF <- function(exonicParts, file) {
	message(paste(Sys.time(), "makeGFF: Saving", file))
	export.gff2(exonicParts, file)
}

makeGFF(GenomicState.Hsapiens.ensembl.GRCh37.p11$fullGenome[ GenomicState.Hsapiens.ensembl.GRCh37.p11$fullGenome$theRegion == "exon"], "GenomicState.Hsapiens.ensembl.GRCh37.p11.exons.gff")
makeGFF(GenomicState.Hsapiens.UCSC.hg19.knownGene$fullGenome[ GenomicState.Hsapiens.UCSC.hg19.knownGene$fullGenome$theRegion == "exon"], "GenomicState.Hsapiens.UCSC.hg19.knownGene.exons.gff")

proc.time()
sessionInfo()
