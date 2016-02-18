## Adapted from /dcs01/ajaffe/Brain/derRuns/derSoftware/figure1/figure1.R

library('derfinder')
library('derfinderHelper')
library('GenomicRanges')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('RColorBrewer')
library('scales')
library('GenomeInfoDb')
library("GenomicFeatures")

## Define paths
mainPath <- '/dcs01/ajaffe/Brain/derRuns/derSupplement'
covPath <- file.path(mainPath, 'brainspan/CoverageInfo')
resPath <- file.path(mainPath, 'brainspan/derAnalysis/run4-v1.0.10')
dataPath <- '/nexsan2/disk3/ajaffe/BrainSpan/RNAseq/bigwig'

## Load data
load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")
files <- pdSpan$wig
names(files) <- pdSpan$lab

load(file.path(mainPath, 'brainspan', 'regionMatrix', 'regionMat-cut0.25.Rdata'))
load(file.path(resPath, 'groupInfo.Rdata'))
load(file.path(resPath, 'models.Rdata'))

## Use same names as Jaffe
rename <- data.frame(ori = c('Neo.F', 'Neo.A', 'notNeo.F', 'notNeo.A', 'CBC.F', 'CBC.A'), new = c('NCX.F', 'NCX.P', 'NonNCX.F', 'NonNCX.P', 'CBC.F', 'CBC.P'))
levels(groupInfo) <- sapply(levels(groupInfo), function(x) { rename$new[rename$ori == x]})

## Some options
pad <- 30
scalefac <- 1

## Selected region
selected <- range(regionMat$chr5$regions[subjectHits(findOverlaps(GRanges('chr5', IRanges(161110000, 161129000), '*'), regionMat$chr5$regions))])
selected <- resize(selected, width(selected) + 2 * pad, fix = 'center')

## Load coverage
chr <- as.character(seqnames(selected))
chr
cov <- loadCoverage(files = files, which = selected, chr = chr)


## Bases
pos <- start(selected):end(selected)

## Log2 transform coverage
cov.log <- cov$coverage[pos, ]
for(i in seq_len(ncol(cov.log))) {
    cov.log[[i]] <- log2(cov.log[[i]] + scalefac)
}

## Calculate F-stats
fstats <- fstats.apply(data = cov.log, mod = models$mod, mod0 = models$mod0, scalefac = scalefac)
fstats.num <- as.numeric(fstats)
summary(fstats)

## Misc
covDat <- as.data.frame(cov$coverage[pos, ])
covDat.log <- as.data.frame(cov.log)


## F-stat panel
pdf("fstat_panel.pdf", h= 6,w=14)
plot(fstats.num ~ pos, type="l", xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8)
cutoff=2.86420435076022
abline(h=cutoff, lty=2)

fstat.pl <- brewer.pal(3, "Greys")
sl <- slice(fstats.num, lower = cutoff)
for(i in seq(along=sl)) {
	Ind = start(sl)[i]:end(sl)[i]
	polygon(x = c(pos[Ind], rev(pos[Ind])),
		y = c(fstats.num[Ind], rep(cutoff, length(Ind))),
		col = fstat.pl[3], density =60)
}
dev.off()


## Calculate overall mean
mean.ov <- log2(rowMeans(covDat) + scalefac)
y.axis <- 0:7/10

## Mean panel
pdf("mean_panel.pdf", h= 6,w=14)
plot(mean.ov ~ pos, type="l", xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8, yaxt="n")
axis(2, at = log2(y.axis + scalefac), labels = y.axis, cex.axis = 1.2)
mean.cutoff <- log2(0.25 + scalefac)
abline(h= mean.cutoff, lty=2)

mean.sl <- slice(mean.ov, lower = mean.cutoff)
pl <- brewer.pal(length(mean.sl), "Paired")
palette(pl)
for(i in seq(along = mean.sl)) {
	Ind = start(mean.sl)[i]:end(mean.sl)[i]
	polygon(x = c(pos[Ind], rev(pos[Ind])),
		y = c(mean.ov[Ind], rep(mean.cutoff, length(Ind))),
		col = i, density =60)
}
dev.off()




## coverage panel
y.axis.sample <- c(0, 0.5, 2^(0:6))
group.pl <- brewer.pal(6, "Dark2")


pdf("fullCov_panel.pdf", h= 6,w=14)

sample.pl <- mapply(function(col, n) {
    ## Need 1/3 of lines for full saturation
    alpha(col, 7 / n)
}, group.pl, table(groupInfo))

palette(sample.pl)
matplot(pos, covDat.log, yaxt="n",
	col=as.numeric(groupInfo), lty=1, type="l",
	xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8)
axis(2, at = log2(y.axis.sample + scalefac), labels = y.axis.sample, cex.axis = 1.3)
#m = max(covDat.log)
m <- log2(32 + scalefac)
for(i in seq(along=mean.sl)) {
	Ind = start(mean.sl)[i]:end(mean.sl)[i]
	rect(xleft=min(pos[Ind]), xright = max(pos[Ind]),
		ybot = 0, ytop = m, col=pl[i], density=10)
}
palette(group.pl)
legend("top", levels(groupInfo), col=seq_len(length(levels(groupInfo))), cex=1.4,pch=15, ncol = 6, bty = 'n')
dev.off()


## annotate
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
ensemblAnno <- annotateRegions(selected,
    GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome)
ensemblCount <- ensemblAnno$countTable

### gene plot
a = as.data.frame(ensemblAnno$annotationList)
Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))

pdf("gene_anno.pdf", h=3,w=14)
plot(0,0, type="n", xlim=range(pos),ylim=c(-1.5,1.5),yaxt="n",ylab="",
	xlab=paste("Chromosome", mapSeqlevels(chr, 'NCBI')),	cex.axis = 1.5, cex.lab =1.8)
axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
abline(h=0,lty=3)
for (j in seq_len(nrow(a))) {
	polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]),
	  Strand[j]/2 + c(-0.3, -0.3, 0.3, 0.3) * Lwd[j],
	  col = Col[j])
}
e <- a[a$theRegion == "exon", ]
s2 <- Strand[a$theRegion == "exon"]
g = unlist(e$symbol)
g[is.na(g)] = ""
if (length(g) > 0) {
	text(x = e$start + e$width/2, y = s2 * 0.8, g,
	  font = 2, pos = s2 + 2,
      cex = c(1.2, 0.01, 0.5, 0.5, 0.5, 0.01, 1.2, 1.2, 1.2, 0.01, 0.01))
}
dev.off()

#### extra tx info

txdb <- loadDb("/home/epi/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12/inst/extdata/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12.sqlite")
txdb <- keepSeqlevels(txdb, mapSeqlevels(chr, 'NCBI'))
seqlevelsStyle(txdb) <- 'UCSC'
tx=exonsBy(txdb)
eList = tx[subjectHits(findOverlaps(selected, tx) )]

e.strand <- unlist(unique(strand(eList)))
e.n.neg <- sum(e.strand == '-')
e.n.pos <- sum(e.strand == '+')
ylim <- c(-1 * e.n.neg + ifelse(e.n.neg > 0, -0.5, 0.5), e.n.pos + 0.5)

pdf("trans_anno.pdf", h=4.5,w=14)
plot(0,0, type="n", xlim=range(pos), ylim=ylim,
	yaxt="n",ylab="", xlab=paste("Chromosome", mapSeqlevels(chr, 'NCBI'), '(161.1 mb)'), xaxt='n', cex.lab = 1.8)
axis(1, at = c(161115000, 161120000, 161125000, 161130000), labels = c('+15k', '+20k', '+25k', '+30k'), cex.axis = 1.5)
axis(2, c(- ifelse(e.n.neg, median(seq_len(e.n.neg)), NA), ifelse(e.n.pos, median(seq_len(e.n.pos)), NA)), c(ifelse(e.n.neg, '-', NA), ifelse(e.n.pos, "+", NA)), tick=FALSE,las=1,cex.axis = 3)
abline(h=0,lty=3)
for(i in seq(along=eList)) {
	a = as.data.frame(eList[[i]])
    i.strand <- sum(e.strand[ seq_len(length(e.strand)) <= i] == e.strand[i]) * ifelse(e.strand[i] == "+", 1, -1)
	for (j in seq_len(nrow(a))) {
		polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]), 
			c(i.strand - 0.25, i.strand -0.25, i.strand +0.25, i.strand +0.25), col="blue")
	}
	
	int = gaps(eList[[i]])
	int = int[seqnames(int) == unique(seqnames(eList[[i]]))]
    int <- int[ end(int) < seqlengths(int) & start(int) > 1]
	end(int) = end(int)+1
	int = as.data.frame(int[start(int) != 1])
	
    
	for (j in seq_len(nrow(int))) {
		polygon(c(int$start[j], int$end[j], int$end[j], int$start[j]), 
			c(i.strand - 0.15, i.strand -0.15, i.strand + 0.15, i.strand +0.15), col="lightblue")
	}
    
}
dev.off()


