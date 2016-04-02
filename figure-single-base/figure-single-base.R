## Adapted from /home/epi/ajaffe/Lieber/Projects/derfinderPaper/figure1_stem.R

## Usage:
# qrsh -l mem_free=80G,h_vmem=90G
# module load R/3.3
# mkdir -p logs
# Rscript figure-single-base.R > logs/figure-single-base_log.txt 2>&1


library('derfinder')
library('derfinderHelper')
library('derfinderPlot')
library('GenomicRanges')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')
library('RColorBrewer')
library('scales')
library('GenomeInfoDb')
library("GenomicFeatures")

## Define paths
#mainPath <- '/dcl01/lieber/ajaffe/Brain/derRuns/derSupplement/'
mainPath <- '/dcl01/lieber/ajaffe/Brain/derRuns/derSoftware/'
covPath <- file.path(mainPath, 'brainspan/CoverageInfo/')
resPath <- file.path(mainPath, 'brainspan/derAnalysis/run4-v1.0.10')
dataPath <- '/nexsan2/disk3/ajaffe/BrainSpan/RNAseq/bigwig/'

## Load data
load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")

## Remove bad samples
bad_samples <- which(rownames(pdSpan) %in% c('216', '218', '219'))
pdSpan[bad_samples, ]
if(nrow(pdSpan) == 487) pdSpan <- pdSpan[-bad_samples, ]
stopifnot(nrow(pdSpan) == 484)

files <- pdSpan$wig
names(files) <- pdSpan$lab

load(file.path(resPath, 'fullRegions.Rdata'))
ders <- fullRegions

load(file.path(resPath, 'groupInfo.Rdata'))
if(length(groupInfo) == 487) groupInfo <- groupInfo[-bad_samples]
stopifnot(length(groupInfo) == 484)

load(file.path(resPath, 'models.Rdata'))
stopifnot(nrow(models$mod) == nrow(models$mod0))
if(unique(sapply(models, nrow)) == 487) {
    models$mod <- models$mod[-bad_samples, ]
    models$mod0 <- matrix(models$mod0[-bad_samples, ], ncol = 1)
}
stopifnot(nrow(models$mod) == 484)

## Use same names as Jaffe
rename <- data.frame(ori = c('Neo.F', 'Neo.A', 'notNeo.F', 'notNeo.A', 'CBC.F', 'CBC.A'), new = c('NCX.F', 'NCX.P', 'NonNCX.F', 'NonNCX.P', 'CBC.F', 'CBC.P'))
levels(groupInfo) <- sapply(levels(groupInfo), function(x) { rename$new[rename$ori == x]})
# levels(groupInfo) <- gsub('notNeo', '-Neo', gsub('\\.', '-', levels(groupInfo)))


## Some options
scalefac <- 1
pad <- 600

## Get best cluster of ders
cluster <- data.frame(area = ders$area,
    clusterChr = paste0(as.integer(ders$cluster),
    chr = as.character(seqnames(ders))))
regionClustAreas <- tapply(cluster$area, cluster$clusterChr, sum)
bestArea <- sapply(names(head(sort(regionClustAreas, decreasing=TRUE),
    20)), function(y) { which(cluster$clusterChr == y)[[1]]})
    
## Explore some clusters
skip <- TRUE
if(!skip) {
    library('epivizr')
    mgr <- startEpiviz()
    ders_dev <- mgr$addDevice(ders[!as.logical(ders$significantFWER)], "DERs no sig FWER")
    ders_sig_dev <- mgr$addDevice(ders[as.logical(ders$significantFWER)], "DERs sig FWER")

    ## Explore some clusters
    for(i in seq_len(20)) {
        i.cluster <- bestArea[i]
        r.cluster <- range(ders[which(cluster$clusterChr == names(i.cluster))])
        print(i)
        print(width(r.cluster))
        readline('Continue?')
        mgr$navigate(as.character(seqnames(ders[i.cluster])), start(r.cluster) - pad, end(r.cluster) + pad)
    }
}



## Selected DER
s <- bestArea[16]
s.cluster <- range(ders[which(cluster$clusterChr == names(s))])
selected <- resize(s.cluster, width(s.cluster) + 2 * pad, fix = 'center')

## Load coverage and annotation
chr <- as.character(seqnames(selected))
chr
cov <- loadCoverage(files = files, which = selected, chr = chr)
load(file.path(resPath, chr, 'annotation.Rdata'))

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


## Plot cluster
p.cluster <- plotCluster(idx = s, regions = ders, annotation = annotation, titleUse = 'fwer', groupInfo = groupInfo, coverageInfo = cov$coverage, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, maxExtend = pad)
pdf(paste0('plotCluster-', s, '.pdf'))
print(p.cluster)
dev.off()
save(p.cluster, file = paste0('plotCluster-', s, '.Rdata'))

## Panel 1
bpInd <- start(ders[s]) + 100:103
ind <- start(ders[s]) - start(selected) + 100:103
covDat <- as.data.frame(cov$coverage[pos, ])
covDat.log <- as.data.frame(cov.log)
ylim <- log2(c(0, ceiling(max(covDat[ind + 1, ]))) + scalefac)
y.axis <- c(0, 0.5, 2^(0:3))

group.pl <- brewer.pal(6, "Dark2")

for(i in seq(along=ind)) {
	pdf(paste0("cov_part", i,".pdf"), h=6, w=4.75)
	palette(group.pl)

	y <- as.numeric(log2(covDat[ind[i] + 1,] + scalefac))
	boxplot(y ~ groupInfo, outline=FALSE, yaxt="n", main="", ylim = ylim,
        xaxt="n")
	mtext(paste0(chr, ":", bpInd[i]), line=0.5, cex=1.6)
	axis(2, at = log2(y.axis + scalefac), labels = y.axis, cex.axis = 1.3)
	points(y ~ jitter(as.numeric(groupInfo), amount=0.15),
		pch = 21, bg = as.numeric(groupInfo), cex=0.8)
    legend("topright", paste0("F=", round(fstats.num[ind[i] + 1], 2)), cex=1.3)
    if(i <= 3) {
        par(xpd = TRUE)
        legend("bottom", levels(groupInfo)[c(2 * i - 1, 2 * i)], col = group.pl[c(2 * i - 1, 2 * i)], cex=1.7, inset = -0.13, ncol = 2, bty = 'n', pch = 16)
    }	
	dev.off()
}


## panel 2
pdf("fstat_panel.pdf", h= 6,w=14)
plot(fstats.num ~ pos, type="l", xlab=chr, ylab="", cex.axis=1.4, cex.lab=1.8)
cutoff=2.86420435076022
abline(h=cutoff, lty=2)

pl = brewer.pal(9, "Paired")
palette(pl)
sl = slice(fstats.num, lower = cutoff)
for(i in seq(along=sl)) {
	Ind = start(sl)[i]:end(sl)[i]
	polygon(x = c(pos[Ind], rev(pos[Ind])),
		y = c(fstats.num[Ind], rep(cutoff, length(Ind))),
		col = i, density =60)
}
abline(v=range(bpInd), col="red")
dev.off()


## coverage panel
y.axis.sample <- c(y.axis, 2^4)
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
m <- log2(8 + scalefac)
for(i in seq(along=sl)) {
	Ind = start(sl)[i]:end(sl)[i]
	rect(xleft=min(pos[Ind]), xright = max(pos[Ind]),
		ybot = 0, ytop =m, col=pl[i], density=10)
}
palette(group.pl)
legend("topright", levels(groupInfo), col=seq_len(length(levels(groupInfo))), cex=1.4,pch=15, ncol = 6)
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
	  font = 2, pos = s2 + 2, cex = rep(c(1.2, 0.5, 1.2, 0.01), c(3, 3, 1, 1)))
}
dev.off()

#### extra tx info

txdb <- loadDb("/home/epi/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12/inst/extdata/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12.sqlite")
txdb <- keepSeqlevels(txdb, mapSeqlevels(chr, 'NCBI'))
seqlevelsStyle(txdb) <- 'UCSC'
tx=exonsBy(txdb)
eList = tx[subjectHits(findOverlaps(selected, tx) )]

pdf("trans_anno.pdf", h=2.5,w=14)
plot(0,0, type="n", xlim=range(pos),ylim=c(0.5,length(eList)+0.5),
	yaxt="n",ylab="", xlab=paste("Chromosome", mapSeqlevels(chr, 'NCBI')), cex.axis = 1.5, cex.lab =1.8)
for(i in seq(along=eList)) {
	a = as.data.frame(eList[[i]])
	for (j in seq_len(nrow(a))) {
		polygon(c(a$start[j], a$end[j], a$end[j], a$start[j]), 
			c(i-0.25, i-0.25, i+0.25, i+0.25), col="blue")
	}
	
	int = gaps(eList[[i]])
	int = int[seqnames(int) == unique(seqnames(eList[[i]]))]
    int <- int[ end(int) < seqlengths(int) & start(int) > 1]
	end(int) = end(int)+1
	int = as.data.frame(int[start(int) != 1])
	
	for (j in seq_len(nrow(int))) {
		polygon(c(int$start[j], int$end[j], int$end[j], int$start[j]), 
			c(i-0.15, i-0.15, i+0.15, i+0.15), col="lightblue")
	}
}
dev.off()

## Reproducibility info
library('devtools')
options(width = 120)
session_info()
Sys.time()
proc.time()
