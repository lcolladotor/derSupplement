## Original script: /home/epi/ajaffe/Lieber/Projects/derfinderPaper/analyze_gtex.R
##
## Usage:
# mkdir -p logs
# Rscript analyze_gtex.R > logs/analyze_gtex_log.txt 2>&1
library('derfinder')
library('derfinderPlot')
library('GenomicRanges')
library('rafalib')
library('GenomeInfoDb')
library('devtools')
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

## Create dir for saving rdas and plots to preserve structure from the original script
dir.create('rdas', recursive = TRUE, showWarnings = FALSE)
dir.create('plots', recursive = TRUE, showWarnings = FALSE)

# get the f statistic from 2 lmFit objects
getF = function(fit, fit0, theData) {
	
	rss1 = rowSums((fitted(fit)-theData)^2)
	df1 = ncol(fit$coef)
	rss0 = rowSums((fitted(fit0)-theData)^2)
	df0 = ncol(fit0$coef)

	fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
	f_pval = pf(fstat, df1-1, ncol(theData)-df1,lower.tail=FALSE)
	fout = cbind(fstat,df1-1,ncol(theData)-df1,f_pval)
	colnames(fout)[2:3] = c("df1","df0")
	fout = data.frame(fout)
	return(fout)
}

# load processed data
load('/dcl01/lieber/ajaffe/derRuns/derSupplement/gtex/regionMat-cut5.Rdata')
load("/dcl01/lieber/ajaffe/derRuns/derSupplement/gtex/gtex_pheno_with_mapped.Rdata")
## Should be the same:
stopifnot(identical(match(colnames(regionMat[[1]]$coverageMatrix), pd2$sra_accession), seq_len(ncol(regionMat[[1]]$coverageMatrix))))
pd2$Tissue <- pd2$SMTS
pd2$SubjectID <- pd2$SUBJID

## Check that we have the same number of samples per tissue
stopifnot(table(pd2$Tissue) - max(table(pd2$Tissue)) == 0)

### extract coverage data
regions = unlist(GRangesList(lapply(regionMat, '[[', 'regions')))
names(regions) = NULL
regionMat = do.call("rbind", lapply(regionMat, '[[', 'coverageMatrix'))
rownames(regionMat) = names(regions) = paste0("er", 1:nrow(regionMat))

### filter out short DERs
keepIndex = width(regions) > 8
regionMat = regionMat[keepIndex,]
regions = regions[keepIndex]

#################
#### analysis ###

# transform and offset
y = log2(regionMat+1)

## width of regions
quantile(width(regions))

### annotate
## Genomic state created by https://github.com/nellore/runs/blob/master/gtex/DER_analysis/coverageMatrix/genomicState/hg38-genomicState.R
load('/dcl01/leek/data/gtex_work/runs/gtex/DER_analysis/coverageMatrix/genomicState/genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5.Rdata')
gs_raw <- genomicState.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5$fullGenome
gs <- renameSeqlevels(gs_raw, paste0('chr', seqlevels(gs_raw)))

## Do the seqlengths match?
stopifnot(max(abs(seqlengths(regions) - seqlengths(gs)[names(seqlengths(regions))])) == 0)

ensemblAnno = annotateRegions(regions,gs)
countTable = ensemblAnno$countTable

pdf(file = 'plots/venn-GRCh38.p5.pdf')
vennRegions(ensemblAnno, main = 'GTEx expressed regions by GRCh38.p5', counts.col = 'blue')
dev.off()

## annotation ####
dim(countTable)
annoClassList = list(strictExonic = 
	which(countTable[,"exon"] > 0 & countTable[,"intron"] == 0 &
		countTable[,"intergenic"] == 0),
	strictIntronic = 
	which(countTable[,"intron"] > 0 & countTable[,"exon"] == 0 &
		countTable[,"intergenic"] == 0),
	strictIntergenic = which(countTable[,"intergenic"] > 0 & countTable[,"exon"] == 0 &
    countTable[,"intron"] == 0),
	exonIntron = which(countTable[,"exon"] > 0 & countTable[,"intron"] > 0 &
		countTable[,"intergenic"] == 0))

annoClassList$All = 1:nrow(regionMat) # add all

## Explore numbers
sapply(annoClassList, length)
100 * sapply(annoClassList, length) / nrow(countTable)
cumsum(100 * sapply(annoClassList, length)[-5] / nrow(countTable))

# width by annotation
t(sapply(annoClassList, function(ii) quantile(width(regions[ii]))))

### PCA ###
pcList = lapply(annoClassList, function(ii) {
	cat(".")
	pc = prcomp(t(y[ii,]))
	pc$rot = NULL # drop rotations
	return(pc) 
})
pcVarMat = sapply(pcList, getPcaVars)
rownames(pcVarMat) = paste0("PC", 1:nrow(pcVarMat))
pc1Mat = sapply(pcList, function(x) x$x[,1])
pc2Mat = sapply(pcList, function(x) x$x[,2])

## plots
ind = c(1:3,5)
pdf(file = 'plots/pca-simple.pdf')
rafalib::mypar(2,2,cex.axis=1)
for(i in ind) boxplot(pc1Mat[,i] ~ pd2$Tissue, 
	main = colnames(pc1Mat)[i],
	ylab=paste0("PC1: ", pcVarMat[1,i], "% of Var Explain"))

for(i in ind) boxplot(pc2Mat[,i] ~ pd2$Tissue, 
	main = colnames(pc2Mat)[i],
	ylab=paste0("PC2: ", pcVarMat[2,i], "% of Var Explain"))
dev.off()
    
    
## Simple plots for PC1 and PC2 with some added color and text
pdf(file = 'plots/pca-plots-gtex.pdf', width = 14, height = 7)
cnames <- c('Strictly exonic ERs', 'Strictly intronic ERs')
rafalib::mypar(1,2,cex.axis=1)
for(i in ind[1:2]) {
    boxplot(pc1Mat[,i] ~ pd2$Tissue, main = cnames[i], ylab = '', cex.axis = 1.5, cex.main = 2)
    text(3, sum(range(pc1Mat[,i])) / 2, labels = paste0("PC1: ", pcVarMat[1,i], "%\nof Var Explain"), col = 'dodgerblue2', cex = 2, font = 2)
}

for(i in ind[1:2]) {
    boxplot(pc2Mat[,i] ~ pd2$Tissue, main = cnames[i], ylab = '', cex.axis = 1.5, cex.main = 2)
    text(3, sum(range(pc2Mat[,i])) / 2, labels = paste0("PC2: ", pcVarMat[2,i], "%\nof Var Explain"), col = 'dodgerblue2', cex = 2, font = 2)
}
dev.off()


library('RColorBrewer')
colors <- brewer.pal(3, 'Set1')
names(colors) <- c('Liver', 'Heart', 'Testis')
titles <- c('strictExonic' = 'Strictly exonic ERs', 'strictIntronic' = 'Strictly intronic ERs', 'strictIntergenic' = 'Strictly intergenic ERs', 'All' = 'All ERs')

pdf(file = 'plots/pca-PC1-vs-PC2_all.pdf')
rafalib::mypar(2,2,cex.axis=1)
for(i in ind) {
    plot(x = pc1Mat[,i], y = pc2Mat[, i], col = colors[pd2$Tissue], pch = 20, xlab = paste0("PC1: ", pcVarMat[1,i], "% of variance explained"), ylab = paste0("PC2: ", pcVarMat[2, i], "% of variance explained"), main = titles[colnames(pc2Mat)[i]])
    if(i == 1)
        legend(0.5, 0.5, names(colors), bty = 'n', lwd = 4, col = colors, cex = 2)
}
dev.off()

pdf(file = 'plots/pca-PC1-vs-PC2.pdf', width = 12, height = 6)
rafalib::mypar(1,2,cex.axis=1, cex.lab = 1.5, mar = c(3, 3, 1.6, 1.1))
for(i in ind) {
    plot(x = pc1Mat[,i], y = pc2Mat[, i], col = colors[pd2$Tissue], pch = 20, xlab = paste0("PC1: ", pcVarMat[1,i], "% of variance explained"), ylab = paste0("PC2: ", pcVarMat[2, i], "% of variance explained"), main = titles[colnames(pc2Mat)[i]], cex = 2, cex.main = 2)
    if(i == 1 | i == 5)
        legend(0.5, 0.5, names(colors), bty = 'n', lwd = 4, col = colors, cex = 2)
}
dev.off()




	
#################
## DE analysis ##
library('limma')
mod = model.matrix(~pd2$Tissue)
mod0 = model.matrix(~1, data=pd2)
fit = lmFit(y, mod)
eb = ebayes(fit)
fit0 = lmFit(y,mod0)
ff = getF(fit, fit0, y)

outStats = data.frame(log2FC_LiverVsHeart = fit$coef[,2],
	log2FC_TestesVsHeart = fit$coef[,3],
	pval_LiverVsHeart = eb$p[,2],
	pval_TestesVsHeart = eb$p[,3])
	
outStats$fstat = ff$fstat
outStats$fPval = ff$f_pval
outStats$fBonf = p.adjust(outStats$fPval, "bonferroni")

sapply(annoClassList, function(ii) sum(outStats$fBonf[ii] < 0.05))
sapply(annoClassList, function(ii) mean(outStats$fBonf[ii] < 0.05))

################################### 
### conditional analysis ##########
## intron ~ tissue + nearExon #####
## intergenic ~ tissue + nearExon #

# extract introns
intronMat = y[annoClassList[["strictIntronic"]],]
intronRegions = regions[annoClassList[["strictIntronic"]]]

# extract exons
exonMat = y[annoClassList[["strictExonic"]],]
exonRegions = regions[annoClassList[["strictExonic"]]]

# get exon nearest to each intron
ooExon = distanceToNearest(intronRegions, exonRegions)
exonMatMatch = exonMat[subjectHits(ooExon),]
exonRegionsMatch = exonRegions[subjectHits(ooExon)]

## Some might not be matching: potential flag for not using the correct annotation!
length(ooExon) == length(intronRegions)
max(queryHits(ooExon)) == length(ooExon)

# PC1 versus distance
pdf(file = "plots/PC1vsDistance.pdf")
plot(x = mcols(ooExon)$distance, y = pcList$strictIntronic$rot[, 1][queryHits(ooExon)], ylab = 'PC1', xlab = 'Distance to nearest exon', cex = 0.5)
reg1 <- lm(pcList$strictIntronic$rot[, 1][queryHits(ooExon)] ~ mcols(ooExon)$distance)
abline(reg1, col = 'orange')
plot(x = log(mcols(ooExon)$distance + 1), y = pcList$strictIntronic$rot[, 1][queryHits(ooExon)], ylab = 'PC1', xlab = 'Distance to nearest exon: log(x + 1)', cex = 0.5)
reg2 <- lm(pcList$strictIntronic$rot[, 1][queryHits(ooExon)] ~ log(mcols(ooExon)$distance + 1))
abline(reg2, col = 'orange')
dev.off()

# conditional regression
outStatsExon = matrix(NA, ncol = 2, nrow = nrow(intronMat))
for(i in 1:nrow(outStatsExon)) {
	if(i %% 1000 == 0) cat(".")
	f = lm(intronMat[i,] ~ 
			pd2$Tissue + exonMatMatch[i,])
	f0 = lm(intronMat[i,] ~ exonMatMatch[i,])
	outStatsExon[i,] = as.numeric(anova(f,f0)[2,5:6])
}
colnames(outStatsExon) = c("Fstat", "pval")
rownames(outStatsExon) = rownames(intronMat)
outStatsExon=as.data.frame(outStatsExon)

outStatsExon$nearExon <- outStatsExon$nearDist <- NA
outStatsExon$nearExon[queryHits(ooExon)] = rownames(exonMatMatch)
outStatsExon$nearDist[queryHits(ooExon)] = mcols(ooExon)$distance

## get gene symbol
library('GenomicFeatures')

## Get data from Biomart
#system.time(xx <- makeTxDbPackageFromBiomart(version = '0.99', maintainer = 'Leonardo Collado-Torres <lcollado@jhu.edu>', author = 'Leonardo Collado-Torres <lcollado@jhu.edu>', destDir = '~/'))

## Load info
sql_file <- "/home/bst/student/lcollado/TxDb.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5/inst/extdata/TxDb.Hsapiens.BioMart.ENSEMBLMARTENSEMBL.GRCh38.p5.sqlite"
TranscriptDb <- loadDb(sql_file)

## Fix seqlevels
seqlevels(TranscriptDb,force=TRUE) = c(1:22,"X","Y","MT")
seqlevels(TranscriptDb) = paste0("chr", c(1:22,"X","Y","M"))
ensGene = genes(TranscriptDb)

library('biomaRt')
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 83, hg38
	dataset="hsapiens_gene_ensembl",
	host="dec2015.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), 
	values=names(ensGene), mart=ensembl)
ensGene$Symbol = sym$hgnc_symbol[match(names(ensGene), sym$ensembl_gene_id)]
ensGene = ensGene[!grepl("^MIR[0-9]", ensGene$Symbol)] # drop mirs

## match
oo1 = findOverlaps(intronRegions, ensGene)
outStatsExon$intronSym = NA
outStatsExon$intronSym[queryHits(oo1)] = ensGene$Symbol[subjectHits(oo1)]
oo2 = findOverlaps(exonRegionsMatch, ensGene)
outStatsExon$exonSym = NA
outStatsExon$exonSym[queryHits(oo2)] = ensGene$Symbol[subjectHits(oo2)]
outStatsExon$sameGene = outStatsExon$intronSym == outStatsExon$exonSym 	

####### take signif
outStatsExon$bonf = p.adjust(outStatsExon$pval, "bonf")
outStatsExonSig = outStatsExon[outStatsExon$bonf < 0.05 & 
	outStatsExon$intronSym!="" & outStatsExon$sameGene,]
outStatsExonSig = outStatsExonSig[order(outStatsExonSig$pval),]

## boxplots
intronToPlot = intronMat[match(rownames(outStatsExonSig), rownames(intronMat)),]
exonToPlot = exonMat[match(outStatsExonSig$nearExon, rownames(exonMat)),]


conditionalIntron <- function(i, subset = FALSE) {
    if(subset) {
        if(i == 2) {
            ylim <- c(0, 9)
        } else if (i == 22) {
            ylim <- c(0, 6)
        } else if (i == 25) {
            ylim <- c(0, 10)
        }
    } else {
        ylim <- c(0, 12)
    }
	boxplot(intronToPlot[i,] ~ pd2$Tissue, ylim = ylim,
		ylab="Log2(Adjusted Coverage)",cex.axis=2, cex.lab=2, 
		main="Intronic ER", cex.main=2, outline=FALSE)
    points(x = jitter(tissueToNum[pd2$Tissue]), y = intronToPlot[i,], col = colors[pd2$Tissue], pch = 20, cex = 1.5)
	legend("top", paste0("p=",signif(outStatsExonSig$pval[i], 3)),cex=1.4)
	par(mar = c(5,3,3,2))
	boxplot(exonToPlot[i,] ~ pd2$Tissue, ylim = ylim,
		ylab="",cex.axis=2, cex.lab=2, 
		main="Nearest Exonic ER", cex.main=2, outline=FALSE)
    points(x = jitter(tissueToNum[pd2$Tissue]), y = exonToPlot[i,], col = colors[pd2$Tissue], pch = 20, cex = 1.5)
	mtext(paste(outStatsExonSig$intronSym[i], "-",
		round(outStatsExonSig$nearDist[i]/1000), "kb away"), 
		side=1, outer=TRUE, line = -2, cex=2)
}


tissueToNum <- c('Heart' = 1, 'Liver' = 2, 'Testis' = 3)

pdf("plots/conditional_intronic_ERs_subset.pdf", h=6, w=12)
par(mfrow = c(1,2))
for(i in c(2, 22, 25)) {
	if(i %% 100 == 0) cat(".")
	par(mar = c(5,6,3,0))
	conditionalIntron(i, subset = TRUE)
}
dev.off()

pdf("plots/conditional_intronic_ERs.pdf", h=6, w=12)
par(mfrow = c(1,2))
for(i in seq_len(700)) {
	if(i %% 100 == 0) cat(".")
	par(mar = c(5,6,3,0))
	conditionalIntron(i)
}
dev.off()


save(outStatsExon, outStatsExonSig, file="rdas/conditionalIntronicERs.rda")


## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()

