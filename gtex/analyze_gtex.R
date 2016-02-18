## Original script: /home/epi/ajaffe/Lieber/Projects/derfinderPaper/analyze_gtex.R
##
library(derfinder)
library(GenomicRanges)
library(rafalib)
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

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
load('/dcs01/ajaffe/Brain/derRuns/railDER/gtex36/regionMat-cut5.Rdata')
load('/dcs01/ajaffe/Brain/derRuns/derSoftware/gtex/mappedInfo.Rdata')
xx = load("/home/epi/ajaffe/Lieber/Projects/derfinderPaper/gtex_pheno.rda")
mappedInfo$Tissue = pd1$SMTS[match(mappedInfo$sample, pd1$sra_accession)]
mappedInfo$SubjectID = pd1$SUBJID[match(mappedInfo$sample, pd1$sra_accession)]

### extract coverage data
regions = unlist(GRangesList(lapply(regionMat, '[[', 'regions')))
names(regions) = NULL
regionMat = do.call("rbind", lapply(regionMat, '[[', 'coverageMatrix'))
regionMat = regionMat[,mappedInfo$sample] # put in order
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
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
ensemblAnno = annotateRegions(regions,gs)
countTable = ensemblAnno$countTable

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
sapply(annoClassList, length)
100*sapply(annoClassList, length)/nrow(countTable)

# width by annotation
t(sapply(annoClassList, function(ii) quantile(width(regions[ii]))))

### PCA ###
annoClassList$All = 1:nrow(regionMat) # add all
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
rafalib::mypar(2,2,cex.axis=1)
for(i in ind) boxplot(pc1Mat[,i] ~ mappedInfo$Tissue, 
	main = colnames(pc1Mat)[i],
	ylab=paste0("PC1: ", pcVarMat[1,i], "% of Var Explain"))

for(i in ind) boxplot(pc2Mat[,i] ~ mappedInfo$Tissue, 
	main = colnames(pc2Mat)[i],
	ylab=paste0("PC2: ", pcVarMat[2,i], "% of Var Explain"))
	
#################
## DE analysis ##
library(limma)
mod = model.matrix(~mappedInfo$Tissue)
mod0 = model.matrix(~1, data=mappedInfo)
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

# conditional regression
outStatsExon = matrix(NA, ncol = 2, nrow = nrow(intronMat))
for(i in 1:nrow(outStatsExon)) {
	if(i %% 1000 == 0) cat(".")
	f = lm(intronMat[i,] ~ 
			mappedInfo$Tissue + exonMatMatch[i,])
	f0 = lm(intronMat[i,] ~ exonMatMatch[i,])
	outStatsExon[i,] = as.numeric(anova(f,f0)[2,5:6])
}
colnames(outStatsExon) = c("Fstat", "pval")
rownames(outStatsExon) = rownames(intronMat)
outStatsExon=as.data.frame(outStatsExon)

outStatsExon$nearExon = rownames(exonMatMatch)
outStatsExon$nearDist = mcols(ooExon)$distance

## get gene symbol
library(GenomicFeatures)
TranscriptDb=loadDb("/home/epi/ajaffe/Lieber/Projects/RNAseq/Ribozero_Compare/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12/inst/extdata/TxDb.Hsapiens.BioMart.ensembl.GRCh37.p12.sqlite")
seqlevels(TranscriptDb,force=TRUE) = c(1:22,"X","Y","MT")
seqlevels(TranscriptDb) = paste0("chr", c(1:22,"X","Y","M"))
ensGene = genes(TranscriptDb)

library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", # VERSION 75, hg19
	dataset="hsapiens_gene_ensembl",
	host="feb2014.archive.ensembl.org")
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

pdf("plots/conditional_intronic_ERs.pdf",h=6,w=12)
par(mfrow = c(1,2))
for(i in 1:1000) {
	if(i %% 100 == 0) cat(".")
	par(mar = c(5,6,3,0))
	boxplot(intronToPlot[i,] ~ mappedInfo$Tissue, ylim = c(0,12),
		ylab="Log2(Adjusted Coverage)",cex.axis=2, cex.lab=2, 
		main="Intronic ER", cex.main=2)
	legend("top", paste0("p=",signif(outStatsExonSig$pval[i], 3)),cex=1.4)
	par(mar = c(5,3,3,2))
	boxplot(exonToPlot[i,] ~ mappedInfo$Tissue, ylim = c(0,12),
		ylab="",cex.axis=2, cex.lab=2, 
		main="Nearest Exonic ER", cex.main=2)
	mtext(paste(outStatsExonSig$intronSym[i], "-",
		round(outStatsExonSig$nearDist[i]/1000), "kb away"), 
		side=1, outer=TRUE, line = -2, cex=2)
}
dev.off()

save(outStatsExon, outStatsExonSig, file="rdas/conditionalIntronicERs.rda")
