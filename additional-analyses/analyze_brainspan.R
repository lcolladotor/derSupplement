## Usage:
# qrsh -l mem_free=80G,h_vmem=150G
# module load R/3.3
# mkdir -p logs
# Rscript analyze_brainspan.R > logs/analyze_brainspan_log.txt 2>&1

###
library(limma)
library(GenomicRanges)
library(derfinder)
library(bumphunter)
library(RColorBrewer)

ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), function(y) y[slot])
splitit = function(x) split(seq(along=x),x) # splits into list
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100

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

# load data
load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")
load("/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.25.Rdata")
load("/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/derAnalysis/run5-v1.5.30/models.Rdata")

## Remove bad samples
bad_samples <- which(rownames(pdSpan) %in% c('216', '218', '219'))
pdSpan[bad_samples, ]
if(nrow(pdSpan) == 487) pdSpan <- pdSpan[-bad_samples, ]
stopifnot(nrow(pdSpan) == 484)
stopifnot(nrow(models$mod) == nrow(models$mod0))
if(unique(sapply(models, nrow)) == 487) {
    models$mod <- models$mod[-bad_samples, ]
    models$mod0 <- matrix(models$mod0[-bad_samples, ], ncol = 1)
}
stopifnot(nrow(models$mod) == 484)

## add pheno info
pdSpan$fetal = ifelse(pdSpan$Age < 0, "Fetal", "Postnatal")
pdSpan$fetal = factor(pdSpan$fetal,levels=c("Postnatal","Fetal"))
pdSpan$struct = factor(pdSpan$structure_acronym, levels = c("DFC","VFC","MFC",
	"OFC","M1C","S1C", "IPC", "A1C", "STC", "ITC", "V1C", "HIP",
	"AMY", "STR", "MD", "CBC"))
ncx = as.character(pdSpan$struct)
ncx[ncx %in% c("DFC","VFC","MFC",
	"OFC","M1C","S1C", "IPC", "A1C", "STC", "ITC", "V1C")] = "NCX"
pdSpan$NCX = factor(ncx, levels = c("NCX",  "HIP",
	"AMY", "STR", "MD", "CBC"))
pdSpan$Group = with(pdSpan, paste0(NCX, ":", fetal))
pdSpan$Group = factor(pdSpan$Group, levels = 
	paste0(rep(levels(pdSpan$NCX), each=2), ":", 
		rep(c("Fetal","Postnatal"), times=6)))

### extract coverage data
regions = unlist(GRangesList(lapply(regionMat, '[[', 'regions')))
names(regions) = NULL
regionMat = do.call("rbind", lapply(regionMat, '[[', 'coverageMatrix'))
## Samples are in order
if(ncol(regionMat) == 487) {
    stopifnot(identical(colnames(regionMat)[-bad_samples], pdSpan$lab))
} else {
    stopifnot(identical(colnames(regionMat), pdSpan$lab))
}
rownames(regionMat) = names(regions) = paste0("er", seq_len(nrow(regionMat)))

if(ncol(regionMat) == 487) regionMat <- regionMat[, -bad_samples]

stopifnot(ncol(regionMat) == 484)

### filter out short DERs
keepIndex = width(regions) > 8
regionMat = regionMat[keepIndex, ]
regions = regions[keepIndex]

# total coverage
sum(width(regions))/1e6

## log transform
y = log2(regionMat + 1)

####################
## annotate regions based on transcriptome databases
####################

## make genomic state
if(!file.exists('rdas/GenomicState_knownGene_derfinderPaper.rda')) {
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    GenomicState_knownGene = makeGenomicState(
    	TxDb.Hsapiens.UCSC.hg19.knownGene,
    	chrs = paste0("chr", c(1:22, "X", "Y","M")))$fullGenome
    dir.create('rdas', showWarnings = FALSE)
    save(GenomicState_knownGene, compress=TRUE,
    	file="rdas/GenomicState_knownGene_derfinderPaper.rda")
} else {
    load('rdas/GenomicState_knownGene_derfinderPaper.rda')
}


#### annotate
ucscAnno = annotateRegions(regions,GenomicState_knownGene)
countTable = ucscAnno$countTable

dir.create('plots', showWarnings = FALSE)
pdf("plots/venn_counts_analyze_brainspan.pdf",h=5,w=6)
vennDiagram(vennCounts(countTable > 0))
mtext("UCSC", line=1,cex=2)
dev.off()

## annotation breakdown ####
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

#########################
### PCA ###

## by annotation class
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

### plots
pdf("plots/brainspan_regionMatrix_PCA_byAnno.pdf")
palette(brewer.pal(3,"Set1"))
name = c("Exonic", "Intronic", "Intergenic","Exon+Intron", "All")
par(mar=c(5,6,2,2))
for(i in 1:ncol(pc1Mat)) {
	plot(x=pc1Mat[,i], y=pc2Mat[,i], 
		bg = as.numeric(pdSpan$fetal),
		pch = rep(c(21,22,24),times=c(11,4,1))[as.numeric(pdSpan$struct )],
		xlab = paste0("PC1: ",pcVarMat[1,i],"% of Var Expl"),
		ylab = paste0("PC2: ",pcVarMat[2,i],"% of Var Expl"),
		cex.axis=2,cex.lab=2, cex.main=1.8,
		main = paste0("PCA of Expressed Regions (", name[i],")"))
	legend("topleft", c("NCX","Non-NCX","CBC"), 
		pch=c(19,15,17), cex=1.5)
	legend("bottomleft", c("Fetal","Postnatal"), 
		col = 1:2, lwd=5,cex=1.5)
}
dev.off()

pdf("plots/brainspan_regionMatrix_PCsByRegion_byAnno.pdf",w=11)
palette(brewer.pal(3,"Set1"))
par(mar=c(11,6,2,2))
for(i in 1:ncol(pc1Mat)) {
	## PC1
	boxplot(pc1Mat[,i] ~ pdSpan$Group, las=3,
		ylab = paste0("PC1: ",pcVarMat[1,i],"% of Var Expl"),
		cex.axis=1.7,cex.lab=2,cex.main=1.8, xlab="",outline=FALSE,
		main = paste0("PCA of Expressed Regions (", name[i],")"))
	points(pc1Mat[,i] ~ jitter(as.numeric(pdSpan$Group), amount=0.2),
		bg = as.numeric(pdSpan$fetal),cex=1.3,
		pch = rep(c(21,22,24),times=c(11,4,1))[as.numeric(pdSpan$struct)])
	# PC2 
	boxplot(pc2Mat[,i] ~ pdSpan$Group, las=3,
		ylab = paste0("PC2: ",pcVarMat[2,i],"% of Var Expl"),
		cex.axis=1.7,cex.lab=2,xlab="",outline=FALSE)
	points(pc2Mat[,i] ~ jitter(as.numeric(pdSpan$Group), amount=0.2),
		bg = as.numeric(pdSpan$fetal),cex=1.3,
		pch = rep(c(21,22,24),times=c(11,4,1))[as.numeric(pdSpan$struct)])
}
dev.off()

#################
## DE analysis ##

fit = lmFit(y, models$mod)
fit0 = lmFit(y, models$mod0)
ff = getF(fit,fit0, y)

sum(p.adjust(ff$f_pval,"bonf") < 0.05)
mean(p.adjust(ff$f_pval,"bonf") < 0.05)

sigIndex=which(p.adjust(ff$f_pval,"bonf") < 0.05)

xx = ucscAnno$annotationList[sigIndex]
theGenes = unlist(unlist(xx)$gene)
theGenes = sapply(theGenes, function(x) x[!is.na(x)])
theGenes = theGenes[!is.na(theGenes)]
length(unique(theGenes))


if(!file.exists('rdas/summarized_BrainSpan_DERs.rda')) {
    stop("Run characterize_brainspan_DERs.R first")
} else {
    pdSpan_more <- pdSpan
    load("rdas/summarized_BrainSpan_DERs.rda")
}


an = annotateNearest(regions, sigSpan)
sum(an$dist==0)

sum(ff$f_pval[an$dist==0] < 0.05/nrow(y))
mean(ff$f_pval[an$dist==0] < 0.05/nrow(y))

an2 = annotateNearest(sigSpan, regions)
table(an2$dist==0)

### subset analysis ###

Index1 = which(pdSpan_more$fetal == "Fetal" & 
	pdSpan_more$NCX %in% c("HIP", "STR"))
mod1 = model.matrix(~as.character(pdSpan_more$NCX[Index1]))
colnames(mod1)[2] = "STR"
fit1 = lmFit(y[,Index1], mod1)
eb1 = ebayes(fit1)
sigIndex1 = order(eb1$p[,2])[1:sum(p.adjust(eb1$p[,2], "bonf") < 0.05)]
length(sigIndex1)

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
         
theGenes1 = matchGenes(regions[sigIndex1], genes)
theGenes1$annotation = ss(theGenes1$annotation, " ")

geneList1 = split(theGenes1$name, sign(eb1$t[sigIndex1, 2]))
geneList1 = lapply(geneList1, unique)
lapply(geneList1, head, 50)



##### lower cutoff
load("/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.1.Rdata")
regList1 = lapply(regionMat, function(x) x$regions)
regions1 = unlist(GRangesList(regList1))
fullRegionMat1 = do.call("rbind",
	lapply(regionMat, function(x) x$coverageMatrix))
if(ncol(fullRegionMat1) == 487) fullRegionMat1 <- fullRegionMat1[, -bad_samples]
stopifnot(ncol(fullRegionMat1) == 484)
keepIndex1=which(width(regions1) >= 6)
regions1 = regions1[keepIndex1]
fullRegionMat1 = fullRegionMat1[keepIndex1,]

## Reproducibility info
library('devtools')
options(width = 120)
session_info()
Sys.time()
proc.time()
