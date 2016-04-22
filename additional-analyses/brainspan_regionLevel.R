## Usage:
# qrsh -l mem_free=130G,h_vmem=150G
# module load R/3.3
# mkdir -p logs
# Rscript brainspan_regionLevel.R > logs/brainspan_regionLevel_log.txt 2>&1

###
library(limma)
library(GenomicRanges)
library(derfinder)
library(bumphunter)
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), function(y) y[slot])
splitit = function(x) split(seq(along=x),x) # splits into list

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
load("/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.25.Rdata")
load("/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/derAnalysis/run5-v1.5.30/models.Rdata")

## Load single base results
if(!file.exists('rdas/summarized_BrainSpan_DERs.rda')) {
    stop("Run characterize_brainspan_DERs.R first")
} else {
    load("rdas/summarized_BrainSpan_DERs.rda")
}

## Fix models
stopifnot(nrow(models$mod) == nrow(models$mod0))
if(unique(sapply(models, nrow)) == 487) {
    models$mod <- models$mod[-bad_samples, ]
    models$mod0 <- matrix(models$mod0[-bad_samples, ], ncol = 1)
}
stopifnot(nrow(models$mod) == 484)
stopifnot(nrow(pdSpan) == 484)

sigSpan$annotation = ss(sigSpan$annotation, " ")

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

# regions
regList = lapply(regionMat, function(x) x$regions)
fullRegionGR = unlist(GRangesList(regList))
print('Number of regions and MB covered: cut 0.25')
length(fullRegionGR)
sum(width(fullRegionGR))/1e6

# coverage matrix
fullRegionMat = do.call("rbind",
	lapply(regionMat, function(x) x$coverageMatrix))
## Drop bad samples
if(ncol(fullRegionMat) == 487) fullRegionMat <- fullRegionMat[, -bad_samples]

stopifnot(ncol(fullRegionMat) == 484)

## Explore single base-level DERs widths (signifcant onlhy)
summary(width(sigSpan))

## drop regions shorter than 6 bp
keepIndex=which(width(fullRegionGR) >= 6)
fullRegionGR = fullRegionGR[keepIndex]
fullRegionMat = fullRegionMat[keepIndex,]

print('Number of regions and MB covered: cut 0.25, >= 6bp')
length(fullRegionGR)
sum(width(fullRegionGR))/1e6

##### lower cutoff
load("/dcl01/lieber/ajaffe/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.1.Rdata")
regList1 = lapply(regionMat, function(x) x$regions)
fullRegionGR1 = unlist(GRangesList(regList1))
fullRegionMat1 = do.call("rbind",
	lapply(regionMat, function(x) x$coverageMatrix))
    
print('Number of regions and MB covered: cut 0.10')
length(fullRegionGR1)
sum(width(fullRegionGR1))/1e6

keepIndex1=which(width(fullRegionGR1) >= 6)
fullRegionGR1 = fullRegionGR1[keepIndex1]
fullRegionMat1 = fullRegionMat1[keepIndex1,]
if(ncol(fullRegionMat1) == 487) fullRegionMat1 <- fullRegionMat1[, -bad_samples]
stopifnot(ncol(fullRegionMat1) == 484)

print('Number of regions and MB covered: cut 0.10, >= 6bp')
length(fullRegionGR1)
sum(width(fullRegionGR1))/1e6

## log transform
y = log2(fullRegionMat + 1)
rownames(y) = NULL

## annotate regions based on transcriptome databases
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
ensemblAnno = annotateRegions(fullRegionGR,gs)
countTable = ensemblAnno$countTable

mean(countTable$exon > 0 & countTable$intergenic == 0 & 
	countTable$intron== 0) * 100
mean(countTable$exon == 0 & (countTable$intergenic > 0 | 
	countTable$intron > 0)) * 100

## with 0.1 cutoff	
ensemblAnno1 = annotateRegions(fullRegionGR1,gs)
countTable1 = ensemblAnno1$countTable

dir.create('plots', showWarnings = FALSE)
pdf("plots/venn_counts_brainspan_regionLevel.pdf",h=5,w=6)
vennDiagram(vennCounts(countTable > 0)); mtext("Ensembl, cutoff = 0.25, ERs >= 6bp", line=1,cex=2)
vennDiagram(vennCounts(countTable1 > 0)); mtext("Ensembl, cutoff = 0.01, ERs >= 6bp", line=1,cex=2)
dev.off()

mean(countTable1$exon > 0 & countTable1$intergenic == 0 & 
	countTable1$intron== 0) * 100
mean(countTable1$exon == 0 & (countTable1$intergenic > 0 | 
	countTable1$intron > 0)) * 100

## compare to single base
fit = lmFit(y, models$mod)
fit0 = lmFit(y, models$mod0)
ff = getF(fit,fit0, y)

sum(p.adjust(ff$f_pval,"bonf") < 0.05)
mean(p.adjust(ff$f_pval,"bonf") < 0.05) * 100

sigIndex=which(p.adjust(ff$f_pval,"bonf") < 0.05)

xx = ensemblAnno$annotationList[sigIndex]
theGenes = unlist(unlist(xx)$gene)
theGenes = sapply(theGenes, function(x) x[!is.na(x)])
theGenes = theGenes[!is.na(theGenes)]

theSymbols <- unlist(unlist(xx)$symbol)
theSymbols <- sapply(theSymbols, function(x) x[!is.na(x)])
theSymbols <- theSymbols[!is.na(theSymbols)]

print("Number of unique ensembl genes, then unique genes with symbols")
length(unique(theGenes))
length(unique(theSymbols))

an = annotateNearest(fullRegionGR, sigSpan)
table(an$dist == 0)
table(an$dist == 0) / length(an$dist) * 100

print("Number of sig ER-level DERs overlapping sig single base-level DERs")
sum(ff$f_pval[an$dist==0] < 0.05/nrow(y))
print("Percent of sig ER-level DERs overlapping sig single base-level DERs")
sum(ff$f_pval[an$dist==0] < 0.05/nrow(y)) / sum(ff$f_pval < 0.05/nrow(y)) * 100
print("Percent of ER-level DERs overlapping sig single base-level DERs that are significant")
mean(ff$f_pval[an$dist==0] < 0.05/nrow(y)) * 100

print("sig SB-level DERs overlapping ER-level DERs")
an2 = annotateNearest(sigSpan, fullRegionGR)
table(an2$dist == 0)
table(an2$dist == 0) / length(an2$dist) * 100

## Load meanCoverage for sig SB-level DERs
if(!file.exists('rdas/summarized_BrainSpan_DERs_meanCov.rda')) {
    stop("Run characterize_brainspan_DERs.R first")
} else {
    load("rdas/summarized_BrainSpan_DERs_meanCov.rda")
}
stopifnot(length(sigSpan) == nrow(meanCoverage))

## Get mean coverage per sig DB-level DER
mean_meanCov <- rowMeans(meanCoverage)

print("Percent of sig SB-level DERs with mean coverage < 0.25")
table(mean_meanCov < 0.25)
table(mean_meanCov < 0.25) / length(mean_meanCov) * 100

print("Percent of sig SB-level DERs not overlapping ER-level DERs with mean coverage < 0.25")
table(mean_meanCov[an2$dist != 0] < 0.25)
mean(mean_meanCov[an2$dist != 0] < 0.25) * 100

print("Width of sig SB-level DERs not overlapping ER-level DERs, then overlapping them, then all of them")
summary(width(sigSpan[an2$dist != 0]))
summary(width(sigSpan[an2$dist == 0]))
summary(width(sigSpan))

### subset analysis ###

Index1 = which(pdSpan$fetal == "Fetal" & 
	pdSpan$NCX %in% c("HIP", "STR"))
mod1 = model.matrix(~as.character(pdSpan$NCX[Index1]))
colnames(mod1)[2] = "STR"
fit1 = lmFit(y[,Index1], mod1)
eb1 = ebayes(fit1)
sigIndex1 = order(eb1$p[,2])[1:sum(p.adjust(eb1$p[,2], "bonf") < 0.05)]
length(sigIndex1)

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
theGenes1 = matchGenes(fullRegionGR[sigIndex1], genes)
theGenes1$annotation = ss(theGenes1$annotation, " ")

geneList1 = split(theGenes1$name, sign(eb1$t[sigIndex1,2]))
geneList1 = lapply(geneList1, unique)
print('Subset analysis results, top 50 in each direction, then number of unique genes')
lapply(geneList1, head, 50)
length(unique(unlist(geneList1)))

## Reproducibility info
library('devtools')
options(width = 120)
session_info()
Sys.time()
proc.time()
