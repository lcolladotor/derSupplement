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
load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")
load("/dcs01/ajaffe/Brain/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.25.Rdata")
load("/dcs01/ajaffe/Brain/derRuns/derSupplement/brainspan/derAnalysis/run4-v1.0.10/models.Rdata")
load("rdas/summarized_BrainSpan_DERs.rda") # single base
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

# coverage matrix
fullRegionMat = do.call("rbind",
	lapply(regionMat, function(x) x$coverageMatrix))

## drop regions shorter than 6 bp
keepIndex=which(width(fullRegionGR) >= 6)
fullRegionGR = fullRegionGR[keepIndex]
fullRegionMat = fullRegionMat[keepIndex,]

# total coverage
sum(width(fullRegionGR))/1e6

##### lower cutoff
load("/dcs01/ajaffe/Brain/derRuns/derSupplement/brainspan/regionMatrix/regionMat-cut0.1.Rdata")
regList1 = lapply(regionMat, function(x) x$regions)
fullRegionGR1 = unlist(GRangesList(regList1))
fullRegionMat1 = do.call("rbind",
	lapply(regionMat, function(x) x$coverageMatrix))
keepIndex1=which(width(fullRegionGR1) >= 6)
fullRegionGR1 = fullRegionGR1[keepIndex1]
fullRegionMat1 = fullRegionMat1[keepIndex1,]

## log transform
y = log2(fullRegionMat + 1)
rownames(y) = NULL

## annotate regions based on transcriptome databases
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
ensemblAnno = annotateRegions(fullRegionGR,gs)
countTable = ensemblAnno$countTable
vennDiagram(vennCounts(countTable > 0)); mtext("Ensembl", line=1,cex=2)

mean(countTable$exon == 1 & countTable$intergenic == 0 & 
	countTable$intron== 0)
mean(countTable$exon == 0 & (countTable$intergenic > 0 | 
	countTable$intron > 0))

## with 0.1 cutoff	
ensemblAnno1 = annotateRegions(fullRegionGR1,gs)
countTable1 = ensemblAnno1$countTable
vennDiagram(vennCounts(countTable1 > 0)); mtext("Ensembl", line=1,cex=2)

mean(countTable1$exon == 1 & countTable1$intergenic == 0 & 
	countTable1$intron== 0)
mean(countTable1$exon == 0 & (countTable1$intergenic > 0 | 
	countTable1$intron > 0))

## compare to single base
fit = lmFit(y, models$mod)
fit0 = lmFit(y, models$mod0)
ff = getF(fit,fit0, y)

sum(p.adjust(ff$f_pval,"bonf") < 0.05)
mean(p.adjust(ff$f_pval,"bonf") < 0.05)

sigIndex=which(p.adjust(ff$f_pval,"bonf") < 0.05)

xx = ensemblAnno$annotationList[sigIndex]
theGenes = unlist(unlist(xx)$gene)
theGenes = sapply(theGenes, function(x) x[!is.na(x)])
theGenes = theGenes[!is.na(theGenes)]
length(unique(theGenes))

an = annotateNearest(fullRegionGR, sigSpan)
sum(an$dist==0)

sum(ff$f_pval[an$dist==0] < 0.05/nrow(y))
mean(ff$f_pval[an$dist==0] < 0.05/nrow(y))

an2 = annotateNearest(sigSpan, fullRegionGR)
table(an2$dist==0)

### subset analysis ###

Index1 = which(pdSpan$fetal == "Fetal" & 
	pdSpan$NCX %in% c("HIP", "STR"))
mod1 = model.matrix(~as.character(pdSpan$NCX[Index1]))
colnames(mod1)[2] = "STR"
fit1 = lmFit(y[,Index1], mod1)
eb1 = ebayes(fit1)
sigIndex1 = order(eb1$p[,2])[1:sum(p.adjust(eb1$p[,2], "bonf") < 0.05)]
length(sigIndex1)

theGenes1 = matchGenes(fullRegionGR[sigIndex1])
theGenes1$annotation = ss(theGenes1$annotation, " ")

geneList1 = split(theGenes1$name, sign(eb1$t[sigIndex1,2]))
geneList1 = lapply(geneList1, unique)
lapply(geneList1, head, 50)

