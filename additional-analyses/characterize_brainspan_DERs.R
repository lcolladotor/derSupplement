###
## Title: Calculate final DERs and explore
## By Andrew Jaffe
## needs: R

## Usage:
## note that the initial run uses quite a bit of memory, but the second time
## not as much is needed
# qrsh -l mem_free=200G,h_vmem=300G
# module load R/devel
# mkdir -p logs
# Rscript characterize_brainspan_DERs.R > logs/characterize_brainspan_DERs_log.txt 2>&1

source("/home/epi/ajaffe/Lieber/lieber_functions_aj.R") 

library(derfinder)
library(GenomicRanges)
 load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")

#path = "/dcs01/ajaffe/Brain/derRuns/derSupplement/brainspan/derAnalysis/run4-v1.0.10/"
path = "/dcs01/ajaffe/Brain/derRuns/derSoftware/brainspan/derAnalysis/run4-v1.0.10/"

# load in DERs from the prep file
load(paste0(path,"groupInfo.Rdata"))
load(paste0(path,"fullAnnotatedRegions.Rdata"))

## Remove bad samples
bad_samples <- which(rownames(pdSpan) %in% c('216', '218', '219'))
pdSpan[bad_samples, ]
if(nrow(pdSpan) == 487) pdSpan <- pdSpan[-bad_samples, ]
stopifnot(nrow(pdSpan) == 484)
if(length(groupInfo) == 487) groupInfo <- groupInfo[-bad_samples]
stopifnot(length(groupInfo) == 484)

dir.create('rdas', showWarnings = FALSE)
if(!file.exists('rdas/summarized_BrainSpan_DERs.rda')) {
    load(paste0(path, "fullRegions.Rdata"))
    # # load coverage
    

    #####################
    ### significant
    sigSpan = fullRegions[fullRegions$significantFWER == "TRUE"]

    
    save(pdSpan, sigSpan, bad_samples, file = "rdas/summarized_BrainSpan_DERs.rda")
} else {
    load("rdas/summarized_BrainSpan_DERs.rda")
}

if(!file.exists('rdas/summarized_BrainSpan_DERs_meanCov.rda')) {
    #load("/dcs01/ajaffe/Brain/derRuns/derSupplement/brainspan/CoverageInfo/fullCov.Rdata")
    load("/dcs01/ajaffe/Brain/derRuns/derSoftware/brainspan/CoverageInfo/bioarxiv_version/fullCov.Rdata")
     
     coverList = getRegionCoverage(fullCov, sigSpan, mc.cores=1)
     meanCoverage = t(sapply(coverList, colMeans))
     if(ncol(meanCoverage) == 487) meanCoverage <- meanCoverage[, -bad_samples]
     stopifnot(ncol(meanCoverage) == 484)
     colnames(meanCoverage) = pdSpan$lab
     
     save(meanCoverage, file = "rdas/summarized_BrainSpan_DERs_meanCov.rda")
} else {
    load('rdas/summarized_BrainSpan_DERs_meanCov.rda')
}



pdSpan$groupInfo= groupInfo
sigSpan$annotation = ss(sigSpan$annotation, " ")

## how many DERs?
length(sigSpan)

# how much genome covered?
sum(width(sigSpan))/1e6
#### pca ####

pdSpan$fetal = ifelse(pdSpan$Age < 0, "Fetal", "Postnatal")
pdSpan$fetal = factor(pdSpan$fetal,levels=c("Postnatal","Fetal"))

## highest by group
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
		
groupMeans = sapply(splitit(pdSpan$Group), function(i) rowMeans(meanCoverage[,i]))
highGroup = colnames(groupMeans)[apply(groupMeans, 1, which.max)]
table(highGroup)
tab=table(ss(highGroup,":"), ss(highGroup,":",2))
tab = tab[levels(pdSpan$NCX),]

## GO?
gIndexes=splitit(factor(highGroup,levels=levels(pdSpan$Group)))
nullgenes =  read.delim("/home/epi/ajaffe/Lieber/Projects/450k/grant/ref_gene_hg19.txt", 
	header=TRUE,as.is=TRUE)
goByGroup = mclapply(gIndexes, function(ii) {
	cat(".")
	sig2 = sigSpan[ii]
	g = sig2$annotation[!(sig2$description %in% c("upstream","downstream") & 
		sig2$distance > 500)]
	g = g[!is.na(g)]
	go = dogo(g, nullgenes[,2])
	go[,-8]
},mc.cores=12)

save(goByGroup,file="rdas/go_output.rda")


######## PCA
if(!file.exists('rdas/brainspan_der_pca.rda')) {
    pca = prcomp(t(log2(meanCoverage+1)))
    pca$rot = pca$rot[,1:10]
    save(pca, file="rdas/brainspan_der_pca.rda")
} else {
    load("rdas/brainspan_der_pca.rda")  
}


pcaVars = getPcaVars(pca)

levels(groupInfo) = c("NCX.F", "NCX.P",
	"NonNCX.F", "NonNCX.P","CBC.F","CBC.P")

dir.create('plots', showWarnings = FALSE)
pdf("plots/brainspan_pcs_ders.pdf")
palette(brewer.pal(6,"Dark2"))
par(mar=c(5,6,2,2))
for(i in 1:10) {
	plot(pca$x[,i], pca$x[,i+1], 
		bg = as.numeric(groupInfo),
		pch = 21,	xlab = paste0("PC",i,": ",pcaVars[i],"% of Variance Explained"),
		ylab = paste0("PC",i+1,": ",pcaVars[i+1],"% of Variance Explained"),
		cex.axis=2,cex.lab=2)
	if(i == 1) {
		legend("topleft", levels(groupInfo), 
			col = seq(along=levels(groupInfo)), 
			lwd=5,cex=1.2,nc=2)
	}
}
dev.off()
	

pdf("plots/brainspan_pcs_ders_boxplot.pdf",w=11)
palette(brewer.pal(6,"Dark2"))
par(mar=c(11,6,2,2))
for(i in 1:10) {
	boxplot(pca$x[,i] ~ pdSpan$Group, las=3,
		ylab = paste0("PC",i,": ",pcaVars[i],"% of Variance"),
		cex.axis=1.7,cex.lab=2,xlab="",outline=FALSE)
	points(pca$x[,i] ~ jitter(as.numeric(pdSpan$Group), amount=0.2),
		bg = as.numeric(groupInfo),	pch = 21)

}
dev.off()

if(interactive()) plot(pca$x[,2] ~ as.numeric(pdSpan$RIN))

## annotate regions based on transcriptome databases
countTable = fullAnnotatedRegions$countTable[seq(along=sigSpan),]
colnames(countTable)[2] = "intergenic"

### numbers for the paper
sum(countTable[,"intron"] > 0)
mean(countTable[,"intron"] > 0)
sum(countTable[,"intergenic"] > 0 &	
	countTable[,"exon"] == 0 & 	countTable[,"intron"] == 0)
mean(countTable[,"intergenic"] > 0 & 
	countTable[,"exon"] == 0 & countTable[,"intron"] == 0)
sum(countTable[,"exon"] > 0)
mean(countTable[,"exon"] > 0)

## compare
cols = rep("Intergenic", nrow(countTable))
cols[countTable[,"intron"] > 0] = "Intronic"
cols[countTable[,"exon"] > 0] = "Exonic"

tab = table(highGroup,cols)[levels(pdSpan$Group),]
tab = cbind(tab, rowSums(tab))
colnames(tab)[4] = "Total"
write.csv(tab, file="brainspan_der_expression.csv")

type=c("Intergenic", "Intronic","Exonic")

#######
## venn diagram of counts
library(limma)
pdf("plots/venn_counts.pdf",h=5,w=6)
vennDiagram(vennCounts(ensemblCount > 0)); mtext("Ensembl", line=1,cex=2)
vennDiagram(vennCounts(ucscCount > 0)); mtext("UCSC", line=1,cex=2)
vennDiagram(vennCounts(gencodeCount > 0)); mtext("Gencode", line=1,cex=2)
dev.off()

##########

## load libd data
xx=load("/dcs01/ajaffe/Brain/derRuns/libd_n36/derCoverageInfo/fullCov.Rdata")
names(fullCov) = paste0("chr", names(fullCov))

coverListLibd = getRegionCoverage(fullCov,sigSpan,mc.cores=4)
meanCoverageLibd = t(sapply(coverListLibd, colMeans))
save(meanCoverageLibd, file = "rdas/mean_LIBD_cover_BrainSpan_DERs.rda")

## Reproducibility info
library('devtools')
options(width = 120)
session_info()
Sys.time()
proc.time()
