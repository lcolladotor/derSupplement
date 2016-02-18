###
# qsub -V -pe local 24 -l jabba,mf=80G,h_vmem=16G,h_stack=256M -cwd -b y R CMD BATCH --no-save feature_counts.R
source("/home/epi/ajaffe/Lieber/lieber_functions_aj.R") 

library(Rsubread)

##  stem
xx=load("/home/epi/ajaffe/Lieber/Projects/RNAseq/UCSD_samples/UCSD_stemcell_pheno.rda")
pdStem = pd
pdStem$bamFile = paste0("/dcs01/ajaffe/UCSC_Epigenome/RNAseq/TopHat/",	
	pdStem$sample, "_out/accepted_hits.bam")
	
## hippo
load("/home/epi/ajaffe/Lieber/Projects/RNAseq/HippoPublic/sra_phenotype_file.rda")
pdHippo = sra
pdHippo = pdHippo[-grep("ED", pdHippo$SampleID),]

## snyder
bam = list.files("/dcs01/ajaffe/Snyder/RNAseq/TopHat",
		pattern="accepted_hits.bam$", recur=TRUE,full.names=TRUE)
names(bam) = list.files("/dcs01/ajaffe/Snyder/RNAseq/TopHat",
		pattern="accepted_hits.bam$", recur=TRUE)
	
pheno  = data.frame(sampleName = c(pdStem$run, pdHippo$SampleID, ss(names(bam), "_")), 
	Study = c(rep("Stem", nrow(pdStem)), rep("Hippo", nrow(pdHippo)),
		rep("Snyder", length(bam))),
	bamFile = c(pdStem$bamFile, pdHippo$bamFile, bam),
	stringsAsFactors=FALSE)
libSize = getTotalMapped(pheno$bamFile,mc.cores=12)
pheno$totalMapped = libSize$totalMapped
pheno$mitoMapped = libSize$mitoMapped

## count genes
geneCounts = featureCounts(pheno$bamFile, annot.inbuilt="hg19",
	useMetaFeatures = TRUE,	nthreads=24)
exonCounts = featureCounts(pheno$bamFile, annot.inbuilt="hg19", 
	useMetaFeatures = FALSE,nthreads=24)
save(geneCounts,exonCounts, pheno, file="featureCounts_output.rda")

#####
load("featureCounts_output.rda")

pheno$geneAlign = as.numeric(geneCounts$stat[1,-1] / colSums(geneCounts$stat[,-1]))
pheno$geneAmbig = as.numeric(geneCounts$stat[2,-1] / colSums(geneCounts$stat[,-1]))
pheno$exonAlign = as.numeric(exonCounts$stat[1,-1] / colSums(exonCounts$stat[,-1]))
pheno$exonAmbig = as.numeric(exonCounts$stat[2,-1] / colSums(exonCounts$stat[,-1]))

signif(tapply(pheno$geneAlign, pheno$Study, mean),3)
signif(tapply(pheno$geneAmbig, pheno$Study, mean),3)
signif(tapply(pheno$exonAlign, pheno$Study, mean),3)
signif(tapply(pheno$exonAmbig, pheno$Study, mean),3)

### mito map
signif(tapply(pheno$mitoMapped/(pheno$mitoMapped+pheno$totalMapped), pheno$Study, mean),3)
