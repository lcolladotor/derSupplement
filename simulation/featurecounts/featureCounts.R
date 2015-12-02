library('Rsubread')
library('devtools')

## Identify BAM files
files <- dir('/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/hisat', pattern = '[0-9]\\.bam$', full.names = TRUE)
names(files) <- gsub('\\.bam', '', dir('/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/hisat', pattern = '[0-9]\\.bam$'))


message(paste(Sys.time(), 'running featureCounts'))
featCounts <- featureCounts(files = files, annot.ext = '/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/gtf/chr17.gtf', isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, isPairedEnd = TRUE, nthreads = 4)

## Fix sample names
colnames(featCounts$counts) <- gsub('.bam', '', gsub('X.dcs01.ajaffe.Brain.derRuns.derSupplement.simulation.hisat.', '', colnames(featCounts$counts)))

colnames(featCounts$stat) <- gsub('.bam', '', gsub('X.dcs01.ajaffe.Brain.derRuns.derSupplement.simulation.hisat.', '', colnames(featCounts$stat)))

## Print stats
stat <- t(featCounts$stat[, -1])
colnames(stat) <- featCounts$stat[, 1]
rownames(stat) <- colnames(featCounts$stat)[-1]

## Summary first, then all the info
summary(stat)
options(width = 120)
stat

message(paste(Sys.time(), 'saving featureCounts output'))
save(featCounts, stat, file = 'featureCounts.Rdata')

## Reproducibility info
Sys.time()
proc.time()
session_info()
