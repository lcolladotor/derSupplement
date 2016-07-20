## Print to Latex

library("xtable")
load('results-exons_one.Rdata')

empirical_exons_one_sum$AnnotationComplete <- ifelse(empirical_exons_one_sum$AnnotationComplete, 'Yes', 'No')
#empirical_exons_one_sum$SummaryMethod[empirical_exons_one_sum$SummaryMethod == 'featureCounts'] <- 'Rsubread'
colnames(empirical_exons_one_sum)[4] <- 'Annotation complete'
colnames(empirical_exons_one_sum)[6] <- 'Summary method'
colnames(empirical_exons_one_sum)[7] <- 'Stat. method'

empirical_exons_one_sum[which(!empirical_exons_one_sum[,'Stat. method'] %in% c('edgeR', 'DESeq2'))[c(7:8, 1:2, 3:6)],]
print.xtable(xtable(empirical_exons_one_sum[which(!empirical_exons_one_sum[,'Stat. method'] %in% c('edgeR', 'DESeq2'))[c(7:8, 1:2, 3:6)],], caption = '\\textbf{Minimum and maximum empirical power, false positive rate (FPR) and false discovery rate (FDR) observed from the three simulation replicates for each analysis pipeline.} \\texttt{ballgown} analyses were done at either the exon or transcript levels. Pipelines that rely on annotation were run with the full annotation or with 20\\% of the transcripts missing (8.28\\% exons missing). Count matrices were analyzed with \\texttt{limma}, \\texttt{DESeq2} and \\texttt{edgeR}-robust (Supplementary Table \\ref{tab:sim2edgeR}). FDR of 5\\% was targeted.', label = 'tab:sim2', digits = 2), include.rownames = FALSE, table.placement = '!ht')

empirical_exons_one_sum[which(empirical_exons_one_sum[,'Stat. method'] %in% c('edgeR', 'DESeq2'))[c(5:8, 1:4)],]
print.xtable(xtable(empirical_exons_one_sum[which(empirical_exons_one_sum[,'Stat. method'] %in% c('edgeR', 'DESeq2'))[c(5:8, 1:4)],], caption = '\\textbf{Simulation results for pipelines that used \\texttt{DESeq2} or \\texttt{edgeR}-robust for the statistical tests.}
Minimum and maximum empirical power, false positive rate (FPR) and false discovery rate (FDR) observed from the three simulation replicates for each analysis pipeline that resulted in a count matrix analyzed with \\texttt{DESeq2} or \\texttt{edgeR}-robust. \\texttt{ballgown} analyses were done at either the exon or transcript levels. Pipelines that rely on annotation were run with the full annotation or with 20\\% of the transcripts missing (8.28\\% exons missing).', label = 'tab:sim2edgeR', digits = 2), include.rownames = FALSE, table.placement = '!ht')
