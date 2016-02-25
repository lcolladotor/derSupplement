library('knitr')
load('peaks.Rdata')

colnames(peaks) <- c('Max memory by core (GB)', 'Time (minutes)', 'Memory (GB)', 'Peak cores', 'Pipeline', 'Analysis step')
peaks$Pipeline[nrow(peaks)] <- 'HISAT'
kable(peaks[, c(1, 2, 4, 5, 6)], format = 'latex', row.names = FALSE, digits = 1, caption = '\\textbf{Summary of computing resources required for each analysis step for the different simulation pipelines.} This table shows the maximum memory (GB) per core, the time in minutes to run the analysis with all jobs running sequentially and the maximum number of cores used in any step of the simulation analysis for the different pipelines. Note that the expressed-regions (HISAT), the feature counts and ballgown pipelines rely on HISAT alignments.')