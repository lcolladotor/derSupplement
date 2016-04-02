## Create manifest file for running rail
sampleNames <- paste0(rep(paste0('sample', 1:10, '-G', rep(1:2, each = 5)), 3), '-R', rep(1:3, each = 10))

foo <- function(i) {
    cat(paste0("/dcl01/lieber/ajaffe/derRuns/derSupplement/simulation/simulated_reads/sample_", sprintf('%02d', i), "_1.fasta.gz\t0\t/dcl01/lieber/ajaffe/derRuns/derSupplement/simulation/simulated_reads/sample_", sprintf('%02d', i), "_2.fasta.gz\t0\t", sampleNames[i], "\n"))
}


sink("rail-manifest-R1.txt")
for(i in seq_len(10)) foo(i)
sink()

sink("rail-manifest-R2.txt")
for(i in 10 + seq_len(10)) foo(i)
sink()

sink("rail-manifest-R3.txt")
for(i in 20 + seq_len(10)) foo(i)
sink()
