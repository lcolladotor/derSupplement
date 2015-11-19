## Create manifest file for running rail
sampleNames <- paste0(rep(paste0('sample', 1:10, 'G', rep(1:2, each = 5)), 3), 'R', rep(1:3, each = 10))

{
sink("rail-manifest.txt")
for(i in seq_len(30)) {
	cat(paste0("/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/simulated_reads/sample_", sprintf('%02d', i), "_1.fasta.gz\t0\t/dcs01/ajaffe/Brain/derRuns/derSupplement/simulation/simulated_reads/sample_", sprintf('%02d', i), "_2.fasta.gz\t0\t", sampleNames[i], "\n"))
}
sink()
}
