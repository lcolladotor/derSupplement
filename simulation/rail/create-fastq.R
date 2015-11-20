library('Biostrings')
library('devtools')


read_dir <- file.path('..', 'simulated_reads')
out_dir <- 'simulated_fastq'
dir.create(out_dir)
files <- dir(read_dir, pattern = 'fasta')
qual <- BStringSet(paste(rep('I', 100), collapse = ''))
for(f in files) {
    message(paste(Sys.time(), 'processing file', f))
    reads <- readDNAStringSet(file.path(read_dir, f))
    writeXStringSet(reads, file.path(out_dir, f), format = 'fastq', qualities = rep(qual, length(reads)), compress = TRUE)
}


## Reproducibility info
proc.time()
options(width = 120)
session_info()
