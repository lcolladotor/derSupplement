library('Biostrings')
library('devtools')


read_dir <- file.path('..', 'simulated_reads')
out_dir <- 'simulated_fastq'
dir.create(out_dir)
files <- dir(read_dir, pattern = 'fasta')
qual <- BStringSet(paste(rep('I', 100), collapse = ''))
for(f in files) {
    message(paste(Sys.time(), 'processing file', f))
    
    ## Load data
    reads <- readDNAStringSet(file.path(read_dir, f))
    
    ## Sort by read length
    read_l <- elementLengths(reads)
    reads <- reads[order(read_l, decreasing = TRUE)]
    
    ## Create fake quality strings (using I)
    tab_l <- table(read_l)
    tab_l <- tab_l[order(as.integer(names(tab_l)), decreasing = TRUE)]
    quals <- rep(BStringSet(sapply(as.integer(names(tab_l)), function(x) { paste(rep('I', x), collapse = '') })), tab_l)
    
    ## Check before writing
    stopifnot(identical(elementLengths(quals), elementLengths(reads)))
    
    ## Write fastq file
    writeXStringSet(reads, file.path(out_dir, gsub('fasta', 'fastq', f)), format = 'fastq', qualities = quals, compress = TRUE)
}

## Reproducibility info
proc.time()
options(width = 120)
session_info()
