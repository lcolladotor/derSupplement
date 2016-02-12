library('GenomicRanges')
files <- dir(pattern = 'Rdata')

check <- sapply(files, function(i) {
    print(i)
    load(i)
    if(grepl('DESeq2', i)) 
        res <- is(deseq, 'GRanges')
    else
        res <- is(edger, 'GRanges')
    return(res)
})

all(check)
