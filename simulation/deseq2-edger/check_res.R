library('GenomicRanges')
files <- dir(pattern = 'Rdata')

check <- sapply(files, function(i) {
    print(i)
    load(i)
    if(grepl('DESeq2', i))
        res <- is(deseq, 'GRanges')
    else if(grepl('edgeR', i))
        res <- is(edger, 'GRanges')
    else
        res <- is(limma, 'GRanges')
    return(res)
})

all(check)
