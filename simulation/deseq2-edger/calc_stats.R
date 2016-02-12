library('GenomicRanges')
library('devtools')
library('getopt')

## Specify parameters
spec <- matrix(c(
    'replicate', 'r', 1, 'integer', 'Replicate number. Either 1, 2 or 3.',
    'complete', 'c', 2, 'character', "'yes' or 'no' for whether the GTF was complete or not",
    'pipeline', 'p', 1, 'character', "'featureCounts', 'regionMatrix' or 'railMatrix'",
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## For testing
if(FALSE) opt <- list(replicate = 1, complete = 'yes', pipeline = 'featureCounts')
    

## Load the necessary data
if(opt$pipeline == 'featureCounts') {
    load(file.path('..', ifelse(opt$complete == 'yes', 'featurecounts', 'featurecounts-inc'), paste0('featureCounts-R', opt$replicate, '.Rdata')))
    
    counts <- featCounts$counts
    
    ## Create the appropriate GRangesList if necessary
    if(any(grepl(';', featCounts$annotation))) {
        regions <- with(featCounts$annotation, makeGRangesListFromFeatureFragments(seqnames = sapply(strsplit(Chr, ';'), '[[', 1), fragmentStart = Start, fragmentEnds = End, strand = substr(Strand, 1, 1), sep = ';'))
    } else {
        regions <- with(featCounts$annotation, GRanges(seqnames = Chr, ranges = IRanges(start = Start, end = End), strand = Strand))
    }
    names(regions) <- featCounts$annotation$GeneID
    
    file <- paste0('stats-featureCounts-R', opt$replicate, '-', ifelse(opt$complete == 'yes', 'comp', 'inc'))
    
} else if (opt$pipeline == 'regionMatrix') {
    load(file.path('..', 'regionMatrix', paste0('regionMat-R', opt$replicate, '.Rdata')))
    
    counts <- regionMat$chr17$coverageMatrix
    regions <- regionMat$chr17$regions
    file <- paste0('stats-regionMatrix-R', opt$replicate)
    
} else if (opt$pipeline == 'railatrix') {
    load(file.path('..', 'railMatrix', paste0('regionMat-R', opt$replicate, '.Rdata')))
    
    counts <- regionMat$chr17$coverageMatrix
    regions <- regionMat$chr17$regions
    file <- paste0('stats-railMatrix-R', opt$replicate)
    
} else {
    stop("Invalid 'pipeline' argument")
}





## DESeq2 analysis
run_deseq <- function(counts, regions, file = NULL) {
    counts <- round(counts, 0)
    groupInfo <- as.factor(as.integer(gsub('.*G|R.*', '', colnames(counts))))
    nonzero <- sapply(rowSums(counts), function(x) {x > 0})
    
    ## Round matrix and specify design
    library('DESeq2')
    dse <- DESeqDataSetFromMatrix(counts[nonzero, ], data.frame(group = groupInfo), ~ group)

    ## Perform DE analysis
    system.time( dse <- DESeq(dse, test = 'LRT', reduced = ~ 1) )

    ## Extract results
    deseq <- regions[nonzero]
    mcols(deseq) <- cbind(mcols(deseq), results(dse))

    ## Which are significant?
    mcols(deseq)$sig <- mcols(deseq)$padj < 0.05
    mcols(deseq)$sig[is.na(mcols(deseq)$sig)] <- FALSE

    ## Save results
    if(!is.null(file))
    save(deseq, file = paste0(file, '-DESeq2.Rdata'))
    
    ## End
    return(deseq)
}

## edgeR analysis
run_edger <- function(counts, regions, file = NULL) {
    counts <- round(counts, 0)
    groupInfo <- as.factor(as.integer(gsub('.*G|R.*', '', colnames(counts))))
    nonzero <- sapply(rowSums(counts), function(x) {x > 0})
    
    ## Determine design matrix
    design <- model.matrix(~ groupInfo)

    ## Perform DE analysis
    library('edgeR')
    d <- DGEList(counts = counts[nonzero, ], group = groupInfo)
    d <- calcNormFactors(d)
    system.time(dw <- estimateGLMRobustDisp(d, design = design, prior.df = 10, maxit = 6))
    fw <- glmFit(dw, design = design, coef = 2)
    lrw <- glmLRT(fw, coef = 2)

    ## Extract results
    edger <- regions[nonzero]
    mcols(edger) <- cbind(mcols(edger), DataFrame(lrw$table))
    mcols(edger)$pvalue <-  lrw$table$PValue
    mcols(edger)$padj <- p.adjust(lrw$table$PValue, 'BH')

    ## Which are significant?
    mcols(edger)$sig <- mcols(edger)$padj < 0.05
    mcols(edger)$sig[is.na(mcols(edger)$sig)] <- FALSE

    ## Save results
    if(!is.null(file))
    save(edger, file = paste0(file, '-edgeR.Rdata'))
    
    ## End
    return(edger)
}

## Run DE analyses
message(paste(Sys.time(), 'running DESeq2 analysis'))
deseq <- run_deseq(counts, regions, file)

message(paste(Sys.time(), 'running edgeR analysis'))
edger <- run_edger(counts, regions, file)

## Print some summary info
print('DESeq2 summary - sig FDR 5%')
table(mcols(deseq)$padj < 0.05, useNA = 'ifany')

print('edgeR summary - sig FDR 5%')
table(mcols(edger)$padj < 0.05, useNA = 'ifany')

print('DEseq2 vs edgeR summary')
table(mcols(deseq)$padj < 0.05, mcols(edger)$padj < 0.05, dnn = list('DESeq2 sig FDR 5%', 'edgeR sig FDR 5%'))



## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
