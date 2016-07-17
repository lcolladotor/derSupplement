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

## Error control options
adj.method <- 'BH' ## Can easily change between 'BH' and 'bonferroni'
adj.acron <- ifelse(adj.method == 'BH', 'FDR', 'FWER')

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
    
} else if (opt$pipeline == 'railMatrix') {
    load(file.path('..', 'railMatrix', paste0('regionMat-R', opt$replicate, '.Rdata')))
    
    counts <- regionMat$chr17$coverageMatrix
    regions <- regionMat$chr17$regions
    file <- paste0('stats-railMatrix-R', opt$replicate)
    
} else {
    stop("Invalid 'pipeline' argument")
}





## DESeq2 analysis
run_deseq <- function(counts, regions, file = NULL, use_zero = FALSE) {
    counts <- round(counts, 0)
    groupInfo <- as.factor(as.integer(gsub('.*G|R.*', '', colnames(counts))))
    if(use_zero) {
        nonzero <- TRUE
    } else {
        nonzero <- sapply(rowSums(counts), function(x) {x > 0})
    }
    
    
    ## Round matrix and specify design
    library('DESeq2')
    dse <- DESeqDataSetFromMatrix(counts[nonzero, ], data.frame(group = groupInfo), ~ group)

    ## Perform DE analysis
    dse <- DESeq(dse)

    ## Extract results
    deseq <- regions[nonzero]
    mcols(deseq) <- cbind(mcols(deseq), results(dse, alpha = 0.05, pAdjustMethod = adj.method))

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
run_edger <- function(counts, regions, file = NULL, use_zero = FALSE) {
    counts <- round(counts, 0)
    groupInfo <- as.factor(as.integer(gsub('.*G|R.*', '', colnames(counts))))
    if(use_zero) {
        nonzero <- TRUE
    } else {
        nonzero <- sapply(rowSums(counts), function(x) {x > 0})
    }
    
    ## Determine design matrix
    design <- model.matrix(~ groupInfo)

    ## Perform DE analysis
    library('edgeR')
    d <- DGEList(counts = counts[nonzero, ], group = groupInfo)
    d <- calcNormFactors(d)
    dw <- estimateGLMRobustDisp(d, design = design)
    fw <- glmFit(dw, design = design, coef = 2)
    lrw <- glmLRT(fw, coef = 2)

    ## Extract results
    edger <- regions[nonzero]
    mcols(edger) <- cbind(mcols(edger), DataFrame(lrw$table))
    mcols(edger)$pvalue <- lrw$table$PValue
    mcols(edger)$padj <- p.adjust(lrw$table$PValue, adj.method)

    ## Which are significant?
    mcols(edger)$sig <- mcols(edger)$padj < 0.05
    mcols(edger)$sig[is.na(mcols(edger)$sig)] <- FALSE

    ## Save results
    if(!is.null(file))
    save(edger, file = paste0(file, '-edgeR.Rdata'))
    
    ## End
    return(edger)
}

## limma analysis
# get the f statistic from 2 lmFit objects
getF <- function(fit, fit0, theData) {
	
	rss1 = rowSums((fitted(fit)-theData)^2)
	df1 = ncol(fit$coef)
	rss0 = rowSums((fitted(fit0)-theData)^2)
	df0 = ncol(fit0$coef)

	fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
	f_pval = pf(fstat, df1-1, ncol(theData)-df1,lower.tail=FALSE)
	fout = cbind(fstat,df1-1,ncol(theData)-df1,f_pval)
	colnames(fout)[2:3] = c("df1", "df0")
	fout = data.frame(fout)
	return(fout)
}
run_limma <- function(counts, regions, file = NULL, use_zero = FALSE) {
    counts <- round(counts, 0)
    groupInfo <- as.factor(as.integer(gsub('.*G|R.*', '', colnames(counts))))
    if(use_zero) {
        nonzero <- TRUE
    } else {
        nonzero <- sapply(rowSums(counts), function(x) {x > 0})
    }
    
    ## Perform DE analysis
    y <- log2(counts[nonzero, ] + 1)
    mod <-  model.matrix(~ groupInfo)
    mod0 <- model.matrix(~ 1, data = groupInfo)
    library('limma')
    fit <- lmFit(y, mod)
    eb <- ebayes(fit)
    fit0 <- lmFit(y, mod0)
    ff <- getF(fit, fit0, y)
    
    ## Extract results
    limma <- regions[nonzero]
    mcols(limma) <- cbind(mcols(limma), DataFrame(ff))
    mcols(limma)$pvalue <- ff$f_pval
    mcols(limma)$padj <- p.adjust(ff$pvalue, adj.method)
    
    ## Which are significant?
    mcols(limma)$sig <- mcols(limma)$padj < 0.05
    mcols(limma)$sig[is.na(mcols(limma)$sig)] <- FALSE

    ## Save results
    if(!is.null(file))
    save(limma, file = paste0(file, '-limma.Rdata'))
    
    ## End
    return(limma)
}

## Run DE analyses
message(paste(Sys.time(), 'running DESeq2 analysis'))
deseq <- run_deseq(counts, regions, file, TRUE)

message(paste(Sys.time(), 'running edgeR analysis'))
edger <- run_edger(counts, regions, file, TRUE)

message(paste(Sys.time(), 'running limma analysis'))
limma <- run_limma(counts, regions, file, TRUE)

## Print some summary info
print(paste('DESeq2 summary - sig', adj.acron, '5%'))
table(mcols(deseq)$padj < 0.05, useNA = 'ifany')

print(paste('edgeR summary - sig', adj.acron, '5%'))
table(mcols(edger)$padj < 0.05, useNA = 'ifany')

print(paste('limma summary - sig', adj.acron, '5%'))
table(mcols(limma)$padj < 0.05, useNA = 'ifany')

print('DESeq2 vs edgeR summary')
addmargins(table(mcols(deseq)$padj < 0.05, mcols(edger)$padj < 0.05, dnn = list(paste('DESeq2 sig', adj.acron, '5%'), paste('edgeR sig', adj.acron, '5%'))))

print('DESeq2 vs limma summary')
addmargins(table(mcols(deseq)$padj < 0.05, mcols(limma)$padj < 0.05, dnn = list(paste('DESeq2 sig', adj.acron, '5%'), paste('limma sig', adj.acron, '5%'))))

print('edgeR vs limma summary')
addmargins(table(mcols(edgeR)$padj < 0.05, mcols(limma)$padj < 0.05, dnn = list(paste('edgeR sig', adj.acron, '5%'), paste('limma sig', adj.acron, '5%'))))


## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()
