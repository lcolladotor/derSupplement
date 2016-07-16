## Usage:
# qrsh -l mem_free=130G,h_vmem=150G
# mkdir -p logs
# module load R/3.3
# Rscript region_cuts_plot.R > logs/region_cuts_plot.txt 2>&1
library('GenomicRanges')
library('ggplot2')
library('reshape2')

## Load data
load('region_cuts.Rdata')
load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
gs <- GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome

gs_exon <- gs[gs$theRegion == 'exon']
stopifnot(all(countOverlaps(gs_exon) == 1))

## Extract information of interest
message(paste(Sys.time(), 'extracting information of interest'))
info <- lapply(region_cuts, function(regs) {
    regs6 <- regs[width(regs) >= 6]
    ## Find overlap with known exons
    ov_ex_reg <- round(
        mean(countOverlaps(gs_exon, regs6) > 0) * 100, 2)
    ov_reg_ex <- round(
        mean(countOverlaps(regs6, gs_exon) > 0) * 100, 2)
    res <- list(
        mean = mean(width(regs)),
        mean6 = mean(width(regs6)),
        n = length(regs),
        n6 = length(regs6),
        quantile = quantile(width(regs), seq(0, 1, by = 0.1)),
        quantile6 = quantile(width(regs6), seq(0, 1, by = 0.1)),
        sd = sd(width(regs)),
        sd6 = sd(width(regs6)),
        ov_ex_reg = ov_ex_reg,
        ov_reg_ex = ov_reg_ex
    )
    return(res)
})

message(paste(Sys.time(), 'creating summary data.frame'))
n <- unlist(lapply(info, '[[', 'n'))
n6 <- unlist(lapply(info, '[[', 'n6'))
mean <- unlist(lapply(info, '[[', 'mean'))
mean6 <- unlist(lapply(info, '[[', 'mean6'))
sd <- unlist(lapply(info, '[[', 'sd'))
sd6 <- unlist(lapply(info, '[[', 'sd6'))
quantile <- do.call(rbind, lapply(info, '[[', 'quantile'))
quantile6 <- do.call(rbind, lapply(info, '[[', 'quantile6'))
ov_ex_reg <- unlist(lapply(info, '[[', 'ov_ex_reg'))
ov_reg_ex <- unlist(lapply(info, '[[', 'ov_reg_ex'))

## Arrange information into a single data frame
regInfo <- data.frame(
    cutoff = as.numeric(names(region_cuts)),
    n = n,
    n6 = n6,
    mean = mean,
    mean6 = mean6,
    sd = sd,
    sd6 = sd6,
    ov_ex_reg = ov_ex_reg,
    ov_reg_ex = ov_reg_ex,
    stringsAsFactors = FALSE)
regInfo <- cbind(regInfo, quantile, quantile6)
save(regInfo, file = 'regInfo.Rdata')


pdf(file = 'region_cuts.pdf')
ggplot(data = regInfo, aes(x = cutoff, y = n)) + geom_point() + ylab('Number of ERs (all)') + xlab('Cutoff') + geom_line() + theme_linedraw(base_size = 16)
ggplot(data = regInfo, aes(x = cutoff, y = n6)) + geom_point() + geom_line() + ylab('Number of ERs') + xlab('Cutoff') + theme_linedraw(base_size = 16)


ggplot(data = regInfo, aes(x = cutoff, y = ov_ex_reg)) + geom_point() + ylab('Percent of ENSEMBL exons overlapping at least one ER') + xlab('Cutoff') + geom_line() + theme_linedraw(base_size = 16) + ylim(25, 60)
ggplot(data = regInfo, aes(x = cutoff, y = ov_reg_ex)) + geom_point() + ylab('Percent of ERs overlapping at least one ENSEMBL exon') + xlab('Cutoff') + geom_line() + theme_linedraw(base_size = 16) + ylim(c(60, 100))

df <- melt(regInfo[, c(1, 21:31)], id = 'cutoff')
ggplot(data = df, aes(x = cutoff, y = log10(value + 1), color = variable)) + geom_line() + theme_linedraw(base_size = 16) + xlab('Cutoff') + ylab('ER width (log10 + 1)') + scale_colour_discrete(name = 'Quantile') + geom_point()
dev.off()

## Reproducibility info
Sys.time()
proc.time()
options(width = 120)
session_info()


