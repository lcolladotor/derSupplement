## Remake the region width plot in higher quality
# module load R/3.3
# mkdir -p logs
# Rscript regionWidth.R -p "sra" > logs/regionWidth_log.txt 2>&1


library('GenomicRanges')
library('ggplot2')
library('grid')
library('gridExtra')

load('fullRegions.Rdata')

## For ggplot
tmp <- fullRegions
names(tmp) <- seq_len(length(tmp))
regions.df <- as.data.frame(tmp)
regions.df$width <- width(tmp)
rm(tmp)

## Find which chrs are present in the data set
chrs <- levels(seqnames(fullRegions))

regions.df.sig <- regions.df[regions.df$significantFWER == 'TRUE', ]



theme_set(theme_bw(base_size = 15))

## Special subsets: need at least 3 points for a density plot
keepChr <- table(regions.df$seqnames) > 2
regions.df.plot <- subset(regions.df, seqnames %in% names(keepChr[keepChr]))

pdf('regionWidth.pdf', width=10, height=10)
xrange <- range(log10(regions.df.plot$width)) * c(0.95, 1.05)
p2a <- ggplot(regions.df.plot, aes(x=log10(width), colour=seqnames)) + 
    geom_line(stat='density') + labs(title='Density of region lengths') +
    xlab('Region width (log10)') + scale_colour_discrete(limits=chrs) +
    xlim(xrange) + theme(legend.position="none")
p2b <- ggplot(regions.df.sig, aes(x=log10(width), colour=seqnames)) +
    geom_line(stat='density') +
    labs(title='Density of region lengths (significant only)') +
    xlab('Region width (log10)') + scale_colour_discrete(limits=chrs) +
    xlim(xrange) + theme(legend.position="none")
grid.arrange(p2a, p2b)
dev.off()

## Reproducibility info
proc.time()
message(Sys.time())
options(width = 120)
devtools::session_info()
