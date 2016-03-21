library('ballgown')
sampleDirs <- dir('../stringtie', pattern = 'sample.*R1')
#sampleDirs <- sampleDirs[!grepl('assembly', sampleDirs)]
sampleDirs <- sampleDirs[grepl('assembly', sampleDirs)]
pData <- data.frame(sampleDir = sampleDirs, sample = as.integer(gsub('sample|G.*', '', sampleDirs)), group = as.integer(gsub('sample.*G|R.*', '', sampleDirs)), replicate = 1)
bg <- ballgown(samples = sampleDirs, pData = pData)

stat_results <- stattest(gown = bg, feature = 'gene', getFC = TRUE, meas = 'FPKM', covariate = 'group')
stat_results2 <- stattest(gown = bg, feature = 'trans', getFC = TRUE, meas = 'FPKM', covariate = 'group')

## Gene-level FC
summary(stat_results$fc)
## Transcript-level FC
summary(stat_results2$fc)

table(stat_results$fc < 2)
table(stat_results$fc < 2) / length(stat_results$fc) * 100

table(stat_results2$fc < 2)
table(stat_results2$fc < 2) / length(stat_results2$fc) * 100

table(stat_results$fc < 2 & stat_results$fc > 0.5)
table(stat_results$fc < 2 & stat_results$fc > 0.5) / length(stat_results$fc) * 100

table(stat_results2$fc < 2 & stat_results2$fc > 0.5)
table(stat_results2$fc < 2 & stat_results2$fc > 0.5) / length(stat_results2$fc) * 100

dim(stat_results)
dim(stat_results2)

## Load cuffcompare data
#ccomp <- read.table('cuffcompare/cuffcomp-R1.tracking', sep = '\t', header = FALSE, stringsAsFactors = FALSE, col.names = c('tcons', 'gene_id', 'gene_name', 'code', 'compact'))
ccomp <- read.table('sample1G1R1-no-assembly/cuffcomp.tracking', sep = '\t', header = FALSE, stringsAsFactors = FALSE, col.names = c('tcons', 'gene_id', 'gene_name', 'code', 'compact'))
dim(ccomp)
head(ccomp)
#ccomp$t_name <- gsub('q1:', '', sapply(strsplit(ccomp$compact, '\\|'), '[[', 2))
ccomp$t_name <- sapply(strsplit(ccomp$compact, '\\|'), '[[', 2)

## Load transcript names
#t_data <- read.table('sample1G1R1/t_data.ctab', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
t_data <- read.table('sample1G1R1-no-assembly/t_data.ctab', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
dim(t_data)
head(t_data)

## Match the tables
i <- match(t_data$t_name, ccomp$t_name)
length(i)
sum(is.na(i))
nrow(t_data) - nrow(ccomp)

## Fold change within 0.5 and 2 (test at the transcript-level)
normal_t <- stat_results2$fc < 2 & stat_results2$fc > 0.5
sum(normal_t)
round(sum(normal_t) / length(normal_t) * 100, 2)

## Transcripts by cuffcompare code (all)
table(ccomp$code[i], useNA = 'ifany')
round(table(ccomp$code[i], useNA = 'ifany') / length(ccomp$code[i]) * 100, 2)

## Transcripts by cuffcompare code (fold change within 0.5 and 2)
table(ccomp$code[i[normal_t]], useNA = 'ifany')
round(table(ccomp$code[i[normal_t]], useNA = 'ifany') / sum(normal_t) * 100, 2)

## Transcripts by cuffcompare code (fold change NOT within 0.5 and 2)
table(ccomp$code[i[!normal_t]], useNA = 'ifany')
round(table(ccomp$code[i[!normal_t]], useNA = 'ifany') / sum(!normal_t) * 100, 2)

