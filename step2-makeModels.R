## Load the data without a filter, save it, then filter it for derfinder processing steps

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'experiment', 'e', 1, 'character', 'Experiment. Either brainspan, snyder or hippo',
	'help' , 'h', 0, 'logical', 'Display help'
), byrow=TRUE, ncol=5)
opt <- getopt(spec)


## if help was asked for print a friendly message
## and exit with a non-zero error code
if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

## Check experiment input
stopifnot(opt$experiment %in% c('brainspan', 'snyder', 'hippo'))

if(opt$experiment != 'brainspan') {
    ## Load the coverage information
    load(file.path('..', '..', 'CoverageInfo', 'fullCov.Rdata'))
    load(file.path('..', '..', 'CoverageInfo', 'chr22CovInfo.Rdata'))

    ## Identify the samplefiles
    files <- colnames(chr22CovInfo$coverage)
}

 ## Calculate the library adjustments and build the models
buildModels <- function(fullCov, testvars, colsubset = NULL) {
    ## Determine sample size adjustments
    if(file.exists("sampleDepths.Rdata")) {
    	load("sampleDepths.Rdata")
    } else {
    	if(file.exists("collapsedFull.Rdata")) {
    		load("collapsedFull.Rdata")
    	} else {
    		## Collapse
    		collapsedFull <- collapseFullCoverage(fullCov, colsubset = colsubset, save=TRUE)
    	}

    	## Get the adjustments
    	sampleDepths <- sampleDepth(collapsedFull = collapsedFull, probs = 1,
            nonzero = TRUE, scalefac = 32, center = FALSE)
    	save(sampleDepths, file="sampleDepths.Rdata")
    }
    ## Build the models
    models <- makeModels(sampleDepths = sampleDepths, testvars = testvars,
        adjustvars = NULL, testIntercept = FALSE)
    
    return(models)
}


##### Note that this whole section is for defining the models using makeModels()
##### You can alternatively define them manually and/or use packages such as splines if needed.

if(opt$experiment == 'brainspan') {
    ## Define the groups
    load("/home/epi/ajaffe/Lieber/Projects/Grants/Coverage_R01/brainspan/brainspan_phenotype.rda")

    ## Build the models
    fetal <- ifelse(pdSpan$Age < 0, "fetal", "adult")
    mod <- model.matrix(~ fetal * pdSpan$structure_acronym)
    mod0 <- model.matrix(~1, data=pdSpan)
    models <- list(mod=mod, mod0=mod0)

    ## Save information used for analyzeChr(groupInfo)
    # https://www.dropbox.com/s/nzv8r9rw7xi27vt/boxplots_brainspan.jpg
    # First 11 are neocortical, next four are not, last is cerebellum
    groupInfo <- factor(paste(ifelse(pdSpan$structure_acronym %in% c("DFC", "VFC", "MFC", "OFC", "M1C", "S1C", "IPC", "A1C", "STC", "ITC", "V1C"), "Neo", ifelse(pdSpan$structure_acronym %in% c("HIP", "AMY", "STR", "MD"), "notNeo", "CBC")), toupper(substr(fetal, 1, 1)), sep="."), levels=paste(rep(c("Neo", "notNeo", "CBC"), each=2), toupper(substr(unique(fetal), 1, 1)), sep="."))
    
} else if(opt$experiment == 'snyder') {
    ## The information table
    info <- read.csv("/home/epi/ajaffe/Lieber/Projects/Timecourse_RNAseq/Profile/phenotype.csv")
    info$shortBAM <- gsub(".*/", "", info$bamFile)

    ## Match dirs with actual rows in the info table
    match <- sapply(files, function(x) { which(info$GEO_ID == x)})
    info <- info[match, ]

    ## Test a spline on days
    library("splines")
    testvars <- bs(info$Day, df=5)

    ## Define the groups for plotting
    library("Hmisc")
    groupInfo <- cut2(info$Day, g=4)
    tmp <- groupInfo
    groupInfo <- factor(gsub(",", "to", gsub("\\[| |)|\\]", "", groupInfo)))
    names(groupInfo) <- tmp
    
    ## Build models
    models <- buildModels(fullCov, testvars)
} else if(opt$experiment == 'hippo') {
    ## Define the groups
    load("/home/epi/ajaffe/Lieber/Projects/RNAseq/HippoPublic/sra_phenotype_file.rda")
    info <- sra
    info <- info[complete.cases(info),]
    ## Match dirs with actual rows in the info table
    match <- sapply(files, function(x) { which(info$SampleID == x)})
    info <- info[match, ]
    ## Set the control group as the reference
    groupInfo <- factor(info$Pheno, levels=c("CT", "CO", "ETOH"))

    ## Define colsubset
    colsubset <- which(!is.na(groupInfo))
    save(colsubset, file="colsubset.Rdata")
    
    ## Update the group labels
    groupInfo <- groupInfo[!is.na(groupInfo)]
    
    ## Build models
    models <- buildModels(fullCov, groupInfo, colsubset)
}

## Save models
save(models, file="models.Rdata")

## Save information used for analyzeChr(groupInfo)
save(groupInfo, file="groupInfo.Rdata")

## Done :-)
proc.time()
options(width = 90)
session_info()
