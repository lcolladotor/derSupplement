## Load the data without a filter, save it, then filter it for derfinder processing steps

## Load libraries
library('getopt')

## Available at http://www.bioconductor.org/packages/release/bioc/html/derfinder.html
library('derfinder')
library('devtools')

## Specify parameters
spec <- matrix(c(
	'experiment', 'e', 1, 'character', 'Experiment. Only brainspan',
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
stopifnot(opt$experiment %in% c('brainspan'))

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
    ## Drop bad samples
    bad_samples <- which(rownames(pdSpan) %in% c('216', '218', '219'))
    pdSpan[bad_samples, ]
    pdSpan <- pdSpan[-bad_samples, ]
    stopifnot(nrow(pdSpan) == 484)
    

    ## Build the models
    fetal <- ifelse(pdSpan$Age < 0, "fetal", "adult")
    mod <- model.matrix(~ fetal * pdSpan$structure_acronym)
    mod0 <- model.matrix(~1, data=pdSpan)
    models <- list(mod=mod, mod0=mod0)

    ## Save information used for analyzeChr(groupInfo)
    # https://www.dropbox.com/s/nzv8r9rw7xi27vt/boxplots_brainspan.jpg
    # First 11 are neocortical, next four are not, last is cerebellum
    groupInfo <- factor(paste(ifelse(pdSpan$structure_acronym %in% c("DFC", "VFC", "MFC", "OFC", "M1C", "S1C", "IPC", "A1C", "STC", "ITC", "V1C"), "Neo", ifelse(pdSpan$structure_acronym %in% c("HIP", "AMY", "STR", "MD"), "notNeo", "CBC")), toupper(substr(fetal, 1, 1)), sep="."), levels=paste(rep(c("Neo", "notNeo", "CBC"), each=2), toupper(substr(unique(fetal), 1, 1)), sep="."))
    
}

## Save models
save(models, file="models.Rdata")

## Save information used for analyzeChr(groupInfo)
save(groupInfo, file="groupInfo.Rdata")

## Done :-)
proc.time()
options(width = 120)
session_info()
