config <- yaml::read_yaml(here::here("config.yml"))

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  run_index <- 1
  seed <- 42
#   stop('Requires two arguments. 1: the run index. 2: the random seed')
} else {
  run_index <- as.numeric(args[[1]])
  seed <- as.numeric(args[[2]])
  print(paste0('run_index: ', run_index))
  print(paste0('seed: ', seed))
}

# Wait for 2 * run_index seconds to guarantee runs have a unique timestamp

# Sys.sleep(2 * run_index)

library(MethylPipeR)
# stop()
initLogs(config$methylpiper_logs_path, note = 'Lasso for all continuous traits. 5-fold cross-validation.')


set.seed(seed)

print(paste0('run_index = ', run_index, ', seed = ', seed))

phenotypesW4 <- readRDS(paste0(config$continuous_traits_data_path, 'phenotypesW4.rds'))
methylW4 <- readRDS(paste0(config$continuous_traits_data_path, 'methylW4.rds'))

# Get 450k IDs from annotation file
get450kIDs <- function(filepath = config$epic_annotation_file_path) {
  annotations <- readRDS(filepath)
  rownames(annotations)[annotations[, 'Methyl450_Loci'] == 'TRUE']
}

print('Calculating 450k and methyl column name intersection.')
subset450kCpGs <- intersect(colnames(methylW4), get450kIDs())
print('Subset methylW4 to 450k subset.')
methylW4 <- methylW4[, subset450kCpGs]
gc()

getFilterByVarianceIDs <- function(data, numberOfFeatures) {
  vars <- apply(data, 2, var)
  sortedVars <- sort(vars, index.return = TRUE, decreasing = TRUE)
  cpgIDs <- colnames(data)[sortedVars$ix]
  cpgIDs[1:numberOfFeatures]
}

p <- 200000
topCpGsByVariance <- getFilterByVarianceIDs(methylW4, p)
methylW4 <- methylW4[, topCpGsByVariance]
gc()

# Get family IDs
pedigree <- read.csv(config$pedigree_file_path)
phenotypesW4 <- merge(phenotypesW4, pedigree[, c('volid', 'famid')], by.x = 'ID', by.y = 'volid', all.x = TRUE)

methylW4Index <- match(phenotypesW4$Sample_Sentrix_ID, row.names(methylW4))
methylW4 <- methylW4[methylW4Index, ]

# Get W1 family IDs to filter W4
phenotypesW1 <- readRDS(paste0(config$continuous_traits_data_path, 'phenotypesW1.rds'))
phenotypesW1 <- merge(phenotypesW1, pedigree[, c('volid', 'famid')], by.x = 'ID', by.y = 'volid', all.x = TRUE)

familyToRemoveW4 <- phenotypesW4$famid %in% phenotypesW1$famid
phenotypesW4 <- phenotypesW4[!familyToRemoveW4, ]
methylW4 <- methylW4[!familyToRemoveW4, ]


traits <- c('age', 'glucose', 'cholest', 'HDL', 'sodium', 'potassium', 'urea', 'creatinine', 'bmi', 'whr', 'fat', 'sBP', 'dBP', 'HR', 'FEV', 'FVC', 'Alc', 'PckYrs', 'g')

removeTraitNAs <- function(traitDF, otherDFs, trait) {
  rowsToKeep <- !is.na(traitDF[[trait]])
  traitDF <- traitDF[rowsToKeep, ]
  otherDFs <- lapply(otherDFs, function(df) {
    if (is.data.frame(df) || is.matrix(df)) {
      df[rowsToKeep, ]
    } else if (is.null(df)) {
      # For example, if foldID is NULL in cvTrait
      df
    } else {
      # Assumes df is a vector
      df[rowsToKeep]
    }
  })
  list(traitDF = traitDF, otherDFs = otherDFs)
}

cvTrait <- function(trainMethyl, trainPhenotypes, trait, nFolds, foldID = NULL) {
  print(paste0('Removing rows with missing ', trait, ' from training data.'))
  trainRemoveNAResult <- removeTraitNAs(trainPhenotypes, list(trainMethyl = trainMethyl, foldID = foldID), trait)
  trainPhenotypes <- trainRemoveNAResult$traitDF
  trainMethyl <- trainRemoveNAResult$otherDFs$trainMethyl
  foldID <- trainRemoveNAResult$otherDFs$foldID

  if (is.null(foldID)) {
    foldID <- getGroupCVFoldIDs(trainPhenotypes$famid, nFolds)$foldIDs
  }

  print('Fitting lasso model')
  methylModel <- fitMPRModelCV(type = 'continuous',
    method = 'glmnet',
    trainXs = trainMethyl,
    trainY = trainPhenotypes[[trait]],
    seed = 42,
    alpha = 1,
    nFolds = nFolds,
    foldID = foldID,
    # parallel = TRUE,
    trace.it = 1)
  list(trait = trait, model = methylModel, foldID = foldID)
}

traitResults <- lapply(traits, function(trait) {
  cvTrait(methylW4, phenotypesW4, trait, 5)
})

names(traitResults) <- traits

sessionStartTimestamp <- getOption('mprSessionStartTimestamp')
sessionLogFolder <- getOption('mprSessionLogFolder')

folderPath <- paste0(sessionLogFolder, 'output_', sessionStartTimestamp, '/')
dir.create(paste0(folderPath))
filePath <- paste0(folderPath, 'traitResults.rds')
saveRDS(traitResults, file = filePath)
