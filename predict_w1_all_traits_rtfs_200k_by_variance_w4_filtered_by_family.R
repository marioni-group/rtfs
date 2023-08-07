library(MethylPipeR)

initLogs('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/results/continuous_traits_w4_filtered_by_family/', note = 'Lasso for all continuous traits W4 filtered by family. 5-fold cross-validation. Prediction on W1')

set.seed(42)

phenotypesW1 <- readRDS('~/rtfs_20k/phenotypesW1.rds')
methylW1 <- readRDS('~/rtfs_20k/methylW1.rds')

topCpGsByVariance <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/top200kCpGsByVarianceW4.rds')
methylW1 <- methylW1[, topCpGsByVariance]
gc()

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

traitResults <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/results/continuous_traits_w4_filtered_by_family/output_2022_10_10_11_44_58/traitResults.rds')

testTrait <- function(testMethyl, testPhenotypes, trait) {
  print(paste0('Removing rows with missing ', trait, ' from test data.'))
  testRemoveNAResult <- removeTraitNAs(testPhenotypes, list(testMethyl = testMethyl), trait)
  testPhenotypes <- testRemoveNAResult$traitDF
  testMethyl <- testRemoveNAResult$otherDFs$testMethyl

  print(paste0('Test model for ', trait, '.'))
  methylModel <- traitResults[[trait]]$model
  testPrediction <- predictMPRModel(methylModel, testMethyl, s = 'lambda.min')
  testDF <- data.frame(trait = testPhenotypes[[trait]], methyl = testPrediction)
  nullModel <- lm(trait ~ 1, data = testDF)
  fullModel <- lm(trait ~ testPrediction, data = testDF)
  list(trait = trait, model = methylModel, testPrediction = testPrediction, testDF = testDF, nullModel = nullModel, fullModel = fullModel)
}

testResults <- lapply(traits, function(trait) {
  testTrait(methylW1, phenotypesW1, trait)
})

names(testResults) <- traits

summaries <- lapply(testResults, function(testResult) {
  testResultSummary <- summary(testResult$fullModel)
  list(coefficients = testResultSummary$coefficients, adj.r.squared = testResultSummary$adj.r.squared)
})

gc()
