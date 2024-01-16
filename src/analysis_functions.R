library(dplyr)

getFilterByVarianceIDs <- function(data, numberOfFeatures) {
  vars <- apply(data, 2, var)
  sortedVars <- sort(vars, index.return = TRUE, decreasing = TRUE)
  cpgIDs <- colnames(data)[sortedVars$ix]
  cpgIDs[1:numberOfFeatures]
}

# Get 450k IDs from annotation file
get450kIDs <- function(filepath = "/Cluster_Filespace/Marioni_Group/Yipeng/data/GS/IHD_pipeline/EPIC_AnnotationObject_df.rds") {
  annotations <- readRDS(filepath)
  rownames(annotations)[annotations[, 'Methyl450_Loci'] == 'TRUE']
}

predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
  uniqueTimes <- sort(unique(dataDF$time_to_event))
  thresholdIndex <- match(threshold, uniqueTimes)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$time_to_event, dataDF$Event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
  onsetPredictions <- 1 - survivalPredictions

  # Event should be 0 if tte is > 10
  dataDF$Event <- sapply(1:nrow(dataDF), function(i) {
    if (dataDF$time_to_event[[i]] > threshold) {
      0
    } else {
      dataDF$Event[[i]]
    }
  })

  auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$Event)
  roc <- pROC::roc(response = dataDF$Event, predictor = onsetPredictions)
  list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
}

# Censoring should be "apr_2022" or "oct_2020" or NULL (original phenotype - oct_2020 censoring date but no TTE rounding)
load450kW3W1 <- function(censoring = NULL) {
  if (!is.null(censoring)) {
    # targetW3 <- readRDS(paste0('/Local_Scratch/Yipeng/w3CoxTable_', censoring, '.rds'))
    # methylW3 <- readRDS('/Local_Scratch/Yipeng/w3Methyl_apr_2022.rds')

    # targetW1 <- readRDS(paste0('/Local_Scratch/Yipeng/w1CoxTable_', censoring, '.rds'))
    # methylW1 <- readRDS('/Local_Scratch/Yipeng/w1Methyl_apr_2022.rds')

    targetW3 <- readRDS(paste0("/Local_Scratch/Yipeng/w3CoxTable_gp_smr_", censoring, ".rds"))
    targetW1 <- readRDS(paste0("/Local_Scratch/Yipeng/w1CoxTable_gp_smr_", censoring, ".rds"))

    methylW3 <- readRDS(paste0("/Local_Scratch/Yipeng/w3Methyl_gp_smr_", censoring, ".rds"))
    methylW1 <- readRDS(paste0("/Local_Scratch/Yipeng/w1Methyl_gp_smr_", censoring, ".rds"))
  } else {
    # Stop script as this option should not be used
    # stop()
    targetW3 <- readRDS('/Local_Scratch/Yipeng/w3CoxTable.rds')
    methylW3 <- readRDS('/Local_Scratch/Yipeng/w3Methyl.rds')

    targetW1 <- readRDS('/Local_Scratch/Yipeng/w1CoxTable.rds')
    methylW1 <- readRDS('/Local_Scratch/Yipeng/w1Methyl.rds')
  }
  # Filter out NA BMI rows from targetW1
  targetW1NonNABMIIndex <- !is.na(targetW1$bmi)
  methylW1 <- methylW1[targetW1NonNABMIIndex, ]
  targetW1 <- targetW1[targetW1NonNABMIIndex, ]
  row.names(targetW1) <- NULL
  gc()

  print('Calculating 450k and methyl column name intersection.')
  subset450kCpGs <- intersect(colnames(methylW3), get450kIDs())
  print('Subset to 450k subset.')
  methylW3 <- methylW3[, subset450kCpGs]
  gc()
  methylW1 <- methylW1[, subset450kCpGs]
  gc()

  list(targetW3 = targetW3, methylW3 = methylW3, targetW1 = targetW1, methylW1 = methylW1)
}

addFamilyInformation <- function(target, pedigreeFile = '/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_diabetes_pipeline/using_methylpiper/src/20k/preprocessing/2022-01-17_pedigree.csv') {
  pedigree <- read.csv(pedigreeFile)
  target <- merge(target, pedigree[, c('volid', 'famid')], by.x = 'Sample_Name', by.y = 'volid', all.x = TRUE)
}

addPRS <- function(target, prsFile = '/Cluster_Filespace/Marioni_Group/Yipeng/INTERVENE/flagship/prs_output/T2D_PRS.sscore', scoreCol = 'SCORE1_SUM', sep = '\t', targetIDCol = 'Sample_Name', prsIDCol = 'IID') {
  prs <- read.csv(prsFile, sep = sep)
  target <- target %>% left_join(prs[, c(prsIDCol, scoreCol)], by = setNames(prsIDCol, targetIDCol))
}

cvAlpha <- function(modelType, modelMethod, methyl, target, seed, nFolds, parallel, standardize, searchAlphas) {
  models <- lapply(searchAlphas, function(alpha) {
    fitMPRModelCV('survival', 'glmnet', methyl, target[, c('Event', 'time_to_event')], seed = seed, nFolds = nFolds, parallel = TRUE, standardize = standardize, alpha = alpha)
  })
  # Cross-validated mean errors at lambda.min
  cvms <- sapply(models, function(model) {
    model$model$cvm[[model$model$index['min',]]]
  })
  minCVM <- min(cvms)
  bestAlphaIndex <- match(minCVM, cvms)
  bestAlpha <- searchAlphas[[bestAlphaIndex]]
  bestModel <- models[[bestAlphaIndex]]
  list(bestAlpha = bestAlpha, bestModel = bestModel, models = models)
}

fitAndPredict <- function(trainMethyl, trainTarget, testMethyl, testTarget, alpha = NULL, seed = 42, nFolds = 9, standardize = FALSE, searchAlphas = NULL) {
  if (!is.null(alpha)) {
    cvAlphaResult <- NULL
    model <- fitMPRModelCV('survival', 'glmnet', trainMethyl, trainTarget[, c('Event', 'time_to_event')], seed = seed, nFolds = nFolds, parallel = TRUE, standardize = standardize, alpha = alpha)
  } else { # If alpha == NULL, optimise using k-fold cross-validation
    cvAlphaResult <- cvAlpha('survival', 'glmnet', trainMethyl, trainTarget[, c('Event', 'time_to_event')], seed = seed, nFolds = nFolds, parallel = TRUE, standardize = standardize, searchAlphas = searchAlphas)
    model <- cvAlphaResult$bestModel
  }

  prediction <- predictMPRModel(model, testMethyl, s = 'lambda.min')

  testTarget$dScore <- prediction

  riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes, testTarget)
  riskFactorsPRSCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + SCORE1_SUM, testTarget)
  directCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + dScore, testTarget)
  directPRSCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + SCORE1_SUM + dScore, testTarget)

  models <- list(r = riskFactorsOnlyCoxPH,
                 rPRS = riskFactorsPRSCoxPH,
                 d = directCoxPH,
                 dPRS = directPRSCoxPH)

  predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
    uniqueTimes <- sort(unique(c(dataDF$time_to_event, threshold)))
    thresholdIndex <- match(threshold, uniqueTimes)
    cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$time_to_event, dataDF$Event, predict(coxPHModel), uniqueTimes)
    survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
    onsetPredictions <- 1 - survivalPredictions

    # Event should be 0 if tte is > 10
    dataDF$Event <- sapply(1:nrow(dataDF), function(i) {
      if (dataDF$time_to_event[[i]] > threshold) {
        0
      } else {
        dataDF$Event[[i]]
      }
    })

    auc <- MLmetrics::AUC(y_pred = onsetPredictions, y_true = dataDF$Event)
    prauc <- MLmetrics::PRAUC(y_pred = onsetPredictions, y_true = dataDF$Event)
    roc <- pROC::roc(response = dataDF$Event, predictor = onsetPredictions)
    list(cumulativeBaseHaz = cumulativeBaseHaz, onsetPredictions = onsetPredictions, auc = auc, prauc = prauc, roc = roc)
  }
  testResults <- lapply(models, function(m) {predictCoxPHOnset(testTarget, m)})

  aucs <- sapply(testResults, function(r) {r$auc})
  praucs <- sapply(testResults, function(r) {r$prauc})
  metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)

  row.names(metricsTable) <- c('Risk factors only',
                               'Risk factors + PRS',
                               'Risk factors + direct EpiScore',
                               'Risk factors + PRS + direct EpiScore')

  list(testResults = testResults, metricsTable = metricsTable, testTarget = testTarget, testMethyl = testMethyl, models = models, methylModel = model, alpha = alpha, cvAlphaResult = cvAlphaResult)
}


# Saves each object in a named list of objects with filenames corresponding to those names.
# Files are saved in a folder named 'output_<sessionStartTimestamp>' in the MethylPipeR session log folder.
saveResults <- function(resultsList) {
  sessionStartTimestamp <- getOption('mprSessionStartTimestamp')
  sessionLogFolder <- getOption('mprSessionLogFolder')
  folderPath <- paste0(sessionLogFolder, 'output_', sessionStartTimestamp, '/')

  dir.create(paste0(folderPath))
  

  sapply(names(resultsList), function(x) {
    print(paste0('Saving ', x))
    saveRDS(resultsList[[x]], paste0(folderPath, x, '.rds'))
  })
  
}


