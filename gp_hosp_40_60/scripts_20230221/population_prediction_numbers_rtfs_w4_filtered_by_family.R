library(MethylPipeR)


w1Target <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_05_11_00_44_56/testTarget.rds')

# Event should be 0 if tte is > 10
  w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
    if (w1Target$time_to_event[[i]] > 10) {
      0
    } else {
      w1Target$Event[[i]]
    }
  })

w1Target$Event <- as.factor(w1Target$Event)

coxTestResults <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_05_11_00_44_56/testResults.rds')

nullResponse <- coxTestResults$rPRS$onsetPredictions
coxLassoIncrementalResponse <- coxTestResults$dPRS$onsetPredictions


rocResult <- pROC::roc(w1Target$Event, coxLassoIncrementalResponse)

thresholds <- seq(from = 0, to = 1, by = 0.1)

fullBinaryPredictions <- lapply(thresholds, function(threshold) {
  as.factor(as.numeric(coxLassoIncrementalResponse >= threshold))
})

names(fullBinaryPredictions) <- thresholds

nullBinaryPredictions <- lapply(thresholds, function(threshold) {
  as.factor(as.numeric(nullResponse >= threshold))
})

names(nullBinaryPredictions) <- thresholds

fullConfusionMatrices <- lapply(fullBinaryPredictions, function(x) {
  caret::confusionMatrix(x, w1Target$Event, positive = '1')
})

nullConfusionMatrices <- lapply(nullBinaryPredictions, function(x) {
  caret::confusionMatrix(x, w1Target$Event, positive = '1')
})

fullMetrics <- do.call(rbind, lapply(fullConfusionMatrices, function(confusionMatrix) {confusionMatrix[[4]][1:4]}))
fm <- fullMetrics

colnames(fullMetrics) <- c('RTFS Sensitivity', 'RTFS Specificity', 'RTFS PPV', 'RTFS NPV')

nullMetrics <- do.call(rbind, lapply(nullConfusionMatrices, function(confusionMatrix) {confusionMatrix[[4]][1:4]}))
colnames(nullMetrics) <- c('Null Sensitivity', 'Null Specificity', 'Null PPV', 'Null NPV')

metrics <- as.data.frame(cbind(fullMetrics, nullMetrics))


prevalence <- 0.33

n <- 10000

nCases <- floor(n * prevalence)
nControls <- n - nCases

metrics$diffTP <- metrics[, 'RTFS Sensitivity'] * nCases - metrics[, 'Null Sensitivity'] * nCases
metrics$diffFN <- (1 - metrics[, 'RTFS Sensitivity']) * nCases - (1 - metrics[, 'Null Sensitivity']) * nCases
metrics$diffTN <- metrics[, 'RTFS Specificity'] * nControls - metrics[, 'Null Specificity'] * nControls
metrics$diffFP <- (1 - metrics[, 'RTFS Specificity']) * nControls - (1 - metrics[, 'Null Specificity']) * nControls
metrics$diffCorrect <- metrics$diffTP + metrics$diffTN

print(metrics)

write.csv(metrics, '/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/rtfs_w4_filtered_by_family_cox_lasso_confusion_matrix_comparison_0-1_prev33.csv')

