# source('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_IHD_pipeline/src/metrics.R')
library(MethylPipeR)
library(caret)
library(reshape2)
library(patchwork)

# Load Wave 1 target
w1Target <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_59_28/testTarget.rds')

# Event should be 0 if tte is > 10
w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
  if (w1Target$time_to_event[[i]] > 10) {
    0
  } else {
    w1Target$Event[[i]]
  }
})

rtfsTestResults <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_59_28/testResults.rds')
ewas76TestResults <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_47_09/testResults.rds')

nullResponse <- rtfsTestResults$rPRS$onsetPredictions
rtfsResponse <- rtfsTestResults$dPRS$onsetPredictions
ewas76Response <- ewas76TestResults$dPRS$onsetPredictions

# incrementalCoxLassoResponse <- fullResponse

# nullMetrics0.5 <- MLmetrics::ConfusionMatrix(nullResponse > 0.5, w1Target$Event)
# incrementalCoxLassoMetrics0.5 <- MLmetrics::ConfusionMatrix(incrementalCoxLassoResponse > 0.5, w1Target$Event)

# caretNullMetrics0.5 <- caret::confusionMatrix(as.factor(nullResponse > 0.5), as.factor(w1Target$Event == 1), positive = 'TRUE')
# caretIncrementalCoxLassoMetrics0.5 <- caret::confusionMatrix(as.factor(incrementalCoxLassoResponse > 0.5), as.factor(w1Target$Event == 1), positive = 'TRUE')

calculateMetricsAtThreshold <- function(predicted, actual, threshold) {
  caret::confusionMatrix(as.factor(predicted > threshold), as.factor(actual == 1), positive = 'TRUE')
}

# nullMetrics0.5 <- calculateMetricsAtThreshold(nullResponse, w1Target$Event, 0.5)
# incrementalCoxLassoMetrics0.5 <- calculateMetricsAtThreshold(incrementalCoxLassoResponse, w1Target$Event, 0.5)

thresholds <- seq(0.1, 1, 0.1)

incrementalCoxLassoMetrics <- lapply(thresholds, function(threshold) {list(
  nullModel = calculateMetricsAtThreshold(nullResponse, w1Target$Event, threshold),
  rtfsModel = calculateMetricsAtThreshold(rtfsResponse, w1Target$Event, threshold),
  ewas76Model = calculateMetricsAtThreshold(ewas76Response, w1Target$Event, threshold)
)})

names(incrementalCoxLassoMetrics) <- thresholds



nullTruePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[2,2]
})

rtfsTruePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$rtfsModel$table[2,2]
})

ewas76TruePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$ewas76Model$table[2,2]
})

nullTrueNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[1,1]
})

rtfsTrueNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$rtfsModel$table[1,1]
})

ewas76TrueNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$ewas76Model$table[1,1]
})

nullFalsePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[2,1]
})

rtfsFalsePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$rtfsModel$table[2,1]
})

ewas76FalsePositives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$ewas76Model$table[2,1]
})

nullFalseNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$nullModel$table[1,2]
})

rtfsFalseNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$rtfsModel$table[1,2]
})

ewas76FalseNegatives <- sapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$ewas76Model$table[1,2]
})

tables <- lapply(thresholds, function(threshold) {
  incrementalCoxLassoMetrics[[toString(threshold)]]$rtfsModel$table
})

pdf(here::here("results", "plots", "confusion_matrix", "confusion_matrix_plot.pdf"), width = 8, height = 8)
par(mfrow = c(2,2))

plot(thresholds, rtfsTruePositives, type = 'o', col = 'blue', pch='o',
     main = 'True positives', 
     xlab = '', ylab = '',
     ylim = c(0, 150), las = 1)
mtext('N', side = 2, line = 3, las = 1)
points(thresholds, nullTruePositives, col = 'red', pch='*')
lines(thresholds, nullTruePositives, col = 'red')
points(thresholds, ewas76TruePositives, col = 'brown', pch=3)
lines(thresholds, ewas76TruePositives, col = 'brown')
legend('topright', legend = c('Risk factors only', 'RTFS EpiScore', 'Incident T2D EWAS EpiScore'), col = c('red', 'blue', 'brown'), lty = 1)

plot(thresholds, rtfsFalseNegatives, type = 'o', col = 'blue', pch='o',
     main = 'False negatives',
     xlab = '', ylab = '',
     ylim = c(50, 250), las = 1)
points(thresholds, nullFalseNegatives, col = 'red', pch='*')
lines(thresholds, nullFalseNegatives, col = 'red')
points(thresholds, ewas76FalseNegatives, col = 'brown', pch=3)
lines(thresholds, ewas76FalseNegatives, col = 'brown')

plot(thresholds, rtfsFalsePositives, type = 'o', col = 'blue', pch='o',
     main = 'False positives', 
     xlab = 'Threshold', ylab = '',
     ylim = c(0, 500), las = 1)
mtext('N', side = 2, line = 3, las = 1)
points(thresholds, nullFalsePositives, col = 'red', pch='*')
lines(thresholds, nullFalsePositives, col = 'red')
points(thresholds, ewas76FalsePositives, col = 'brown', pch=3)
lines(thresholds, ewas76FalsePositives, col = 'brown')

plot(thresholds, rtfsTrueNegatives, type = 'o', col = 'blue', pch='o',
     main = 'True negatives',
     xlab = 'Threshold', ylab = '',
     ylim = c(4100, 4600), las = 1)
points(thresholds, nullTrueNegatives, col = 'red', pch='*')
lines(thresholds, nullTrueNegatives, col = 'red')
points(thresholds, ewas76TrueNegatives, col = 'brown', pch=3)
lines(thresholds, ewas76TrueNegatives, col = 'brown')
dev.off()


# png('rtfs_metrics_at_different_thresholds_w4_filtered_by_family.png', width = 8, height = 8, units = 'in', res = 300)
# par(mfrow = c(2,2))

# plot(thresholds, rtfsTruePositives, type = 'o', col = 'blue', pch='o',
#      main = 'True positives',
#      xlab = '', ylab = '',
#      ylim = c(0, 150), las = 1)
# mtext('N', side = 2, line = 3, las = 1)
# points(thresholds, nullTruePositives, col = 'red', pch='*')
# lines(thresholds, nullTruePositives, col = 'red')
# legend('topright', legend = c('Risk factors only', 'Risk factors + RTFS EpiScore'), col = c('red', 'blue'), lty = 1)
# 
# plot(thresholds, rtfsFalseNegatives, type = 'o', col = 'blue', pch='o',
#      main = 'False negatives',
#      xlab = '', ylab = '',
#      ylim = c(50, 250), las = 1)
# points(thresholds, nullFalseNegatives, col = 'red', pch='*')
# lines(thresholds, nullFalseNegatives, col = 'red')
# 
# plot(thresholds, rtfsFalsePositives, type = 'o', col = 'blue', pch='o',
#      main = 'False positives',
#      xlab = 'Threshold', ylab = '',
#      ylim = c(0, 500), las = 1)
# mtext('N', side = 2, line = 3, las = 1)
# points(thresholds, nullFalsePositives, col = 'red', pch='*')
# lines(thresholds, nullFalsePositives, col = 'red')
# 
# plot(thresholds, rtfsTrueNegatives, type = 'o', col = 'blue', pch='o',
#      main = 'True negatives',
#      xlab = 'Threshold', ylab = '',
#      ylim = c(4100, 4600), las = 1)
# points(thresholds, nullTrueNegatives, col = 'red', pch='*')
# lines(thresholds, nullTrueNegatives, col = 'red')
# dev.off()
# 
# fullTPRs <- fullTruePositives / (fullTruePositives + fullFalseNegatives)
# nullTPRs <- nullTruePositives / (nullTruePositives + nullFalseNegatives)

# fullTNRs <- fullTrueNegatives / (fullTrueNegatives + fullFalsePositives)
# nullTNRs <- nullTrueNegatives / (nullTrueNegatives + nullFalsePositives)

# fullPPVs <- fullTruePositives / (fullTruePositives + fullFalsePositives)
# nullPPVs <- nullTruePositives / (nullTruePositives + nullFalsePositives)

# fullNPVs <- fullTrueNegatives / (fullTrueNegatives + fullFalseNegatives)
# nullNPVs <- nullTrueNegatives / (nullTrueNegatives + nullFalseNegatives)


