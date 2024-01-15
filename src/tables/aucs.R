library(dplyr)

readFromResults <- function(filename) {
  readRDS(paste0("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/", filename))
}

# results <- list(
#   lassoEPIC450 = readFromResults("output_2023_02_20_11_44_12/testResults.rds"),
#   lasso200k = readFromResults("output_2023_02_17_11_13_09/testResults.rds"),
#   lasso100k = readFromResults("output_2023_02_20_17_23_30/testResults.rds"),
#   lassoRTFS = readFromResults("output_2023_02_17_12_02_33/testResults.rds"),
#   lassoEWAS76 = readFromResults("output_2023_02_20_18_29_01/testResults.rds"),
#   ridgeEWAS76 = readFromResults("output_2023_02_20_18_40_42/testResults.rds"),
#   lassoEWASIncidentPrevalent = readFromResults("output_2023_02_20_19_01_48/testResults.rds"),
#   ridgeEWASIncidentPrevalent = readFromResults("output_2023_02_20_19_19_47/testResults.rds"),
#   ridgeRTFS = readFromResults("output_2023_02_17_17_07_32/testResults.rds")
# )

metricsTables <- list(
#   pcaEWASIncident = readFromResults("output_2023_02_24_15_54_59/metricsTable.rds"),
  epic450k = readFromResults("output_2023_10_01_18_08_51/metricsTable.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/metricsTable.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/metricsTable.rds"),
#   pcaEWASPrevalent = readFromResults("output_2023_02_27_10_51_43/metricsTable.rds"),
#   pcaEWASPrevalentIncident = readFromResults("output_2023_02_27_11_03_59/metricsTable.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56/metricsTable.rds"),
  pcaTop200kByVariance = readFromResults("output_2023_10_09_11_38_35/metricsTable.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/metricsTable.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/metricsTable.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/metricsTable.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/metricsTable.rds")
#   pcaRTFS = readFromResults("output_2023_05_11_16_33_26/metricsTable.rds")
)

results <- list(
#   pcaEWASIncident = readFromResults("output_2023_02_24_15_54_59/testResults.rds"),
  epic450k = readFromResults("output_2023_10_01_18_08_51/testResults.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/testResults.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/testResults.rds"),
#   pcaEWASPrevalent = readFromResults("output_2023_02_27_10_51_43/testResults.rds"),
#   pcaEWASPrevalentIncident = readFromResults("output_2023_02_27_11_03_59/testResults.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56/testResults.rds"),
  pcaTop200kByVariance = readFromResults("output_2023_10_09_11_38_35/testResults.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/testResults.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/testResults.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/testResults.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/testResults.rds")
#   pcaRTFS = readFromResults("output_2023_05_11_16_33_26/testResults.rds")
)

cvAlphaResults <- list(
#   pcaEWASIncident = readFromResults("output_2023_02_24_15_54_59/cvAlphaResult.rds"),
  epic450k = readFromResults("output_2023_10_01_18_08_51/cvAlphaResult.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/cvAlphaResult.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/cvAlphaResult.rds"),
#   pcaEWASPrevalent = readFromResults("output_2023_02_27_10_51_43/cvAlphaResult.rds"),
#   pcaEWASPrevalentIncident = readFromResults("output_2023_02_27_11_03_59/cvAlphaResult.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56/cvAlphaResult.rds"),
  pcaTop200kByVariance = readFromResults("output_2023_10_09_11_38_35/cvAlphaResult.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/cvAlphaResult.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/cvAlphaResult.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/cvAlphaResult.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/cvAlphaResult.rds")
#   pcaRTFS = readFromResults("output_2023_05_11_16_33_26/cvAlphaResult.rds")
)


modelNames <- c(# "PCA Incident T2D EWAS",
                "EPIC-450k",
                "Top 200k by Variance",
                "Top 100k by Variance",
#                 "PCA Prevalent T2D EWAS",
#                 "PCA Prevalent and Incident T2D EWAS",
                "PCA EPIC-450k",
                "PCA Top 200k by Variance",
                "Incident T2D EWAS",
                "Prevalent T2D EWAS",
                "Prevalent and Incident T2D EWAS",
                "RTFS",
#                 "PCA RTFS",
                "Risk Factors",
                "Risk Factors + PRS",
                "Risk Factors + PRS + Incident T2D EWAS EpiScore")

methylAUCs <- sapply(metricsTables, function(metricsTable) {
  metricsTable[3, 1]
}) # %>% as.data.frame

methylAUCs <- c(methylAUCs, 
                riskFactors = metricsTables$var200k[1,1], 
                riskFactorsPRS = metricsTables$var200k[2,1],
                riskFactorsPRSEpiScore = metricsTables$ewasIncident[4,1])

methylAUCDF <- data.frame(AUC = methylAUCs)

methylPRAUCs <- sapply(metricsTables, function(metricsTable) {
  metricsTable[3, 2]
}) # %>% as.data.frame

methylPRAUCs <- c(methylPRAUCs, 
                  riskFactors = metricsTables$var200k[1,2],
                  riskFactorsPRS = metricsTables$var200k[2,2],
                  riskFactorsPRSEpiScore = metricsTables$ewasIncident[4,2])

methylPRAUCDF <- data.frame(PRAUC = methylPRAUCs)

methylMetricsDF <- data.frame(AUC = methylAUCs, PRAUC = methylPRAUCs) %>% 
  mutate(Model = modelNames) %>% 
  arrange(AUC)

# prsBaselineAUCs <- sapply(metricsTables, function(metricsTable) {
#   metricsTable[4, 1]
# }) # %>% as.data.frame
# 
# prsBaselineAUCs <- c(prsBaselineAUCs, riskFactors = metricsTables$var200k[2,1])
# 
# prsBaselineAUCDF <- data.frame(AUC = prsBaselineAUCs)
# 
# prsBaselinePRAUCs <- sapply(metricsTables, function(metricsTable) {
#   metricsTable[4, 2]
# }) # %>% as.data.frame
# 
# prsBaselinePRAUCs <- c(prsBaselinePRAUCs, riskFactors = metricsTables$var200k[2,2])
# 
# prsBaselinePRAUCDF <- data.frame(PRAUC = prsBaselinePRAUCs)
# 
# prsBaselineBestAlphas <- sapply(row.names(prsBaselinePRAUCDF), function(x) {cvAlphaResults[[x]]$bestAlpha})
# 
# prsBaselineBestAlphas$riskFactors <- NA
# 
# prsBaselineBestLambdas <- sapply(row.names(prsBaselinePRAUCDF), function(x) {cvAlphaResults[[x]]$bestModel$model$lambda.min})
# 
# prsBaselineBestLambdas$riskFactors <- NA

bestAlphas <- sapply(row.names(methylPRAUCDF), function(x) {cvAlphaResults[[x]]$bestAlpha})

bestAlphas$riskFactors <- NA
bestAlphas$riskFactorsPRS <- NA
bestAlphas$riskFactorsPRSEpiScore <- bestAlphas$ewasIncident

bestLambdas <- sapply(row.names(methylPRAUCDF), function(x) {cvAlphaResults[[x]]$bestModel$model$lambda.min})

bestLambdas$riskFactors <- NA
bestLambdas$riskFactorsPRS <- NA
bestLambdas$riskFactorsPRSEpiScore <- bestLambdas$ewasIncident

# prsBaselineMetricsDF <- data.frame(AUC = prsBaselineAUCs, PRAUC = prsBaselinePRAUCs, alpha = unlist(prsBaselineBestAlphas), lambda = unlist(prsBaselineBestLambdas)) %>% 
# mutate(Model = modelNames) %>%
# arrange(AUC)
# write.csv(prsBaselineMetricsDF, 'prsBaselineMetricsDF_gp_smr.csv')
methylMetricsDF <- data.frame(AUC = methylAUCs, PRAUC = methylPRAUCs, Alpha = unlist(bestAlphas), Lambda = unlist(bestLambdas)) %>%
mutate(Model = modelNames) %>%
arrange(AUC)
# select(Model, AUC, PRAUC, Alpha, Lambda)
write.csv(methylMetricsDF, here::here("results", "tables", "aucs_alphas_lambdas.csv"))
