library(dplyr)

config <- yaml::read_yaml(here::here("config.yml"))

readFromResults <- function(filename) {
  readRDS(paste0(config$methylpiper_logs_path, filename))
}

metricsTables <- list(
  epic450k = readFromResults("output_2023_10_01_18_08_51/metricsTable.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/metricsTable.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/metricsTable.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56/metricsTable.rds"),
  pcaTop200kByVariance = readFromResults("output_2023_10_09_11_38_35/metricsTable.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/metricsTable.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/metricsTable.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/metricsTable.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/metricsTable.rds")
)

results <- list(
  epic450k = readFromResults("output_2023_10_01_18_08_51/testResults.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/testResults.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/testResults.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56/testResults.rds"),
  pcaTop200kByVariance = readFromResults("output_2023_10_09_11_38_35/testResults.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/testResults.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/testResults.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/testResults.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/testResults.rds")
)

cvAlphaResults <- list(
  epic450k = readFromResults("output_2023_10_01_18_08_51/cvAlphaResult.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/cvAlphaResult.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/cvAlphaResult.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56/cvAlphaResult.rds"),
  pcaTop200kByVariance = readFromResults("output_2023_10_09_11_38_35/cvAlphaResult.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/cvAlphaResult.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/cvAlphaResult.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/cvAlphaResult.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/cvAlphaResult.rds")
)


modelNames <- c("EPIC-450k",
                "Top 200k by Variance",
                "Top 100k by Variance",
                "PCA EPIC-450k",
                "PCA Top 200k by Variance",
                "Incident T2D EWAS",
                "Prevalent T2D EWAS",
                "Prevalent and Incident T2D EWAS",
                "RTFS",
                "Risk Factors",
                "Risk Factors + PRS",
                "Risk Factors + PRS + Incident T2D EWAS EpiScore")

methylAUCs <- sapply(metricsTables, function(metricsTable) {
  metricsTable[3, 1]
})

methylAUCs <- c(methylAUCs, 
                riskFactors = metricsTables$var200k[1,1], 
                riskFactorsPRS = metricsTables$var200k[2,1],
                riskFactorsPRSEpiScore = metricsTables$ewasIncident[4,1])

methylAUCDF <- data.frame(AUC = methylAUCs)

methylPRAUCs <- sapply(metricsTables, function(metricsTable) {
  metricsTable[3, 2]
})

methylPRAUCs <- c(methylPRAUCs, 
                  riskFactors = metricsTables$var200k[1,2],
                  riskFactorsPRS = metricsTables$var200k[2,2],
                  riskFactorsPRSEpiScore = metricsTables$ewasIncident[4,2])

methylPRAUCDF <- data.frame(PRAUC = methylPRAUCs)

methylMetricsDF <- data.frame(AUC = methylAUCs, PRAUC = methylPRAUCs) %>% 
  mutate(Model = modelNames) %>% 
  arrange(AUC)

bestAlphas <- sapply(row.names(methylPRAUCDF), function(x) {cvAlphaResults[[x]]$bestAlpha})

bestAlphas$riskFactors <- NA
bestAlphas$riskFactorsPRS <- NA
bestAlphas$riskFactorsPRSEpiScore <- bestAlphas$ewasIncident

bestLambdas <- sapply(row.names(methylPRAUCDF), function(x) {cvAlphaResults[[x]]$bestModel$model$lambda.min})

bestLambdas$riskFactors <- NA
bestLambdas$riskFactorsPRS <- NA
bestLambdas$riskFactorsPRSEpiScore <- bestLambdas$ewasIncident

methylMetricsDF <- data.frame(AUC = methylAUCs, PRAUC = methylPRAUCs, Alpha = unlist(bestAlphas), Lambda = unlist(bestLambdas)) %>%
mutate(Model = modelNames) %>%
arrange(AUC)
write.csv(methylMetricsDF, here::here("results", "tables", "aucs_alphas_lambdas.csv"))
