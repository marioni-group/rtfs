library(dplyr)
library(SurvMetrics)
library(survival)

config <- yaml::read_yaml(here::here("config.yml"))

predictCoxPHOnset <- function(dataDF, coxPHModel, threshold = 10) {
  uniqueTimes <- sort(unique(c(dataDF$time_to_event, threshold)))
  thresholdIndex <- match(threshold, uniqueTimes)
  cumulativeBaseHaz <- gbm::basehaz.gbm(dataDF$time_to_event, dataDF$Event, predict(coxPHModel), uniqueTimes)
  survivalPredictions <- exp(-cumulativeBaseHaz[[thresholdIndex]]) ^ exp(predict(coxPHModel))
  onsetPredictions <- 1 - survivalPredictions

  # Event should be 0 if tte is > threshold
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

time_points <- seq(1, 10, 1)

readFromResults <- function(filename) {
  readRDS(paste0(config$methylpiper_logs_path, filename))
}

results <- list(
  # pcaEWASIncident = readFromResults("output_2023_02_24_15_54_59/testResults.rds"),
  epic450k = readFromResults("output_2023_10_01_18_08_51/testResults.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/testResults.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/testResults.rds"),
  pcaVar200k = readFromResults("output_2023_10_09_11_38_35/testResults.rds"),
  # pcaEWASPrevalent = readFromResults("output_2023_02_27_10_51_43/testResults.rds"),
  # pcaEWASPrevalentIncident = readFromResults("output_2023_02_27_11_03_59/testResults.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56/testResults.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/testResults.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/testResults.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/testResults.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/testResults.rds")
  # pcaRTFS = readFromResults("output_2023_05_11_16_33_26/testResults.rds")
)

models <- list(
  # pcaEWASIncident = readFromResults("output_2023_02_24_15_54_59/models.rds"),
  epic450k = readFromResults("output_2023_10_01_18_08_51/models.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/models.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/models.rds"),
  pcaVar200k = readFromResults("output_2023_10_09_11_38_35/models.rds"),
  # pcaEWASPrevalent = readFromResults("output_2023_02_27_10_51_43/models.rds"),
  # pcaEWASPrevalentIncident = readFromResults("output_2023_02_27_11_03_59/models.rds"),
  pcaEPIC450k = readFromResults("output_2023_10_02_22_37_56//models.rds"),
  ewasIncident = readFromResults("output_2023_10_01_19_47_09/models.rds"),
  ewasPrevalent = readFromResults("output_2023_10_02_08_06_31/models.rds"),
  ewasPrevalentIncident = readFromResults("output_2023_10_02_08_36_58/models.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/models.rds")
  # pcaRTFS = readFromResults("output_2023_05_11_16_33_26/models.rds")
)


modelNames <- c(# "PCA Incident T2D EWAS",
                "EPIC-450k",
                "Top 200k by Variance",
                "Top 100k by Variance",
                "PCA Top 200k by Variance",
                # "PCA Prevalent T2D EWAS",
                # "PCA Prevalent and Incident T2D EWAS",
                "PCA EPIC-450k",
                "Incident T2D EWAS",
                "Prevalent T2D EWAS",
                "Prevalent and Incident T2D EWAS",
                "RTFS",
                # "PCA RTFS",
                "Risk Factors",
                "Risk Factors + PRS",
                "Risk Factors + PRS + Incident T2D EWAS EpiScore")

null_model <- models$epic450k$r
prs_model <- models$epic450k$rPRS
prs_top_episcore_model <- models$ewasIncident$dPRS
full_models <- lapply(models, function(x) {x$d})

cox_models <- c(full_models, "riskFactors" = list(null_model), "riskFactorsPRS" = list(prs_model), "riskFactorsPRSEpiScore" = list(prs_top_episcore_model))

target <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_18_08_51/testTarget.rds")

predictions <- lapply(time_points, function(time_point) {
  model_predictions <- lapply(cox_models, function(model) {
    predictCoxPHOnset(target, model, time_point)$onsetPredictions
  })
  names(model_predictions) <- modelNames
  model_predictions
})

brier_scores <- lapply(time_points, function(time_point) {
  sapply(predictions[[time_point]], function(x) {
    Brier(Surv(target$time_to_event, target$Event), 1 - x, t_star = time_point)
  })
})

c_indexes <- lapply(time_points, function(time_point) {
  sapply(predictions[[time_point]], function(x) {
    Cindex(Surv(target$time_to_event, target$Event), 1 - x, t_star = time_point)
  })
})

brier_scores_df <- do.call(cbind, brier_scores)
c_indexes <- c_indexes[[1]]
names(c_indexes) <- modelNames
c_indexes_df <- data.frame(c_indexes)

row.names(brier_scores_df) <- modelNames

write.csv(brier_scores_df, here::here("results", "tables", "brier_scores.csv"))
write.csv(c_indexes_df, here::here("results", "tables", "c_indexes.csv"))
