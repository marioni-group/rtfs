library(survival)

# Assumes that the predicted scores are included as a column in target
bootstrapTest <- function(target, bootstrapSeed) {
  set.seed(bootstrapSeed)
  bootstrapIndex <- sample(1:nrow(target), replace = TRUE)
  target <- target[bootstrapIndex, ]
  
  riskFactorsOnlyCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes, target)
  riskFactorsPRSCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + SCORE1_SUM, target)
  directCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + dScore, target)
  directPRSCoxPH <- coxph(Surv(time_to_event, Event) ~ Age + sex + bmi + high_BP + family_diabetes + SCORE1_SUM + dScore, target)

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
  testResults <- lapply(models, function(m) {predictCoxPHOnset(target, m)})

  aucs <- sapply(testResults, function(r) {r$auc})
  praucs <- sapply(testResults, function(r) {r$prauc})
  metricsTable <- data.frame(AUC = aucs, PRAUC = praucs)

  row.names(metricsTable) <- c('Risk factors only',
                               'Risk factors + PRS',
                               'Risk factors + direct EpiScore',
                               'Risk factors + PRS + direct EpiScore')
  
  list(testResults = testResults, aucs = aucs, praucs = praucs, metricsTable = metricsTable, target = target)
}

