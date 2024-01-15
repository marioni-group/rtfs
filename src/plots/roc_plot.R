# Currently some paths are still absolute as the results files still need to be relocated

library(MethylPipeR)
library(ggplot2)
library(RColorBrewer)
library(MLmetrics)
library(pROC)
library(dplyr)
library(extrafont)
loadfonts()

w1Target <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_59_28/testTarget.rds')
# w1Target$Event <- as.factor(w1Target$Event)

coxTestResultsRTFS <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_59_28/testResults.rds')

coxTestResultsEWAS76 <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_47_09/testResults.rds')

coxTestResultsEPIC450k <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_18_08_51/testResults.rds')

nullResponse <- coxTestResultsRTFS$r$onsetPredictions
rtfsResponse <- coxTestResultsRTFS$d$onsetPredictions
ewas76Response <- coxTestResultsEWAS76$d$onsetPredictions
epic450kResponse <- coxTestResultsEPIC450k$d$onsetPredictions

w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
  if (w1Target$time_to_event[[i]] > 10) {
    0
  } else {
    w1Target$Event[[i]]
  }
})

nullROC <- roc(w1Target$Event, nullResponse)
rtfsROC <- roc(w1Target$Event, rtfsResponse)
ewas76ROC <- roc(w1Target$Event, ewas76Response)
epic450kROC <- roc(w1Target$Event, epic450kResponse)

rocList <- list(nullROC, 
                # epic450kROC, 
                ewas76ROC, 
                rtfsROC)

aucValues <- 
  lapply(rocList, "[[", "auc") %>% 
  as.numeric() %>% 
  round(3)

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


p <- ggroc(rocList) + 
       theme_minimal(base_size = 18, base_family = 'Segoe UI Semilight') + 
       coord_equal() + 
       scale_colour_manual(values=cbPalette, 
                           labels = c(paste0('Risk factors only: ', aucValues[[1]]), 
                                      # paste0('Risk factors + EPIC-450k EpiScore: ', aucValues[[2]]),
                                      paste0('Risk factors + Incident T2D EWAS EpiScore: ', aucValues[[2]]),
                                      paste0('Risk factors + RTFS EpiScore: ', aucValues[[3]]))) +
       guides(color = guide_legend(title = "Model: AUC")) + 
       theme(axis.title.y = element_text(angle = 0, vjust = 0.5), legend.position = c(0.6, 0.15), legend.background = element_rect(fill = "white", color = "black")) +
       ylab('Sensitivity') + xlab('Specificity')
ggsave(here::here("results", "plots", "roc", "roc_plot.png"), dpi = 600, bg = 'white', width = 10, height = 7)

