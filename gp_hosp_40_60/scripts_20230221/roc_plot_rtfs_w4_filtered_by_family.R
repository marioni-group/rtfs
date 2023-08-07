library(MethylPipeR)
library(ggplot2)
library(RColorBrewer)
library(MLmetrics)
library(pROC)
library(dplyr)
library(extrafont)
loadfonts()

w1Target <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_27_11_18_04/testTarget.rds')
# w1Target$Event <- as.factor(w1Target$Event)

coxTestResultsRTFS <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_05_11_00_44_56/testResults.rds')

coxTestResultsEWAS76 <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_27_10_26_05/testResults.rds')

coxTestResultsEPIC450k <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_26_00_39_01/testResults.rds')

nullResponse <- coxTestResultsRTFS$rPRS$onsetPredictions
rtfsResponse <- coxTestResultsRTFS$dPRS$onsetPredictions
ewas76Response <- coxTestResultsEWAS76$dPRS$onsetPredictions
epic450kResponse <- coxTestResultsEPIC450k$dPRS$onsetPredictions

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
                # epic450kROC, ewas76ROC, 
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
                                      # paste0('Risk factors + Incident T2D EWAS EpiScore: ', aucValues[[3]]),
                                      paste0('Risk factors + RTFS EpiScore: ', aucValues[[2]]))) +
       guides(color = guide_legend(title = "Model: AUC")) + 
       theme(axis.title.y = element_text(angle = 0, vjust = 0.5), legend.position = c(0.6, 0.15), legend.background = element_rect(fill = "white", color = "black")) +
       ylab('Sensitivity') + xlab('Specificity')
ggsave('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/rtfs_null_roc_plot_w4_filtered_by_family.png', dpi = 600, bg = 'white', width = 10, height = 7)

