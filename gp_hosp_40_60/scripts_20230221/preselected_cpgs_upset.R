library(UpSetR)
library(MethylPipeR)
library(glmnet)
library(dplyr)

rtfs_model <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_05_19_15_25_26/methylModel.rds")$model

prevalent_ewas_model <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_27_10_30_42/methylModel.rds")$model

incident_ewas_model <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_27_10_26_05/methylModel.rds")$model

var200k_model <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_24_16_39_40/methylModel.rds")$model

var100k_model <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_26_00_36_53/methylModel.rds")$model

epic450k_model <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_02_26_00_39_01/methylModel.rds")$model

rtfs_cpgs <- coef(rtfs_model, s = "lambda.min")[,1] %>% names
prevalent_ewas_cpgs <- coef(prevalent_ewas_model, s = "lambda.min")[,1] %>% names
incident_ewas_cpgs <- coef(incident_ewas_model, s = "lambda.min")[,1] %>% names
var200k_cpgs <- coef(var200k_model, s = "lambda.min")[,1] %>% names
var100k_cpgs <- coef(var100k_model, s = "lambda.min")[,1] %>% names
epic450k_cpgs <- coef(epic450k_model, s = "lambda.min")[,1] %>% names

cpgs_list <- list(rtfs = rtfs_cpgs,
                  prevalent_ewas = prevalent_ewas_cpgs,
                  incident_ewas = incident_ewas_cpgs,
                  var200k = var200k_cpgs,
                  var100k = var100k_cpgs,
                  epic450k = epic450k_cpgs)

png('preselected_cpgs_upset.png', width = 1920, height = 1080, units = "px")
print(upset(fromList(cpgs_list), sets = names(cpgs_list), order.by = 'freq', mb.ratio = c(0.5, 0.5), text.scale = 3))
dev.off()

png('preselected_cpgs_rtfs_ewas_upset.png', width = 1920, height = 1080, units = "px")
print(upset(fromList(cpgs_list), sets = c("rtfs", "prevalent_ewas", "incident_ewas"), order.by = 'freq', mb.ratio = c(0.5, 0.5), text.scale = 3))
dev.off()
