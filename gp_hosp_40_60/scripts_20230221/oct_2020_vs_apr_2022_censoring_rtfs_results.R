library(dplyr)
library(survival)

apr_2022_metrics_table <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_07_20_18_29_36/metricsTable.rds")

oct_2020_metrics_table <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_07_20_18_45_51/metricsTable.rds")

apr_2022_models <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_07_20_18_29_36/models.rds")

oct_2020_models <- readRDS("/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/output_2023_07_20_18_45_51/models.rds")

apr_2022_summary <- apr_2022_models$dPRS %>% summary
oct_2020_summary <- oct_2020_models$dPRS %>% summary

apr_2022_coef_table <- apr_2022_summary$coefficients
oct_2020_coef_table <- oct_2020_summary$coefficients


