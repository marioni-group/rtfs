library(dplyr)
library(survival)

config <- yaml::read_yaml(here::here("config.yml"))


apr_2022_metrics_table <- readRDS(paste0(config$methylpiper_logs_path, "output_2023_10_01_19_59_28/metricsTable.rds"))

oct_2020_metrics_table <- readRDS(paste0(config$methylpiper_logs_path, "output_2023_12_11_14_05_33/metricsTable.rds"))

apr_2022_models <- readRDS(paste0(config$methylpiper_logs_path, "output_2023_10_01_19_59_28/models.rds"))

oct_2020_models <- readRDS(paste0(config$methylpiper_logs_path, "output_2023_12_11_14_05_33/models.rds"))

apr_2022_summary <- apr_2022_models$d %>% summary
oct_2020_summary <- oct_2020_models$d %>% summary

apr_2022_coef_table <- apr_2022_summary$coefficients
oct_2020_coef_table <- oct_2020_summary$coefficients

write.csv(apr_2022_metrics_table, here::here("results", "tables", "rtfs_apr_2022_metrics_table.csv"))
write.csv(oct_2020_metrics_table, here::here("results", "tables", "rtfs_oct_2020_metrics_table.csv"))

write.csv(apr_2022_coef_table, here::here("results", "tables", "rtfs_apr_2022_coef_table.csv"))
write.csv(oct_2020_coef_table, here::here("results", "tables", "rtfs_oct_2020_coef_table.csv"))

