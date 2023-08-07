library(tidyverse)
library(glmnet)
library(MethylPipeR)

traits <- c('age', 'glucose', 'cholest', 'HDL', 'sodium', 'potassium', 'urea', 'creatinine', 'bmi', 'whr', 'fat', 'sBP', 'dBP', 'HR', 'FEV', 'FVC', 'Alc', 'PckYrs', 'g')

# trait_results <- lapply(traits, function(trait) {
#   
# })

# trait_summaries <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/results/methylpiper_logs/output_2022_08_21_13_00_49/testSummaries.rds')

trait_results <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/results/continuous_traits_w4_filtered_by_family/output_2022_10_10_11_44_58/traitResults.rds')
trait_summaries <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/results/continuous_traits_w4_filtered_by_family/output_2023_05_11_15_46_33/testSummaries.rds')

trait_coefficients <- lapply(trait_summaries, function(trait_summary) {trait_summary$coefficients})
trait_adj_r2_values <- sapply(trait_summaries, function(trait_summary) {trait_summary$adj.r.squared})
trait_summaries <- NULL
gc()

trait_episcore_coefficients <- sapply(trait_coefficients, function(df) {df[2,1]})
trait_episcore_std_errs <- sapply(trait_coefficients, function(df) {df[2,2]})
trait_episcore_p_values <- sapply(trait_coefficients, function(df) {df[2,4]})

trait_df <- data.frame(adjusted_r2 = trait_adj_r2_values,
                       episcore_coefficient = trait_episcore_coefficients,
                       episcore_std_err = trait_episcore_std_errs,
                       episcore_p_value = trait_episcore_p_values)
trait_df_formatted <- trait_df
# Format trait_df columns
trait_df_formatted$adjusted_r2 <- sprintf("%.2f", trait_df$adjusted_r2)
trait_df_formatted$episcore_coefficient <- sprintf("%.2f", trait_df$episcore_coefficient)
trait_df_formatted$episcore_std_err <- sprintf("%.3f", trait_df$episcore_std_err)
trait_df_formatted$episcore_p_value <- sprintf("%.5g", trait_df$episcore_p_value)

trait_df_formatted$Trait <- c('Age (Years)', 'Glucose (mmol/L)', 'Total Cholestrol (mmol/L)', 'HDL Cholestrol (mmol/L)', 'Sodium (mmol/L)',
                'Potassium (mmol/L)', 'Urea (mmol/L)', 'Creatinine (umol/L)', 'BMI (kg/m2)', 'Waist-Hip Ratio',
                'Body Fat (%)', 'Systolic Blood Pressure (mmHg)', 'Diastolic Blood Pressure (mmHg)', 'Heart Rate (BPM)', 'Forced Expiratory Volume (L)',
                'Forced Vital Capacity (L)', 'Alcohol Consumption (Units/Week)', 'Smoking (Pack Years)', 'g')

trait_df_formatted <- trait_df_formatted %>% select(Trait, adjusted_r2, episcore_coefficient, episcore_std_err, episcore_p_value)
write.csv(trait_df_formatted, 'continuous_trait_w4_filtered_by_family_test_results.csv', row.names = FALSE)
