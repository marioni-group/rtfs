library(dplyr)
library(purrr)

source('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/script_pheno_prep_20k_separate_waves.R')

originals <- lapply(results, function(wave_result) {
  wave_result$original
})

# Reverse transform traits
originals <- lapply(originals, function(original) {
  original$bmi <- exp(original$bmi)
  original$Alc <- exp(original$Alc) - 1
  original$PckYrs <- exp(original$PckYrs) - 1
  original
})

# summaries <- lapply(originals, function(original) {
#   original %>% summary
# })

# For each wave, calculate statistics over each trait.
originals_stats <- lapply(originals, function(original) {
  mins <- apply(original, 2, function(x) {min(x, na.rm = TRUE)})
  maxes <- apply(original, 2, function(x) {max(x, na.rm = TRUE)})
  means <- apply(original, 2, function(x) {mean(x, na.rm = TRUE)})
  sds <- apply(original, 2, function(x) {sd(x, na.rm = TRUE)})
  
  total_rows <- nrow(original)
  count_nas <- apply(original, 2, function(x){sum(is.na(x))})
  ns <- total_rows - count_nas
  
  df <- data.frame(N = ns, Min = sprintf('%.1f', mins), Mean = sprintf('%.1f', means), Max = sprintf('%.1f', maxes), SD = sprintf('%.1f', sds), NAs = count_nas)
  # Remove ID column
  df <- df[-1,]
  df$Trait <- c('Age (Years)', 'Glucose (mmol/L)', 'Total Cholestrol (mmol/L)', 'HDL Cholestrol (mmol/L)', 'Sodium (mmol/L)',
                'Potassium (mmol/L)', 'Urea (mmol/L)', 'Creatinine (umol/L)', 'BMI (kg/m2)', 'Waist-Hip Ratio',
                'Body Fat (%)', 'Systolic Blood Pressure (mmHg)', 'Diastolic Blood Pressure (mmHg)', 'Heart Rate (BPM)', 'Forced Expiratory Volume (L)',
                'Forced Vital Capacity (L)', 'Alcohol Consumption (Units/Week)', 'Smoking (Pack Years)', 'g')
  df <- df[, c('Trait', 'N', 'Min', 'Mean', 'Max', 'SD', 'NAs')]
  df
})

# w4_nas <- readRDS("w4_family_filtered_nas.rds")
# w1_nas <- readRDS("w1_family_filtered_nas.rds")

# originals_stats$w4$NAs <- w4_nas
# originals_stats$w1$NAs <- w1_nas
# 
# sample_sizes <- lapply(originals, function(original) {
#   total_rows <- nrow(original)
#   count_nas <- apply(original, 2, function(x){sum(is.na(x))})
#   total_rows - count_nas
# })

# write stats to file
map2(originals_stats, names(originals_stats), function(x, y) {write.csv(x, here::here("results", "tables", paste0('continuous_stats_', y, '.csv')), row.names = FALSE)})
