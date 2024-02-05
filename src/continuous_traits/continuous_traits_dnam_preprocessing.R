config <- yaml::read_yaml(here::here("config.yml"))

gs20kTargets <- readRDS(paste0(config$continuous_traits_data_path, 'GS20k_Targets.rds'))

phenotypes1 <- readRDS(paste0(config$continuous_traits_data_path, 'GS_pheno_resids_20k_w1.rds'))
phenotypes3 <- readRDS(paste0(config$continuous_traits_data_path, 'GS_pheno_resids_20k_w3.rds'))
phenotypes4 <- readRDS(paste0(config$continuous_traits_data_path, 'GS_pheno_resids_20k_w4.rds'))

phenotypes <- do.call(rbind, list(phenotypes1, phenotypes3, phenotypes4))

phenotypes <- merge(phenotypes, gs20kTargets, by.x = 'ID', by.y = 'Sample_Name')
phenotypes$age <- phenotypes$age.x
phenotypes$age.x <- NULL
phenotypes$age.y <- NULL

methylW1 <- readRDS(paste0(config$raw_dnam_data_path, 'wave1_no_compression.RDS'))
methylW1PhenotypesIndex <- match(colnames(methylW1), phenotypes$Sample_Sentrix_ID)

phenotypesW1 <- phenotypes[methylW1PhenotypesIndex, ]

methylW3 <- readRDS(paste0(config$raw_dnam_data_path, 'wave3_no_compression.RDS'))
methylW3PhenotypesIndex <- match(colnames(methylW3), phenotypes$Sample_Sentrix_ID)

phenotypesW3 <- phenotypes[methylW3PhenotypesIndex, ]

methylW4 <- readRDS(paste0(config$raw_dnam_data_path, 'wave4_no_compression.RDS'))
methylW4PhenotypesIndex <- match(colnames(methylW4), phenotypes$Sample_Sentrix_ID)

phenotypesW4 <- phenotypes[methylW4PhenotypesIndex, ]

# Remove NA rows
naRowsW1 <- is.na(phenotypesW1$Sample_Sentrix_ID)
naRowsW3 <- is.na(phenotypesW3$Sample_Sentrix_ID)
naRowsW4 <- is.na(phenotypesW4$Sample_Sentrix_ID)

nNARowsW1 <- sum(naRowsW1)
nNARowsW3 <- sum(naRowsW3)
nNARowsW4 <- sum(naRowsW4)

methylW1 <- methylW1[, !naRowsW1]
phenotypesW1 <- phenotypesW1[!naRowsW1, ]

methylW3 <- methylW3[, !naRowsW3]
phenotypesW3 <- phenotypesW3[!naRowsW3, ]

methylW4 <- methylW4[, !naRowsW4]
phenotypesW4 <- phenotypesW4[!naRowsW4, ]

gc()

# Get intersection of CpG sites between W1, W3 and W4
print('Calculating common CpG sites between W1, W3 and W4.')
commonCpGSites <- intersect(intersect(row.names(methylW1), row.names(methylW3)), row.names(methylW4))

print('Filtering W1 to common CpG sites.')
methylW1 <- methylW1[commonCpGSites, ]

print('Filtering W3 to common CpG sites.')
methylW3 <- methylW3[commonCpGSites, ]

print('Filtering W4 to common CpG sites.')
methylW4 <- methylW4[commonCpGSites, ]

gc()

# Function for converting M-values to beta values
mToBeta <- function(m) {
  2^m/(2^m + 1)
}

# Mean imputation function for methylation data
imputeMean <- function(data) {
  for (i in 1:ncol(data)) {
    data[is.na(data[, i]), i] <- mean(data[, i], na.rm = TRUE)
  }
  data
}

# Transpose so that rows correspond to individuals and columns correspond to CpG sites
print('Transposing W1.')
methylW1 <- t(methylW1)
print('Transposing W3.')
methylW3 <- t(methylW3)
print('Transposing W4.')
methylW4 <- t(methylW4)

gc()

print('Converting W1 M values to beta values.')
methylW1 <- mToBeta(methylW1)
print('Converting W3 M values to beta values.')
methylW3 <- mToBeta(methylW3)
print('Converting W4 M values to beta values.')
methylW4 <- mToBeta(methylW4)

gc()

print('Performing mean imputation in W1.')
methylW1 <- imputeMean(methylW1)
print('Performing mean imputation in W3.')
methylW3 <- imputeMean(methylW3)
print('Performing mean imputation in W4.')
methylW4 <- imputeMean(methylW4)

gc()

print('Saving phenotypesW1.')
saveRDS(phenotypesW1, paste0(config$continuous_traits_data_path, 'phenotypesW1.rds'))
print('Saving phenotypesW3.')
saveRDS(phenotypesW3, paste0(config$continuous_traits_data_path, 'phenotypesW3.rds'))
print('Saving phenotypesW4.')
saveRDS(phenotypesW4, paste0(config$continuous_traits_data_path, 'phenotypesW4.rds'))

print('Saving methylW1.')
saveRDS(methylW1, paste0(config$continuous_traits_data_path, 'methylW1.rds'))
print('Saving methylW3.')
saveRDS(methylW3, paste0(config$continuous_traits_data_path, 'methylW3.rds'))
print('Saving methylW4.')
saveRDS(methylW4, paste0(config$continuous_traits_data_path, 'methylW4.rds'))
