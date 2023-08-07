library(doMC)
library(doParallel)
registerDoMC(cores = 3)
registerDoParallel(3)

library(MethylPipeR)
library(survival)

source('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/analysis_functions.R')

startTimestamp <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')

initLogs('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results/methylpiper_logs/', note = 'Cox lasso protein and direct EpiScore predictors. Direct trained on w3 only')

set.seed(42)

cpgs <- read.csv('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/ewas_76_cpgs.csv')$cpg

loadResult <- load450kW3W1()
targetW3 <- loadResult$targetW3
methylW3 <- loadResult$methylW3
targetW1 <- loadResult$targetW1
methylW1 <- loadResult$methylW1


cpgs <- intersect(colnames(methylW3), cpgs)
methylW3 <- methylW3[, cpgs]
gc()

methylW1 <- methylW1[, cpgs]
gc()

# Scale methylation data
methylW3 <- scale(methylW3)
methylW1 <- scale(methylW1)
gc()

# Add family information
targetW3 <- addFamilyInformation(targetW3)
targetW1 <- addFamilyInformation(targetW1)

# Add T2D PRS
targetW1 <- addPRS(targetW1)
naPRSIndex <- is.na(targetW1$SCORE1_SUM)
targetW1 <- targetW1[!naPRSIndex, ]
methylW1 <- methylW1[!naPRSIndex, ]

writeLines('Loaded data')

set.seed(42)

w1ShuffleIndex <- sample(1:nrow(targetW1))
methylW1 <- methylW1[w1ShuffleIndex, ]
targetW1 <- targetW1[w1ShuffleIndex, ]
row.names(targetW1) <- NULL

set.seed(42)

writeLines('Fitting direct T2D Lasso')

results <- fitAndPredict(methylW3, targetW3, methylW1, targetW1, seed = 42, nFolds = 9, standardize = FALSE, searchAlphas = seq(0, 1, 0.1))

saveResults(results)

gc()

