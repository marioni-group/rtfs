library(MethylPipeR)
library(survival)

config <- yaml::read_yaml(here::here("config.yml")) 

source(here::here("src", "analysis_functions.R"))

startTimestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

initLogs(config$methylpiper_logs_path, note = 'Cox elastic-net EpiScore predictor, CpGs filtered to EPIC-450k intersection. Trained on w3 only')

set.seed(42)

loadResult <- load450kW3W1(censoring = "apr_2022")
targetW3 <- loadResult$targetW3
methylW3 <- loadResult$methylW3
targetW1 <- loadResult$targetW1
methylW1 <- loadResult$methylW1

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

