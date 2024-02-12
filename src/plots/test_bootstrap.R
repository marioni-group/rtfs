library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(extrafont)
loadfonts()

config <- yaml::read_yaml(here::here("config.yml"))

source(here::here("src", "plots", "bootstrap", "test_bootstrap_functions.R"))

readFromResults <- function(filename) {
  readRDS(paste0(config$methylpiper_logs_path, filename))
}

targets <- list(
  epic450k = readFromResults("output_2023_10_01_18_08_51/testTarget.rds"),
  var200k = readFromResults("output_2023_10_01_18_35_10/testTarget.rds"),
  var100k = readFromResults("output_2023_10_01_19_00_59/testTarget.rds"),
  ewas76 = readFromResults("output_2023_10_01_19_47_09/testTarget.rds"),
  ewas58 = readFromResults("output_2023_10_02_08_06_31/testTarget.rds"),
  ewas5876 = readFromResults("output_2023_10_02_08_36_58/testTarget.rds"),
  rtfs = readFromResults("output_2023_10_01_19_59_28/testTarget.rds")
)

nBootstraps <- 1000

# bootstrapResults <- lapply(targets, function(target) {
#   lapply(1:nBootstraps, function(bootstrapSeed) {
#     bootstrapTest(target, bootstrapSeed)
#   })
# })

# saveRDS(bootstrapResults, paste0(config$bootstrap_results_path, "bootstrapResultsGPSMR.rds"))

# Load previously saved results
bootstrapResults <- readRDS(paste0(config$bootstrap_results_path, "bootstrapResultsGPSMR.rds"))

# stop()

bootstrapResults <- bootstrapResults[c("epic450k", "var200k", "var100k", "ewas76", "ewas58", "ewas5876", "rtfs")]

fullAUCsByMethod <- sapply(bootstrapResults, function(result) {
  sapply(result, function(r) {
    r$aucs[['d']]
  })
}) %>% data.frame

fullAUCsByMethod$riskFactors <- sapply(bootstrapResults$var200k, function(r) {
  r$aucs[['r']]
})

fullPRAUCsByMethod <- sapply(bootstrapResults, function(result) {
  sapply(result, function(r) {
    r$praucs[['d']]
  })
}) %>% data.frame

fullPRAUCsByMethod$riskFactors <- sapply(bootstrapResults$var200k, function(r) {
  r$praucs[['r']]
})


aucQuantiles <- sapply(fullAUCsByMethod, function(aucs) {quantile(aucs, c(0.05, 0.5, 0.95))})

aucDiffRiskFactors <- lapply(names(targets), function(m) {
  fullAUCsByMethod[[m]] - fullAUCsByMethod$riskFactors
})
names(aucDiffRiskFactors) <- names(targets)

columnNames <- names(fullAUCsByMethod)
rowNames <- names(fullAUCsByMethod)
matrixNames <- paste0('bootstrap', 1:nBootstraps)

diffArray <- array(1:(nBootstraps*length(fullAUCsByMethod)^2), dim = c(length(fullAUCsByMethod), length(fullAUCsByMethod), nBootstraps), dimnames = list(rowNames, columnNames, matrixNames))

for (i in names(fullAUCsByMethod)) {
  for (j in names(fullAUCsByMethod)) {
    for (k in 1:nBootstraps) {
      diffArray[i, j, paste0('bootstrap', k)] <- fullAUCsByMethod[[i]][[k]] - fullAUCsByMethod[[j]][[k]]
    }
  }
}

meanDiffArray <- apply(diffArray, 3, mean)

models <- row.names(meanDiffArray)

# Heatmap 
meanDiffPlot <- meanDiffArray %>%
  
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Z", -1) %>%
  mutate(X = sapply(X, function(i) {models[[i]]})) %>%
  
  # Change Y to numeric
  # mutate(Y = match(Y, models)) %>%

# Viz
  ggplot(aes(X, Y, fill= -Z)) + 
    geom_tile() +
    theme_minimal() +
    # theme(legend.position="none") +
    scale_fill_gradientn(colours = c("red", "white", "green")) +
    xlab("") + ylab("") +
    guides(fill=guide_legend(title="Mean Bootstrap AUC Difference (Y over X)"))

aucRankOrdersByBootstrap <- apply(-fullAUCsByMethod, 1, function(x){rank(x, ties.method = "min")})
aucRankOrdersByBootstrapFreq <- apply(aucRankOrdersByBootstrap, 1, table)
names(aucRankOrdersByBootstrapFreq) <- names(fullPRAUCsByMethod)

praucRankOrdersByBootstrap <- apply(-fullPRAUCsByMethod, 1, function(x) {rank(x, ties.method = "min")})
praucRankOrdersByBootstrapFreq <- apply(praucRankOrdersByBootstrap, 1, table)
names(praucRankOrdersByBootstrapFreq) <- names(fullPRAUCsByMethod)

# x is the vector of rank order frequencies (by method)
# r is the rank (from 1 to the number of methods)
aucRankFrequencies <- lapply(aucRankOrdersByBootstrapFreq, function(x) {
  sapply(1:length(aucRankOrdersByBootstrapFreq), function(r) {
    if (as.character(r) %in% names(x)) {
      x[[as.character(r)]]
    } else {
      0
    }
  })
})

aucRankFrequencyDFs <- lapply(names(aucRankFrequencies), function(method) {
  frequencies <- aucRankFrequencies[[method]]
  as.data.frame(list(frequency = frequencies, rank = 1:length(frequencies), method = method))
})

aucRankFrequencyDF <- do.call(rbind, aucRankFrequencyDFs)

# x is the vector of rank order frequencies (by method)
# r is the rank (from 1 to the number of methods)
praucRankFrequencies <- lapply(praucRankOrdersByBootstrapFreq, function(x) {
  sapply(1:length(praucRankOrdersByBootstrapFreq), function(r) {
    if (as.character(r) %in% names(x)) {
      x[[as.character(r)]]
    } else {
      0
    }
  })
})

praucRankFrequencyDFs <- lapply(names(praucRankFrequencies), function(method) {
  frequencies <- praucRankFrequencies[[method]]
  as.data.frame(list(frequency = frequencies, rank = 1:length(frequencies), method = method))
})

praucRankFrequencyDF <- do.call(rbind, praucRankFrequencyDFs)


# Colour-blind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

aucRankFrequencyPlot <- aucRankFrequencyDF %>% filter(method %in% c("epic450k", "var200k", "ewas76", "rtfs", "riskFactors")) %>%
                          ggplot(aes(x = rank, y = frequency, colour = as.factor(method))) +
                            geom_point() +
                            geom_line() +
                            theme_minimal() +
                            scale_colour_manual(values = cbPalette)

pdf(here::here("results", "plots", "bootstrap", "auc_rank_frequency_plot.pdf"))
print(aucRankFrequencyPlot)
dev.off()

praucRankFrequencyPlot <- praucRankFrequencyDF %>% filter(method %in% c("epic450k", "var200k", "ewas76", "rtfs", "riskFactors")) %>%
                            ggplot(aes(x = rank, y = frequency, colour = as.factor(method))) +
                              geom_point() +
                              geom_line() +
                              theme_minimal() +
                              scale_colour_manual(values = cbPalette)

pdf(here::here("results", "plots", "bootstrap", "prauc_rank_frequency_plot.pdf"))
print(praucRankFrequencyPlot)
dev.off()

aucRankCumulativeDF <- aucRankFrequencyDF %>% 
  group_by(method) %>% 
  mutate(cumulative = cumsum(frequency)) %>% 
  ungroup

praucRankCumulativeDF <- praucRankFrequencyDF %>%
  group_by(method) %>%
  mutate(cumulative = cumsum(frequency)) %>%
  ungroup

plotRankCumulative <- function(df, plotCols, plotColLegendNames, xLabel) {
  # Colour-blind friendly palette
  cbPalette <- c(# "#999999", 
                 "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  p <- df %>% filter(method %in% plotCols) %>%
         ggplot(aes(x = rank, y = cumulative, colour = as.factor(method))) +
           geom_point(alpha = 0.8) +
           geom_line(alpha = 0.8) +
           theme_minimal(base_size = 18) + # , base_family = 'Segoe UI SemiLight') +
           scale_colour_manual(values = cbPalette,
                               breaks = plotCols,
                               labels = plotColLegendNames) +
           theme(legend.title=element_blank()) +
           ylab("Cumulative frequency in 1000 bootstrap samples") +
           xlab(xLabel)
}

pdf(here::here("results", "plots", "bootstrap", "auc_rank_cumulative_plot.pdf"), width = 10, height = 7) 
aucRankCumulativePlot <- plotRankCumulative(aucRankCumulativeDF, 
                                            plotCols = c("rtfs", "ewas76", "var200k", "epic450k", "riskFactors"), 
                                            plotColLegendNames = c("RTFS", "Incident T2D EWAS", "Top 200k by Variance", "EPIC-450k Intersection", "Risk Factors Only"),
                                            xLabel = "Top n Rank by AUC")
print(aucRankCumulativePlot)
dev.off()

pdf(here::here("results", "plots", "bootstrap", "prauc_rank_cumulative_plot.pdf"), width = 10, height = 7) 
praucRankCumulativePlot <- plotRankCumulative(praucRankCumulativeDF,
                                              plotCols = c("rtfs", "ewas76", "var200k", "epic450k", "riskFactors"),
                                              plotColLegendNames = c("RTFS", "Incident T2D EWAS", "Top 200k by Variance", "EPIC-450k Intersection", "Risk Factors Only"),
                                              xLabel = "Top n Rank by PRAUC")
print(praucRankCumulativePlot)
dev.off()

# pngs
png(here::here("results", "plots", "bootstrap", "auc_rank_cumulative_plot.png"), width = 10, height = 7) 
aucRankCumulativePlot <- plotRankCumulative(aucRankCumulativeDF,
                                            plotCols = c("rtfs", "ewas76", "var200k", "epic450k", "riskFactors"),
                                            plotColLegendNames = c("RTFS", "Incident T2D EWAS", "Top 200k by Variance", "EPIC-450k Intersection", "Risk Factors Only"),
                                            xLabel = "Top n Rank by AUC")
print(aucRankCumulativePlot)
dev.off()

png(here::here("results", "plots", "bootstrap", "prauc_rank_cumulative_plot.png"), width = 10, height = 7)
praucRankCumulativePlot <- plotRankCumulative(praucRankCumulativeDF,
                                              plotCols = c("rtfs", "ewas76", "var200k", "epic450k", "riskFactors"),
                                              plotColLegendNames = c("RTFS", "Incident T2D EWAS", "Top 200k by Variance", "EPIC-450k Intersection", "Risk Factors Only"),
                                              xLabel = "Top n Rank by PRAUC")
print(praucRankCumulativePlot)
dev.off()
