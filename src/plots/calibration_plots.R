library(MethylPipeR)
library(CalibrationCurves)
library(patchwork)

config <- yaml::read_yaml(here::here("config.yml"))

w1Target <- readRDS(paste0(config$methylpiper_logs_path, 'output_2023_10_01_19_59_28/testTarget.rds'))

coxTestResultsRTFS <- readRDS(paste0(config$methylpiper_logs_path, 'output_2023_10_01_19_59_28/testResults.rds'))

coxTestResultsEWAS76 <- readRDS(paste0(config$methylpiper_logs_path, 'output_2023_10_01_19_47_09/testResults.rds'))

nullResponse <- coxTestResultsRTFS$r$onsetPredictions
rtfsResponse <- coxTestResultsRTFS$d$onsetPredictions
ewas76Response <- coxTestResultsEWAS76$d$onsetPredictions

w1Target$Event <- sapply(1:nrow(w1Target), function(i) {
  if (w1Target$time_to_event[[i]] > 10) {
    0
  } else {
    w1Target$Event[[i]]
  }
})

set.seed(42)
ewasCurve <- CalibrationCurves::valProbggplot(ewas76Response, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlab = "")
set.seed(42)
rtfsCurve <- CalibrationCurves::valProbggplot(rtfsResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlab = '')
set.seed(42)
nullCurve <- CalibrationCurves::valProbggplot(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE)

formatPlot <- function(p) {
  p + coord_fixed() + theme_minimal(base_size = 12)
}

ewasPlot <- formatPlot(ewasCurve$ggPlot)
rtfsPlot <- formatPlot(rtfsCurve$ggPlot)
nullPlot <- formatPlot(nullCurve$ggPlot)

pdf(file = here::here("results", "plots", "calibration", "calibration_curves_ggplot.pdf"))
print(ewasPlot / rtfsPlot / nullPlot + plot_layout(guides = "collect"))

dev.off()

