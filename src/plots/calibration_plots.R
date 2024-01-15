library(MethylPipeR)
# library(rms)
library(CalibrationCurves)


source(here::here("src", "plots", "calibration", "val.prob.ci.2.R"))
source(here::here("src", "plots", "calibration", "auc.nonpara.mw.R"))
source(here::here("src", "plots", "calibration", "ci.auc.R"))
source(here::here("src", "plots", "calibration", "BT.samples.R"))

w1Target <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_59_28/testTarget.rds')
# w1Target$Event <- as.factor(w1Target$Event)

coxTestResultsRTFS <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_59_28/testResults.rds')

coxTestResultsEWAS76 <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/results_gp_smr/methylpiper_logs/output_2023_10_01_19_47_09/testResults.rds')

# w1Target <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_diabetes_pipeline/using_methylpiper/src/20k/final_episcores/results/methylpiper_logs/output_2022_03_17_11_17_24/targetW1.rds')

# coxTestResults <- readRDS('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/incident_diabetes_pipeline/using_methylpiper/src/20k/final_episcores/results/methylpiper_logs/output_2022_03_17_11_17_24/testResults.rds')

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

pdf(file = here::here("results", "plots", "calibration", "calibration_curve_rtfs.pdf"))
set.seed(42)
rtfsModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(rtfsResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlab = '')
dev.off()

pdf(file = here::here("results", "plots", "calibration", "calibration_curve_rtfs_zoomed.pdf"))
set.seed(42)
rtfsModelCalibrationCurveZoomed <- val.prob.ci.2(rtfsResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlim = c(-0.02, 0.2), ylim = c(-0.15, 0.4), xlab = '', ylab = '')
dev.off()

pdf(file = here::here("results", "plots", "calibration", "calibration_curve_risk_factors.pdf"))
set.seed(42)
nullModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE)
dev.off()

pdf(file = here::here("results", "plots", "calibration", "calibration_curve_risk_factors_zoomed.pdf"))
set.seed(42)
nullModelCalibrationCurveZoomed <- val.prob.ci.2(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlim = c(-0.02, 0.2), ylim = c(-0.15, 0.4), ylab = '')
dev.off()

pdf(file = here::here("results", "plots", "calibration", "calibration_curve_ewas76.pdf"))
set.seed(42)
ewas76ModelCalibrationCurve <- CalibrationCurves::val.prob.ci.2(ewas76Response, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE)
dev.off()

pdf(file = here::here("results", "plots", "calibration", "calibration_curve_ewas76_zoomed.pdf"))
set.seed(42)
ewas76ModelCalibrationCurveZoomed <- val.prob.ci.2(ewas76Response, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlim = c(-0.02, 0.2), ylim = c(-0.15, 0.4), ylab = '')
dev.off()

pdf(file = here::here("results", "plots", "calibration", "calibration_curves_combined.pdf"))
par(mfrow=c(3,1))
set.seed(42)
CalibrationCurves::val.prob.ci.2(ewas76Response, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlab = "")
set.seed(42)
CalibrationCurves::val.prob.ci.2(rtfsResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE, xlab = '')
set.seed(42)
CalibrationCurves::val.prob.ci.2(nullResponse, w1Target$Event, statloc = FALSE, CL.BT = TRUE, d0lab = '', d1lab = '', legendloc = FALSE)
dev.off()
