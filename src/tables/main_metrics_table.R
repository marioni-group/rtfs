library(dplyr)

auc_alpha_lambda_table <- read.csv(here::here("results", "tables", "aucs_alphas_lambdas.csv"))
c_index_table <- read.csv(here::here("results", "tables", "c_indexes.csv")) %>% rename(Model = X, `C-index` = c_indexes)
combined_table <- auc_alpha_lambda_table %>% left_join(c_index_table) 

# Change to clearer model labels
new_model_labels <- c("Risk Factors-only",
                      "PCA EPIC-450k EpiScore",
                      "PCA Top 200k by Variance EpiScore",
                      "EPIC-450k EpiScore",
                      "PRS",
                      "Top 100k by Variance EpiScore",
                      "Prevalent T2D EWAS EpiScore",
                      "Top 200k by Variance EpiScore",
                      "Prevalent and Incident T2D EWAS EpiScore",
                      "RTFS EpiScore",
                      "Incident T2D EWAS EpiScore",
                      "PRS + Incident T2D EWAS EpiScore")

X_column <- c("riskFactors",
              "pcaEPIC450k",
              "pcaTop200kByVariance",
              "epic450k",
              "riskFactorsPRS",
              "var100k",
              "ewasPrevalent",
              "var200k",
              "ewasPrevalentIncident",
              "rtfs",
              "ewasIncident",
              "riskFactorsPRSEpiScore")

new_label_df <- data.frame(X = X_column, Incremental_Model_Variables = new_model_labels)

combined_table <- combined_table %>% left_join(new_label_df) %>% select(Incremental_Model_Variables, AUC, PRAUC, `C-index`, Alpha, Lambda)
write.csv(combined_table, here::here("results", "tables", "main_metrics_table.csv"), row.names = FALSE)
