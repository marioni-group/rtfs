library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsurvfit)

source('/Cluster_Filespace/Marioni_Group/Yipeng/prediction-pipelines/rtfs_20k/gp_hosp_40_60/scripts_20230221/analysis_functions.R')

# load_result <- load450kW3W1()
# target_w3 <- load_result$targetW3
# target_w1 <- load_result$targetW1

target_w3 <- readRDS("/Local_Scratch/Yipeng/w3CoxTable.rds")
target_w1 <- readRDS("/Local_Scratch/Yipeng/w1CoxTable.rds")

# load_result <- NULL
gc()

# target_w1 <- addPRS(target_w1)
# na_prs_index <- is.na(target_w1$SCORE1_SUM)
# target_w1 <- target_w1[!na_prs_index, ]

w3_survfit <- survfit2(Surv(time_to_event, Event) ~ 1, data = target_w3)
w1_survfit <- survfit2(Surv(time_to_event, Event) ~ 1, data = target_w1)

target_w3$set <- "Training"
target_w1$set <- "Test"

kms <- lapply(list(w3 = w3_survfit, w1 = w1_survfit), function(sf) {
  sf %>% ggsurvfit() +
  labs(x = "Time-to-event in years", y = "Survival probability (1 - T2D onset probability)") +
  add_confidence_interval() +
  add_risktable()
})


combined_target <- rbind(target_w3, target_w1)
combined_target$set <- as.factor(combined_target$set)
combined_survfit <- survfit2(Surv(time_to_event, Event) ~ set, data = combined_target)

combined_km <- combined_survfit %>% ggsurvfit() + labs(x = "Time-to-event in years", y = "Survival probability\n(1 - T2D onset probability)") + add_confidence_interval() + add_risktable(
  size = 2, # increase font size of risk table statistics
    theme =   # increase font size of risk table title and y-axis label
      list(
        theme_risktable_default(axis.text.y.size = 8,
                                plot.title.size = 8),
        theme(plot.title = element_text(face = "bold"))
      )
) +
coord_cartesian(xlim = c(0, 15)) +
  scale_y_continuous(
    limits = c(0.85, 1),
    labels = scales::percent,
    expand = c(0.01, 0)
  ) +
  scale_x_continuous(breaks = 0:15, expand = c(0.02, 0)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none") +
  guides(color = guide_legend(ncol = 1))

ggsave("train_test_km.pdf", combined_km, width = 3840, height = 2160, units = "px", dpi = 600)


style_plot <- function(p) {
  p + geom_density(alpha = 0.5) + 
      ylab("") + 
      scale_color_manual(name = "", values = c("#E69F00", "#56B4E9")) +
      # scale_colour_manual(name = "", labels = c("Censored", "Case"), values = c("#E69F00", "#56B4E9")) + 
      scale_linetype_discrete(name = "", labels = c("Censored", "Case")) +
      theme_minimal(base_size = 16) + 
      theme(axis.text.y = element_blank())
}

# w3_tte_densities <- target_w3 %>% ggplot(aes(x = time_to_event, color = as.factor(Event))) + xlab("Time-to-event/censoring") + ggtitle("Set 2 (training set) time-to-event/censoring distribution")
# w3_tte_densities <- w3_tte_densities %>% style_plot

# ggsave("w3_tte_densities.pdf", w3_tte_densities)

# w1_tte_densities <- target_w1 %>% ggplot(aes(x = time_to_event, color = as.factor(Event))) + xlab("Time-to-event/censoring") + ggtitle("Set 1 (test set) time-to-event/censoring distribution")
# w1_tte_densities <- w1_tte_densities %>% style_plot

# ggsave("w1_tte_densities.pdf", w1_tte_densities)


# w3_bmi_densities <- target_w3 %>% ggplot(aes(x = bmi, color = as.factor(Event))) + xlab("BMI") + ggtitle("Set 2 (training set) BMI distribution")
# w3_bmi_densities <- w3_bmi_densities %>% style_plot
# ggsave("w3_bmi_densities.pdf", w3_bmi_densities)


# w1_bmi_densities <- target_w1 %>% ggplot(aes(x = bmi, color = as.factor(Event))) + xlab("BMI") + ggtitle("Set 1 (test set) BMI distribution")
# w1_bmi_densities <- w1_bmi_densities %>% style_plot
# ggsave("w1_bmi_densities.pdf", w1_bmi_densities)


# w3_age_densities <- target_w3 %>% ggplot(aes(x = Age + time_to_event, color = as.factor(Event))) + xlab("Age at event/censoring") + ggtitle("Set 2 (training set) age at event/censoring distribution")
# w3_age_densities <- w3_age_densities %>% style_plot
# ggsave("w3_age_densities.pdf", w3_age_densities)


# w1_age_densities <- target_w1 %>% ggplot(aes(x = Age + time_to_event, color = as.factor(Event))) + xlab("Age at event/censoring") + ggtitle("Set 1 (test set) age at event/censoring distribution")
# w1_age_densities <- w1_age_densities %>% style_plot
# ggsave("w1_age_densities.pdf", w1_age_densities)

combined_bmi_densities <- combined_target %>% ggplot(aes(x = bmi, color = set, linetype = as.factor(Event))) + xlab("BMI") + ggtitle("")
combined_bmi_densities <- combined_bmi_densities %>% style_plot

combined_age_densities <- combined_target %>% ggplot(aes(x = Age + time_to_event, color = set, linetype = as.factor(Event))) + xlab("Age at event/censoring") + ggtitle("")
combined_age_densities <- combined_age_densities %>% style_plot

# combined <- (w3_tte_densities + w1_tte_densities) / (w3_bmi_densities + w1_bmi_densities) / (w3_age_densities + w1_age_densities)
combined_plot <- ggarrange(combined_bmi_densities,
                           combined_age_densities,
                           ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
#                       w3_bmi_densities, w1_bmi_densities,
#                       w3_age_densities, w1_age_densities,
#                       ncol=2, nrow=3, common.legend = TRUE, legend="bottom")

combined_plot <- ggarrange(combined_km, combined_plot, ncol = 1, nrow = 2)

ggsave("summary_densities.pdf", combined_plot, width = 210, height = 297, units = "mm")
