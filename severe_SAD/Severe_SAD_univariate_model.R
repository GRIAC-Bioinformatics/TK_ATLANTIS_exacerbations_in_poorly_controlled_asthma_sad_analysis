# This code will use univariate cox regression to predict time to exacerbation in 3 asthma control groups 
# using 3 parameters to define severe SAD

#########################################################
#packages inladen
library(survival)
library(survminer)
library(dplyr)
library(yaml)
library(here)

cfg <- yaml::read_yaml(here("config.yaml"))

#########################################################
### load the clinical table + exacerbation + + severe sad definition table 
db <- read.csv(cfg$paths$patient_data, header = TRUE) %>%
  filter(VISIT == "VISIT 1", ASTHEA == "A", f_eval == "Y")

load(cfg$paths$exacerbations_data) # from Tessa K.
db_exacerbations <- joined %>%
  mutate(PT = as.numeric(PT)) %>%
  dplyr::select(c("PT", "date_baseline", "date_exacerbation", "time_to_exacerbation",
                  "exacerbation", "time")) %>%
  mutate(exacerbation = as.integer(exacerbation)) %>%
  mutate(time = if_else(time>420, 420, time))


## add severe SAD info- columns "R520_higher_3rsd", "X5_lower_3rsd", "AX_higher_3rsd"
db_exacerbations <- db_exacerbations %>%
  left_join(read.csv(cfg$paths$severe_sad_table), by = c("PT"="PT")) %>%
  left_join(db, by = c("PT" = "PT"))
  
#########################################################
## prepare table: 
db_exacerbations <- db_exacerbations %>% 
  mutate(
  ACQ6group = case_when(
    acq6_score < 0.75                  ~ "well controlled asthma",
    acq6_score >= 0.75 & acq6_score <= 1.5 ~ "partially controlled asthma",
    acq6_score > 1.5                   ~ "poorly controlled asthma",
    TRUE                                ~ NA_character_
  ),
  ACQ6group = factor(ACQ6group)  # convert to factor
)

#########################################################
# Function: 
# Survival analyse: severe sad - time to exacerbation - in 3 groups of asthma control

sad_survival_uni <- function(sad_variable, db_exacerbations, subgroup = NULL) {
  
  # no NAs for the sad variable of interest:
  filter_df <- db_exacerbations %>%
    filter(!is.na(.data[[sad_variable]])) %>%
    mutate(!!sym(sad_variable) := as.factor(!!sym(sad_variable)))
  # asthma control group - if defined 
  if (!is.null(subgroup)) {
    filter_df <- filter_df %>%
      filter(.data[["ACQ6group"]] == subgroup)
    
    message(paste("Filtering for", subgroup))
  }
  
  # Define the Formula
  formula_obj <- as.formula(paste0("Surv(time, exacerbation) ~ ", sad_variable))
  # KM fit
  km <- survfit(formula_obj, data = filter_df)
  
  # Log-rank test
  logrank_test <- survdiff(formula_obj, data = filter_df)
  p_value_logrank <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  
  # HR + P-value from the Cox-model
  cox_model <- coxph(formula_obj, data = filter_df)
  cox_summary <- summary(cox_model)
  HR <- exp(cox_summary$coefficients[1, 1])
  p_value_cox <- cox_summary$coefficients[1, 5]
  
  # Manually assign the formula to the km object for the plot
  km$call$formula <- formula_obj
  # Plot voor figuur
  p <- ggsurvplot(km, 
                  data = filter_df,
                  pval = FALSE, 
                  ylab = "Exacerbation rate", 
                  xlab = "Time (days)",
                  risk.table = "nrisk_cumcensor", 
                  conf.int = TRUE, 
                  fun = "event",
                  legend.title = sad_variable,
                  legend.labs = levels(filter_df[[sad_variable]]))
  p$plot <- p$plot +
    ggplot2::annotate("text", x = 10, y = 0.25,
                      label = paste("HR =", round(HR, 2),
                                    "\nLog-rank p =", signif(p_value_logrank, 3)),
                      size = 4, color = "black", hjust = 0) +
    ggtitle(paste0(subgroup))
  file = file.path(cfg$paths$output_severe_sad, paste0("uni_model_", subgroup, "_", sad_variable, ".png"))
  png(file, width = 1200, height = 1000, res= 150)
  print(p)
  dev.off()
  print(paste0("Saved plot ", file))
  return(NULL)
}

#########################################################
sad_variables <- c("R520_higher_3rsd", "X5_lower_3rsd", "AX_higher_3rsd")
asthma_control <- unique(db_exacerbations$ACQ6group)

for (sad_definition in sad_variables) {
  for (control_subgroup in asthma_control)
  {
    sad_survival_uni(sad_definition,
                     db_exacerbations, 
                     control_subgroup)
  }}
