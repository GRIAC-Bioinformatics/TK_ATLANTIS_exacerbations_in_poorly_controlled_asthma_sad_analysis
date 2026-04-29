# This code will use univariate cox regression to predict time to exacerbation in 3 asthma control groups 
# using 3 parameters to define SAD: R5-20, X5, AX

#########################################################
# load packages
library(survival)
library(survminer)
library(dplyr)
library(broom)
library(tidyr)
library(ggplotify)
library(here)
library(yaml)

cfg <- yaml::read_yaml(here("config.yaml"))

#########################################################
### load the clinical table + exacerbation 
db <-  read.csv(cfg$paths$patient_data, header = TRUE) %>%
  filter(VISIT == "VISIT 1", ASTHEA == "A", f_eval == "Y")

load(cfg$paths$exacerbations_data) # from Tessa K.
db_exacerbations <- joined %>%
  mutate(PT = as.numeric(PT)) %>%
  dplyr::select(c("PT", "date_baseline", "date_exacerbation", "time_to_exacerbation",
                  "exacerbation", "time")) %>%
  mutate(exacerbation = as.integer(exacerbation)) %>%
  mutate(time = if_else(time>420, 420, time)) %>%
  left_join(db, by = c("PT" = "PT"))

#########################################################
## prepare table:
# define asthma control groups 
db_exacerbations <- db_exacerbations %>% 
  mutate(
    ACQ6group = case_when(
      acq6_score < 0.75                  ~ "well controlled asthma",
      acq6_score >= 0.75 & acq6_score <= 1.5 ~ "partially controlled asthma",
      acq6_score > 1.5                   ~ "poorly controlled asthma",
      TRUE                                ~ NA_character_
    ),
    ACQ6group = factor(ACQ6group), # convert to factor
    GINA45 = GINA %in% c("4", "5"),
    GINA12 = GINA %in% c("1", "2"),
    ex_smoker = SMOKE == "Ex-smoker",
    current_smoker = SMOKE == "Current Smoker",
    more_than_1_exac_last_year = NUM_EX >= 1,
    B_R520ABNRp = ifelse(B_R520ABNR == "Y", TRUE, ifelse(B_R520ABNR == "N", FALSE, NA)),
    B_AXABNRp   = ifelse(B_AXABNR == "Y", TRUE, ifelse(B_AXABNR == "N", FALSE, NA)),
    B_X5ABNRp   = ifelse(B_X5ABNR == "Y", TRUE, ifelse(B_X5ABNR == "N", FALSE, NA))
  )

# db$B_R520ABNRp <- factor(filter_df$B_R520ABNRp, levels = c(FALSE, TRUE),
#                          labels = c("R5-R20 below ULN", "R5-R20 above ULN"))
# db$B_AXABNRp <- factor(db$B_AXABNRp, levels = c(FALSE, TRUE),
#                        labels = c("AX below ULN", "AX above ULN"))
# Factorlabels (TRUE = below LLN, FALSE = above LLN)
# db$B_X5ABNRp <- factor(db$B_X5ABNRp, levels = c(FALSE, TRUE),
#                        labels = c("X5 above LLN", "X5 below LLN"))
# 

#########################################################
# function to perform univariate cox regression and make a plot
sad_survival_uni <- function(sad_variable, db_exacerbations, subgroup = NULL) {
  print(paste("calculating", subgroup, sad_variable))
  # prepare df
  # no NAs for the sad variable of interest:
  filter_df <- db_exacerbations %>%
    filter(!is.na(.data[[sad_variable]])) %>%
    mutate(!!sym(sad_variable) := as.factor(!!sym(sad_variable)))
  
  # filter for asthma control group - if defined 
  if (!is.null(subgroup)) {
    filter_df <- filter_df %>%
      filter(.data[["ACQ6group"]] == subgroup)
    
    message(paste("Filtering for", subgroup))
  }
  # Print the size of the dataset
  print(paste0("SAD variable: ", sad_variable))
  print(paste0("N = ", as.character(nrow(filter_df))))
  
  # association between SAD and exacerbations
  
  s <- Surv(time = filter_df$time, event = filter_df$exacerbation,
            type = "right") # time to event if TRUE/ time follow-up if false
  formula_obj <- as.formula(paste0("s ~ ", sad_variable))
  logrank_test <- survdiff(formula_obj, data = filter_df) #Tests if there is a difference between two or more survival curves
  p_value_logrank <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)

  # HR + P-value
  cox_model <- coxph(formula_obj, data = filter_df)
  cox_summary <- summary(cox_model)
  HR <- exp(cox_summary$coefficients[1, 1])
  p_value_cox <- cox_summary$coefficients[1, 5]
  

  ## change levels for the plot 
  filter_df <- filter_df %>%
    mutate(
      B_R520ABNRp = as.factor(if_else(B_R520ABNRp == TRUE, "R5-R20 above ULN", "R5-R20 below ULN")),
      B_AXABNRp   = as.factor(if_else(B_AXABNRp == TRUE, "AX above ULN", "AX below ULN")),
      B_X5ABNRp   = as.factor(if_else(B_X5ABNRp == TRUE, "X5 below LLN", "X5 above LLN"))
    )
  # Ensure the factor levels are clean (removes 'ghost' levels not in the data)
  filter_df[[sad_variable]] <- droplevels(filter_df[[sad_variable]])
  # Extract current labels safely
  current_labels <- levels(filter_df[[sad_variable]])
  
  # plot 
  # KM fit for the plot
  km <- survfit(formula_obj, data = filter_df)
  # MANUALLY ASSIGN THE FORMULA - for the plot
  km$call$formula <- formula_obj
  km$call$data <- filter_df
  
  print("passed prepare")
  p <- ggsurvplot(km, pval = FALSE, 
                  ylab = "Cumulative incidence of exacerbation", 
                  xlab = "Time (days)",
                  ylim = c(0, 0.8),
                  risk.table = "nrisk_cumcensor", conf.int = TRUE, fun = "event",
                  legend.labs = levels(filter_df[[sad_variable]]))
  print("passed p")

  p_inset <- ggsurvplot(km, pval = FALSE,
                        conf.int = TRUE,
                        fun = "event",
                        legend = "none",
                        risk.table = FALSE,
                        xlab = "",
                        ylab = "",
                        legend.labs = levels(filter_df[[sad_variable]]))
  inset_grob <- as.grob(p_inset$plot +
                          theme(axis.text = element_text(size = 6),
                                axis.title = element_blank(),
                                plot.background = element_rect(color = "black", linewidth = 0.5)))
  # Add inset to main plot using annotation_custom
  p$plot <- p$plot +
    annotation_custom(
      grob = inset_grob,
      xmin = -Inf, xmax = 200,   # e.g. xmax = 100 if x goes to 300
      ymin = 0.45, ymax = 0.8                               # top-left area
    )

  # add annotation
  p$plot <- p$plot + ggplot2::annotate("text", x = 250, y = 0.85,
                                       label = paste0("HR = ", round(HR, 2),
                                                    " (95% CI ", as.character(round(cox_summary$conf.int[,"lower .95"], 2)),
                                                    "-", as.character(round(cox_summary$conf.int[,"upper .95"], 2)),
                                                    ");",
                                                    "\nLog-rank p = ", round(p_value_logrank, 3)),
                                       size = 4, color = "black", hjust = 0)
  png(file.path(cfg$paths$output_figures, "sad_uni_model", paste0("uni_surv_", subgroup, "_", sad_variable, ".png")), 
      width = 600, 
      height = 400)
  print(p)
  dev.off()
  # save R object for the Figure 1
  if ((subgroup %in% c("well controlled asthma", "poorly controlled asthma")) & (sad_variable == "B_R520ABNRp")) {
    saveRDS(p, file.path(cfg$paths$output_figures, "sad_uni_model", paste0("uni_surv_", subgroup, "_", sad_variable, ".rds")))
  }
  
 
}

#########################################################
sad_variables <- c("B_R520ABNRp", "B_AXABNRp", "B_X5ABNRp")
asthma_control <- unique(db_exacerbations$ACQ6group)

for (sad_definition in sad_variables) {
  for (control_subgroup in asthma_control)
  {
    sad_survival_uni(sad_definition,
                       db_exacerbations, 
                       control_subgroup)
  }}
