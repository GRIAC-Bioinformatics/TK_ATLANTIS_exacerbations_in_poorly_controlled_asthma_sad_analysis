# This code will use multivariate cox regression to predict time to exacerbation in 3 asthma control groups 
# using 3 parameters to define severe SAD

#########################################################
# load packages
library(survival)
library(survminer)
library(dplyr)
library(broom)
library(yaml)
library(here)

cfg <- yaml::read_yaml(here("config.yaml"))

#########################################################
### load the clinical table + exacerbation + + severe sad definition table 
db <-  read.csv(cfg$paths$patient_data, header = TRUE) %>%
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
# define athma control groups 
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
    
  )

#########################################################
# function to perform multivariate cox regression
sad_survival_multi <- function(sad_variable, db_exacerbations, subgroup = NULL) {
  # prepare df
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
  formula_obj <- as.formula(
    paste0("Surv(time, exacerbation) ~ AGE + SEX + ex_smoker + GINA45 + more_than_1_exac_last_year + LABEOSV + B_FEV1PNVG + ", 
           sad_variable)
  )
  # run cox regression
  cox_modelbasis <- coxph(formula_obj, data = filter_df)
  # tidy the table
  model_table <- tidy(cox_modelbasis) %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    mutate(term = recode(term, "more_than_1_exac_last_yearTRUE" = "1+exac last year"),
           term = recode(term, "LABEOSV" = "Blood eosinophils"),
           term = recode(term, "B_FEV1PNVG" = "FEV1 % predicted"))
  # save models: 
  readr::write_delim(model_table,
                     file.path(cfg$paths$output_severe_sad, 
                               paste0(subgroup, "_", sad_variable, "_multivar_cox.csv")), 
                     delim = ";")
}


#########################################################
sad_variables <- c("R520_higher_3rsd", "X5_lower_3rsd", "AX_higher_3rsd")
asthma_control <- unique(db_exacerbations$ACQ6group)

for (sad_definition in sad_variables) {
  for (control_subgroup in asthma_control)
  {
    sad_survival_multi(sad_definition,
                     db_exacerbations, 
                     control_subgroup)
  }}s