# This code will use multivariate Negative Binomial Model to predict exacerbation frequencies in 3 asthma control groups 
# using 3 parameters to define SAD: R5-20, X5, AX

library(dplyr)
library(broom)
library(tidyr)
library(MASS)
library(openxlsx)

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

# time - days of follow up 
# NUM_EX_D

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

#########################################################
# function to perform multivariate negative binomial with Offset
#########################################################

sad_nb_multi <- function(sad_variable, db_exacerbations, subgroup = NULL) {
  print(paste("calculating", subgroup, sad_variable))
  # prepare df
  # no NAs for the sad variable of interest:
  filter_df <- db_exacerbations %>%
    filter(!is.na(.data[[sad_variable]])) %>%
    mutate(!!sym(sad_variable) := as.factor(!!sym(sad_variable)))
  # no NAs in the other variables for the model:
  vars_model <- c("time",
                  "NUM_EX_D",
                  "AGE",
                  "SEX",
                  "ex_smoker",
                  "current_smoker",
                  "GINA45",
                  "more_than_1_exac_last_year",
                  "LABEOSV", 
                  "B_FEV1PNVG",
                  "B_RVTLC")
  filter_df <- filter_df %>%
    drop_na(all_of(vars_model))
  # filter for asthma control group - if defined 
  if (!is.null(subgroup)) {
    filter_df <- filter_df %>%
      filter(.data[["ACQ6group"]] == subgroup)
    
    message(paste("Filtering for", subgroup))
  }

  # Print the size of the dataset
  print(paste0("SAD variable: ", sad_variable))
  print(paste0("N = ", as.character(nrow(filter_df))))
  
  # association between SAD and exacerbations using NB model and offset - time of the follow-up
  # define the formula
  formula_obj <- as.formula(
    paste0("NUM_EX_D ~ AGE + SEX + ex_smoker + current_smoker + GINA45 + more_than_1_exac_last_year + LABEOSV + B_FEV1PNVG + B_RVTLC + ",
           sad_variable, "+ offset(log(time))"))
  # fit the model
  model_nb <- glm.nb(formula_obj, data = filter_df)
  model_summary <- summary(model_nb)
  
  coef_table <- as.data.frame(coef(summary(model_nb)))
  coef_table <- coef_table %>%
    tibble::rownames_to_column("term") %>%
    filter(grepl(sad_variable, term)) %>% 
    mutate(
      aRR     = round(exp(Estimate), 2),
      CI_low  = round(exp(Estimate - 1.96 * `Std. Error`), 2),
      CI_high = round(exp(Estimate + 1.96 * `Std. Error`), 2),
      CI      = paste0(CI_low, " – ", CI_high),
      p_value = round(`Pr(>|z|)`, 3)
    ) %>%
    mutate(
      sad_var = sad_variable,               # which SAD variable
      subgroup     = ifelse(is.null(subgroup), "All", subgroup),  # which subgroup
      n            = nrow(filter_df))             # sample size
  return(coef_table)
}
    
#########################################################
## run the model for 3 asthma groups and 3 SAD parameters: 

sad_variables <- c("B_R520ABNRp", "B_AXABNRp", "B_X5ABNRp")
asthma_control <- unique(db_exacerbations$ACQ6group)

results_list <- list()  # store results here

for (sad_definition in sad_variables) {
  for (control_subgroup in asthma_control) {
    results_list[[paste(sad_definition, control_subgroup, sep = "_")]] <- 
      sad_nb_multi(sad_definition, db_exacerbations, control_subgroup)
  }
  }

# Combine all results into one clean table
results_table <- bind_rows(results_list)
results_table

## change names for the table 
results_table <- results_table %>%
  mutate(term = case_when(
    term == "B_R520ABNRpTRUE" ~ "R5-R20 above ULN",
    term == "B_AXABNRpTRUE" ~ "AX above ULN", 
    term == "B_X5ABNRpTRUE" ~ "X5 below LLN"
  ))

# save table
write.xlsx(results_table,
           file = file.path(cfg$paths$output_tables, "sad_multi_model", 
                            "sad_multi_NB_model_all_sad_all_asthma.xlsx"),
           rowNames = FALSE)


