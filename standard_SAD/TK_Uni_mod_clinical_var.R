library(survival)
library(survminer)
library(dplyr)
library(here)
library(yaml)

cfg <- yaml::read_yaml(here("config.yaml"))

#dataset inladen
load(cfg$paths$exacerbations_data) # from Tessa K.

db1 <- joined %>%
  mutate(PT = as.numeric(PT)) %>%
  dplyr::select(c("PT", "date_baseline", "date_exacerbation", "time_to_exacerbation",
                  "exacerbation", "time")) %>%
  mutate(exacerbation = as.integer(exacerbation)) %>%
  mutate(time = if_else(time>420, 420, time))


db2 <-  read.csv(cfg$paths$patient_data, header = TRUE) %>%
  filter(VISIT == "VISIT 1", ASTHEA == "A", f_eval == "Y") 

db <- merge(db1, db2, by = "PT", all.y = TRUE)

db <- db %>%
  mutate(B_R520ABNRp = ifelse(B_R520ABNR == "Y", TRUE, ifelse(B_R520ABNR == "N", FALSE, NA)),
         B_AXABNRp   = ifelse(B_AXABNR == "Y", TRUE, ifelse(B_AXABNR == "N", FALSE, NA)),
         B_X5ABNRp   = ifelse(B_X5ABNR == "Y", TRUE, ifelse(B_X5ABNR == "N", FALSE, NA)),
         R520_AX     = ifelse(B_R520ABNRp & B_AXABNRp, TRUE, ifelse(!B_R520ABNRp & !B_AXABNRp, FALSE, NA)),
         R520_X5     = ifelse(B_R520ABNRp & B_X5ABNRp, TRUE, ifelse(!B_R520ABNRp & !B_X5ABNRp, FALSE, NA)),
         AX_X5       = ifelse(B_AXABNRp & B_X5ABNRp, TRUE, ifelse(!B_AXABNRp & !B_X5ABNRp, FALSE, NA)),
         AX_X5_R520  = ifelse(B_AXABNRp & B_X5ABNRp & B_R520ABNRp, TRUE,
                              ifelse(!B_AXABNRp & !B_X5ABNRp & !B_R520ABNRp, FALSE, NA)),
         last_year_ex = ifelse(NUM_EX >= 1, "Yes", "No"),
         future_ex = ifelse(NUM_EX_D >= 1, "Yes", "No"),
         GINA45 = ifelse(GINA >= 4, "Yes", "No"),
         ex_smoker = if_else(SMOKE == "Ex-smoker", "Yes", "No"),
         current_smoker = if_else(SMOKE == "Current Smoker", "Yes", "No"),
         never_smoker = if_else(SMOKE == "Non-smoker", "Yes", "No"),
         ACQ6group = case_when(
           acq6_score < 0.75 ~ "well controlled asthma",
           acq6_score >= 0.75 & acq6_score <= 1.5 ~ "partially controlled asthma",
           acq6_score > 1.5 ~ "poorly controlled asthma",
           TRUE  ~ NA_character_),
           ACQ6group = factor(ACQ6group)) %>%
  filter(!is.na(exacerbation))

# Check clinical variables in a univariate model
vars_model <- c("SEX","GINA45","ex_smoker","current_smoker",
                "never_smoker", "B_FEV1PNVG","LABEOSV", "FENRES",
                "last_year_ex",  "B_R520ABNRp", "B_AXABNRp", "B_X5ABNRp")


results <- purrr::map_dfr(vars_model, function(v) {
  fml <- as.formula(paste0("Surv(time, exacerbation) ~ ", v))
  fit <- coxph(fml, data = db)
  broom::tidy(fit) %>% mutate(variable = v)
}) %>%
  mutate(across(where(is.numeric), signif, 2))

write.csv(results, file.path(cfg$paths$output_tables,"all_asthma_uni_clinical.csv"))

# divide asthma controls: 
results_WC <- purrr::map_dfr(vars_model, function(v) {
  fml <- as.formula(paste0("Surv(time, exacerbation) ~ ", v))
  fit <- coxph(fml, data = db[db$ACQ6group == "well controlled asthma",])
  broom::tidy(fit) %>% mutate(variable = v) %>%
  mutate(across(where(is.numeric), signif, 2))
})
write.csv(results_WC, file.path(cfg$paths$output_tables,"WC_asthma_uni_clinical.csv"))


results_partially <- purrr::map_dfr(vars_model, function(v) {
  fml <- as.formula(paste0("Surv(time, exacerbation) ~ ", v))
  fit <- coxph(fml, data = db[db$ACQ6group == "partially controlled asthma",])
  broom::tidy(fit) %>% mutate(variable = v) %>%
  mutate(across(where(is.numeric), signif, 2))
})
write.csv(results_partially, file.path(cfg$paths$output_tables,"partially_asthma_uni_clinical.csv"))

results_poorly <- purrr::map_dfr(vars_model, function(v) {
  fml <- as.formula(paste0("Surv(time, exacerbation) ~ ", v))
  fit <- coxph(fml, data = db[db$ACQ6group == "poorly controlled asthma",])
  broom::tidy(fit) %>% mutate(variable = v) %>%
  mutate(across(where(is.numeric), signif, 2))
})
write.csv(results_poorly, file.path(cfg$paths$output_tables, "poorly_asthma_uni_clinical.csv"))


