
library(dplyr)
library(stringr)
library(tidyr)
library(tableone)
library(yaml)
library(here)

cfg <- yaml::read_yaml(here("config.yaml"))

source(here("R", "baseline_table_source.R"))

db <-  read.csv(cfg$paths$patient_data, header = TRUE) %>%
  filter(VISIT == "VISIT 1", ASTHEA == "A", f_eval == "Y")

## format the table - rename columns, change factors
db <- column_rename(db)

## add severe SAD info- columns "R520_higher_3rsd", "X5_lower_3rsd", "AX_higher_3rsd"
db <- db %>%
  left_join(read.csv(cfg$paths$severe_sad_table), by = c("PT"="PT"))

## prepare for the tableone
variabelen <- c(
  "RACE", "GINA12", 
  "Age, y", 
  "SEX", 
  "Body mass index, kg/m2", 
  "Current or former smokers",
  "Number of pack-years",
  "Asthma duration, y", 
  "AGE_DIAG",
  "GINA",
  "Presence of atopy",
  "Blood eosinophils, 10^9/L",
  "Sputum eosinophils, (% of non-squamous cells)",
  "FeNO, ppb", 
  "FEV1 (% predicted)",
  "FEV1/FVC (% predicted)",
  "Bronchodilator reversibility, % initial",
  "FEF50 (% predicted)",
  "FEF25-75 (% predicted)",
  "Resistance between 5 and 20 Hz",
  "Resistance at 5 Hz",
  "Resistance between 5 and 20 Hz, (% predicted)",
  "R5-R20 above ULN", "AX above ULN", "X5 below LLN",
  "Asthma exacerbations (last year)",
  "More than 1 asthma exacerbations (last year)",
  "Asthma exacerbations (during the study)",
  "More than 1 asthma exacerbations (during study)",
  "Inhaled corticosteroid dosage, (beclomethasone equivalent)*, μg",
  "Inhaled corticosteroid type",
  "Systemic corticosteroids, n (%)")

nietnormaalverdeeld <- c( "FeNO, ppb", "Inhaled corticosteroid dosage, (beclomethasone equivalent)*, μg",
                          "AGE_DIAG",
                          "Asthma exacerbations (last year)",
                          "Asthma exacerbations (during the study)",
                          "Blood eosinophils, 10^9/L",
                          "Sputum eosinophils, (% of non-squamous cells)")

### stratify by SAD 
# severe R5-20

db_R5_20 <- db %>%
  filter(!is.na(R520_higher_3rsd)) %>%
  mutate (severe_SAD_R5_20 = if_else(R520_higher_3rsd == "Y", paste0(ACQ6group, " ", R520_higher_3rsd),
                              paste0(ACQ6group, " ", R520_higher_3rsd)))

R520 <- CreateTableOne(vars = variabelen, strata = "severe_SAD_R5_20", data = db_R5_20)
Final_statistics <- print(R520, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)

write.table(Final_statistics, file = file.path(cfg$paths$output_tables, "asthma_control_baseline_severe_R5_20strata.csv"),
            sep = "\t")

# severe AX
db_AX <- db %>%
  filter(!is.na(AX_higher_3rsd)) %>%
  mutate (severe_SAD_AX = if_else(AX_higher_3rsd == "Y", paste0(ACQ6group, " ", AX_higher_3rsd),
                                     paste0(ACQ6group, " ", AX_higher_3rsd)))

AX <- CreateTableOne(vars = variabelen, strata = "severe_SAD_AX", data = db_AX)
Final_statistics <- print(AX, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)

write.table(Final_statistics, file = file.path(cfg$paths$output_tables, "asthma_control_baseline_severe_AXstrata.csv"),
            sep = "\t")

# severe X5
db_X5 <- db %>%
  filter(!is.na(X5_lower_3rsd)) %>%
  mutate (severe_SAD_X5 = if_else(X5_lower_3rsd == "Y", paste0(ACQ6group, " ", X5_lower_3rsd),
                                  paste0(ACQ6group, " ", X5_lower_3rsd)))

X5 <- CreateTableOne(vars = variabelen, strata = "severe_SAD_X5", data = db_X5)
Final_statistics <- print(AX, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)

write.table(Final_statistics, file = file.path(cfg$paths$output_tables, "asthma_control_baseline_severe_X5strata.csv"),
            sep = "\t")

