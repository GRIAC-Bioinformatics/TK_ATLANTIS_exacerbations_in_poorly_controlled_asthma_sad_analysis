
library(dplyr)
library(stringr)
library(tidyr)
library(yaml)
library(here)

cfg <- yaml::read_yaml(here("config.yaml"))

#Datasets inladen
db <- read.csv(cfg$paths$patient_data, header = TRUE) %>%
  filter(VISIT == "VISIT 1", ASTHEA == "A", f_eval == "Y") %>%
  mutate(extrafine_ics_all = if_else((EXTRAFINE_ICS == "Yes" | EXTRAFINE_ICS_LABA == "Yes"), "Yes", "No"),
         nonextrafine_ics_all = if_else((NONEXTRAFINE_ICS == "Yes" | NONEXTRAFINE_ICS_LABA == "Yes"), "Yes", "No")) %>%
  mutate(ics_treatment = case_when(
    extrafine_ics_all == "Yes" & nonextrafine_ics_all == "No" ~ "extrafine",
    extrafine_ics_all == "No" & nonextrafine_ics_all == "Yes" ~ "nonextrafine",
    extrafine_ics_all == "No" & nonextrafine_ics_all == "No" ~ "none",
    extrafine_ics_all == "Yes" & nonextrafine_ics_all == "Yes" ~ "both",
    .default = NA
  )) %>%
  mutate(GINA45 = GINA %in% c("4", "5"),
         GINA12 = GINA %in% c("1", "2"),
         PHADRES_N = as.factor(PHADRES_N))

# Evt filteren op astma controle
# db <- subset.data.frame(db, db$acq6_score < 0.75)  # well contorlled
# db <- subset.data.frame(db, db$acq6_score >= 0.75 & db$acq6_score <= 1.5)  # partially controlled
# db <- subset.data.frame(db, db$acq6_score > 1.5)   # poorly controlled

#variabelen aanmaken
#bhr
db <- db %>%
  mutate(
    bhr_c = cut(
      PC20, #concentraion 
      breaks = c(0.0000, 0.25, 1, 4, 16, Inf),
      labels = c("Severe", "Moderate", "Mild", "Very mild", "Normal") # concentration at 20% drop
    )
  )

db <- db %>%
  mutate(
    bhr_d = cut(
      PD20, #dose
      breaks = c(0.0000, 0.03, 0.13, 0.5, 2, Inf),
      labels = c("Severe", "Moderate", "Mild", "Very mild", "Normal")
    )
  )

db <- db %>%
  tidyr::unite(bhr, c(bhr_c, bhr_d), sep = "_", remove = TRUE) %>% 
  mutate(
    bhr = str_replace_all(bhr, "^NA_|_NA$", ""),   # remove leading "NA_" or trailing "_NA"
    bhr = na_if(bhr, "NA"),                         # convert "NA" string to real NA
    bhr = factor(bhr, levels = c("Very mild", "Mild", "Moderate", "Severe")), # factor with levels
    bhrseveremoderate = bhr %in% c("Moderate", "Severe"),  # Moderate/Severe flag
    any_ICS = if_else(ICS == "No" & ICS_LABA == "No", "No", "Yes"),
    ICS_dose_sum = if_else(is.na(ICS_DDOSE_EQ) & is.na(ICS_LABA_DDOSE_EQ), NA,
                                    coalesce(ICS_DDOSE_EQ, 0) + coalesce(ICS_LABA_DDOSE_EQ, 0)),
    reversibility_FEV1 = (C_FEV1 - B_FEV1) / B_FEV1 * 100 # reversibility calculation (the same as DFEV1P)
  )  
  

# Asthma controle groep aanmaken
db <- db %>%
  mutate(
    ACQ6group = case_when(
      acq6_score < 0.75                  ~ "well controlled asthma",
      acq6_score >= 0.75 & acq6_score <= 1.5 ~ "partially controlled asthma",
      acq6_score > 1.5                   ~ "poorly controlled asthma",
      TRUE                                ~ NA_character_
    ),
    ACQ6group = factor(ACQ6group)  # convert to factor
  )


# SAD variabelen (zodat NA niet wordt meegenomen)
db$B_R520ABNRp <- ifelse(db$B_R520ABNR == "Y", TRUE,
                         ifelse(db$B_R520ABNR == "N", FALSE, NA))
db$B_AXABNRp   <- ifelse(db$B_AXABNR == "Y", TRUE,
                         ifelse(db$B_AXABNR == "N", FALSE, NA))
db$B_X5ABNRp   <- ifelse(db$B_X5ABNR == "Y", TRUE,
                         ifelse(db$B_X5ABNR == "N", FALSE, NA))


db <- db %>%
  mutate(
    B_R520ABNRp = case_when(
      B_R520ABNR == "Y" ~ TRUE,
      B_R520ABNR == "N" ~ FALSE,
      B_R520ABNR == ""  ~ NA
    ),
    B_AXABNRp = case_when(
      B_AXABNR == "Y" ~ TRUE,
      B_AXABNR == "N" ~ FALSE,
      B_AXABNR == ""  ~ NA
    ),
    B_X5ABNRp = case_when(
      B_X5ABNR == "Y" ~ TRUE,
      B_X5ABNR == "N" ~ FALSE,
      B_X5ABNR == ""  ~ NA
    ),
    `Current or former smokers` = case_when(
      SMOKE == "Current Smoker" | SMOKE == "Ex-smoker"  ~ TRUE,
      SMOKE == "Non-smoker" ~ FALSE,
    ),
    GINA = as.factor(GINA),
    `More than 1 asthma exacerbations (last year)` = if_else(NUM_EX>0, "Yes", "No"),
    `More than 1 asthma exacerbations (during study)` = if_else(NUM_EX_D>0, "Yes", "No"),
  )

# rename columns:
db <- db %>%
  rename("Age, y" = "AGE",
         "Body mass index, kg/m2" = "BMI",
         "Number of pack-years" = "SPACKNO",
         "Asthma duration, y" = "DUR_DIS",
         "Presence of atopy" = "PHADRES",
         "Blood eosinophils, 10^9/L" = "LABEOSV",
         "Sputum eosinophils, (% of non-squamous cells)" = "EOSP",
         "FeNO, ppb" = "FENRES",
         "FEV1 (% predicted)" = "B_FEV1PNVG", #pre bronchodilator
         "FEV1/FVC (% predicted)" = "B_FEV1FPNVG",
         "Bronchodilator reversibility, % initial" = "DFEV1P",
         "FEF50 (% predicted)" = "B_F50PNVR",
         "FEF25-75 (% predicted)" = "B_F2575PNVG",
         "Resistance between 5 and 20 Hz" = "B_R520",
         "Resistance at 5 Hz" = "B_R5",
         "Resistance between 5 and 20 Hz, (% predicted)" = "B_R520PNVR",
         "Asthma exacerbations (last year)" = "NUM_EX",
         "Asthma exacerbations (during the study)" = "NUM_EX_D",
         "Inhaled corticosteroid dosage, (beclomethasone equivalent)*, ug" = "ICS_dose_sum",
         "Inhaled corticosteroid type" = "ics_treatment",
         "Systemic corticosteroids, n (%)" = "SYS_COR",
         "R5-R20 above ULN" ="B_R520ABNRp",
         "AX above ULN" = "B_AXABNRp",
         "X5 below LLN" = "B_X5ABNRp")
         
         
         
#Baseline tabel aanmaken
library(tableone)
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
  "GINA45",
  "Presence of atopy",
  "PHADRES_N",
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
  "any_ICS",
  "Systemic corticosteroids, n (%)")
  
nietnormaalverdeeld <- c( "FeNO, ppb", "Inhaled corticosteroid dosage, (beclomethasone equivalent)*, μg",
                          "AGE_DIAG",
                          "Asthma exacerbations (last year)",
                          "Asthma exacerbations (during the study)",
                          "Blood eosinophils, 10^9/L",
                          "Sputum eosinophils, (% of non-squamous cells)",
                          "B_FEV1FPNVG")


ACQ6baseline <- CreateTableOne(vars = variabelen, strata = "ACQ6group", data = db)
Final_statistics <- print(ACQ6baseline, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)

write.table(Final_statistics, file = file.path(cfg$paths$output_tables, "asthma_control_baseline.csv"),
          sep = "\t")

db_well_poor <- db %>%
  filter(ACQ6group != "partially controlled asthma") %>%
  mutate(ACQ6group = droplevels(ACQ6group))

ACQ6baseline_well_poor <- CreateTableOne(vars = variabelen, strata = "ACQ6group", data = db_well_poor)
Final_statistics_well_poor <- print(ACQ6baseline_well_poor, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)
write.table(Final_statistics_well_poor, file = file.path(cfg$paths$output_tables, "asthma_control_well_poor_baseline.csv"),
            sep = "\t")
  
### stratify by SAD 
# R5-20

bd_R5_20 <- db %>%
  filter(!is.na(`R5-R20 above ULN`)) %>%
  mutate (SAD_R5_20 = if_else(`R5-R20 above ULN` == "TRUE", paste0(ACQ6group, " ", "R5-R20 above ULN"),
                              paste0(ACQ6group, "", "R5-R20 below ULN")))

ACQ6baseline <- CreateTableOne(vars = variabelen, strata = "SAD_R5_20", data = bd_R5_20)
Final_statistics <- print(ACQ6baseline, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)

write.table(Final_statistics, file = file.path(cfg$paths$output_tables, "asthma_control_baseline_R5_20strata.csv"),
            sep = "\t")

  
# AX

bd_AX <- db %>%
  filter(!is.na(`AX above ULN`)) %>%
  mutate (SAD_AX = if_else(`AX above ULN` == "TRUE", paste0(ACQ6group, " ", "AX above ULN"),
                              paste0(ACQ6group, "", "AX below ULN")))

ACQ6baseline <- CreateTableOne(vars = variabelen, strata = "SAD_AX", data = bd_AX)
Final_statistics <- print(ACQ6baseline, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)

write.table(Final_statistics, file = file.path(cfg$paths$output_tables, "asthma_control_baseline_AXstrata.csv"),
            sep = "\t")

# X5

bd_X5 <- db %>%
  filter(!is.na(`X5 below LLN`)) %>%
  mutate (SAD_X5 = if_else(`X5 below LLN` == "TRUE", paste0(ACQ6group, " ", "X5 below LLN"),
                           paste0(ACQ6group, "", "X5 above LLN")))

ACQ6baseline <- CreateTableOne(vars = variabelen, strata = "SAD_X5", data = bd_X5)
Final_statistics <- print(ACQ6baseline, nonnormal = nietnormaalverdeeld, showAllLevels = FALSE, pDigits = 5)

write.table(Final_statistics, file = file.path(cfg$paths$output_tables, "asthma_control_baseline_X5strata.csv"),
            sep = "\t")
