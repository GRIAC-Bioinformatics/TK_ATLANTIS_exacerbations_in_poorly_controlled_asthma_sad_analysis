## source code to prepare ATLANTIS data for the baseline tables 

# rename columns: 
column_rename <- function(db){
  db <- db %>%
  mutate(extrafine_ics_all = if_else((EXTRAFINE_ICS == "Yes" | EXTRAFINE_ICS_LABA == "Yes"), "Yes", "No"),
         nonextrafine_ics_all = if_else((NONEXTRAFINE_ICS == "Yes" | NONEXTRAFINE_ICS_LABA == "Yes"), "Yes", "No")) %>%
  mutate(ics_treatment = case_when(
    extrafine_ics_all == "Yes" & nonextrafine_ics_all == "No" ~ "extrafine",
    extrafine_ics_all == "No" & nonextrafine_ics_all == "Yes" ~ "nonextrafine",
    extrafine_ics_all == "No" & nonextrafine_ics_all == "No" ~ "none",
    extrafine_ics_all == "Yes" & nonextrafine_ics_all == "Yes" ~ "both",
    .default = NA
  )) %>%
    mutate(
      bhr_c = cut(
        PC20, #concentraion 
        breaks = c(0.0000, 0.25, 1, 4, 16, Inf),
        labels = c("Severe", "Moderate", "Mild", "Very mild", "Normal") # concentration at 20% drop
      )
    ) %>%
    mutate(
      bhr_d = cut(
        PD20, #dose
        breaks = c(0.0000, 0.03, 0.13, 0.5, 2, Inf),
        labels = c("Severe", "Moderate", "Mild", "Very mild", "Normal")
      )
    ) %>%
    tidyr::unite(bhr, c(bhr_c, bhr_d), sep = "_", remove = TRUE) %>% 
    mutate(
      bhr = str_replace_all(bhr, "^NA_|_NA$", ""),   # remove leading "NA_" or trailing "_NA"
      bhr = na_if(bhr, "NA"),                         # convert "NA" string to real NA
      bhr = factor(bhr, levels = c("Very mild", "Mild", "Moderate", "Severe")), # factor with levels
      bhrseveremoderate = bhr %in% c("Moderate", "Severe"),  # Moderate/Severe flag
      any_ICS = if_else(ICS == "No" & ICS_LABA == "No", "No", "Yes"),
      ICS_dose_sum = if_else(is.na(ICS_DDOSE_EQ) & is.na(ICS_LABA_DDOSE_EQ), NA,
                             coalesce(ICS_DDOSE_EQ, 0) + coalesce(ICS_LABA_DDOSE_EQ, 0)),
      reversibility_FEV1 = (C_FEV1 - B_FEV1) / B_FEV1 * 100) %>% # reversibility calculation (the same as DFEV1P)
    mutate(
      ACQ6group = case_when(
        acq6_score < 0.75                  ~ "well controlled asthma",
        acq6_score >= 0.75 & acq6_score <= 1.5 ~ "partially controlled asthma",
        acq6_score > 1.5                   ~ "poorly controlled asthma",
        TRUE                                ~ NA_character_
      ),
      ACQ6group = factor(ACQ6group)  # convert to factor
    ) %>%
    mutate(
    `Current or former smokers` = case_when(
      SMOKE == "Current Smoker" | SMOKE == "Ex-smoker"  ~ TRUE,
      SMOKE == "Non-smoker" ~ FALSE,
    ),
    GINA = as.factor(GINA),
    `More than 1 asthma exacerbations (last year)` = if_else(NUM_EX>0, "Yes", "No"),
    `More than 1 asthma exacerbations (during study)` = if_else(NUM_EX_D>0, "Yes", "No"),
    ) %>%
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
      )) %>%
    rename("Age, y" = "AGE",
           "Body mass index, kg/m2" = "BMI",
           "Number of pack-years" = "SPACKNO",
           "Asthma duration, y" = "DUR_DIS",
           "Presence of atopy" = "PHADRES",
           "Blood eosinophils, 10^9/L"= "LABEOSV",
           "Sputum eosinophils, (% of non-squamous cells)"= "EOSP",
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
           "Inhaled corticosteroid dosage, (beclomethasone equivalent)*, μg" = "ICS_dose_sum",
           "Inhaled corticosteroid type" = "ics_treatment",
           "Systemic corticosteroids, n (%)" = "SYS_COR",
           "R5-R20 above ULN" ="B_R520ABNRp",
           "AX above ULN"= "B_AXABNRp", 
           "X5 below LLN" = "B_X5ABNRp")
return(db)
  }