## check the definition of ULL/LLN 
# equations are taken from supplementary to main ATLANTIS paper

library(dplyr)
library(yaml)
library(here)

cfg <- yaml::read_yaml(here("config.yaml"))

db <-  read.csv(cfg$paths$patient_data, header = TRUE) %>%
  filter(VISIT == "VISIT 1", f_eval == "Y")

# R520 ULN -CHECK
RSD_M_R520 <- 0.039838
RSD_F_R520 <- 0.039838
db <- db %>%
  mutate(R520ULN = ifelse(
  SEX == "M",
  (0.270439 + (AGE * 0.000285) - (HT * 0.001378)) +  1.645 * RSD_M_R520,
  (0.270439 + (AGE * 0.000285) - (HT * 0.001378) - 0.037831) + 1.645 * RSD_F_R520)) %>%
  mutate(R520_hoger_ULN = if_else(!is.na(R520ULN) & (B_R520 >= R520ULN), "Y", "N"))

db[which(!(db$R520_hoger_ULN == db$B_R520ABNR)), c("B_R520", "R520ULN", "B_R520ABNR", "R520_hoger_ULN")]

# 1. define severe SAD - R5-20 : 3 * RSD 
RSD_M_R520 <- 0.039838
RSD_F_R520 <- 0.039838
db <- db %>%
  mutate(R520ULN = ifelse(
    SEX == "M",
    (0.270439 + (AGE * 0.000285) - (HT * 0.001378)) +  3 * RSD_M_R520,
    (0.270439 + (AGE * 0.000285) - (HT * 0.001378) - 0.037831) + 3 * RSD_F_R520)) %>%
  mutate(R520_higher_3rsd = if_else(!is.na(R520ULN) & (B_R520 >= R520ULN), "Y", "N"))
table(db$R520_higher_3rsd)

# 2. define severe SAD - X5 LLN : 3 * RSD 
RSD_M_X5 <- 0.040651
RSD_F_X5 <- 0.040651
db <- db %>%
  mutate(X5LLN = ifelse(
    SEX == "M",
    (-0.473195 - (AGE * 0.000867) - (HT * 0.002408)) +  3 * RSD_M_X5,
    (-0.473195 + (AGE * 0.000867) - (HT * 0.002408) - 0.007246 + 3) * RSD_F_X5)) %>%
  mutate(X5_lower_3rsd = if_else(!is.na(X5LLN) & (B_X5 <= X5LLN), "Y", "N"))

db[, c("X5LLN", "B_X5", "X5_lower_3rsd")]
table(db$X5_lower_3rsd)


# 3. define severe SAD - AX ULN : 3 * RSD 
RSD_M_AX <- 0.326243
RSD_F_AX <- 0.326243

db <- db %>%
  mutate(AXULN = ifelse(
    SEX == "M",
    (2.426203 + (AGE * 0.004599) - (HT * 0.013162)) + 3 * RSD_M_AX,
    (2.426203 + (AGE * 0.004599) - (HT * 0.013162) - 0.088393) + 3 * RSD_F_AX)) %>%
  mutate(AX_higher_3rsd = if_else(!is.na(AXULN) & (B_AX >= AXULN), "Y", "N"))

db[, c("AXULN", "B_AX", "AX_higher_3rsd")]
table(db$AX_higher_3rsd)

write.csv(db %>%
            dplyr::select(c("R520_higher_3rsd", "X5_lower_3rsd", "AX_higher_3rsd", "PT")),
          cfg$paths$severe_sad_table, row.names = F)


