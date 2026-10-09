# ============================================================
# CKD Statistical Risk Assessment Using NHANES
# 01_data_preparation.R
#
# Purpose:
#   Prepare an analysis-ready NHANES 2021-2023 dataset for
#   the CKD risk assessment project.
#
# Inputs:
#   data/raw/DEMO_L.xpt
#   data/raw/ALB_CR_L.xpt
#   data/raw/BIOPRO_L.xpt
#
# Output:
#   data/processed/ckd_analysis_2021_2023.rds
#
# Author: Denis Folitse
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(haven)
library(dplyr)
library(readr)


# ------------------------------------------------------------
# 2. File paths
# ------------------------------------------------------------

demo_path <- "data/raw/DEMO_L.xpt"
albcr_path <- "data/raw/ALB_CR_L.xpt"
biopro_path <- "data/raw/BIOPRO_L.xpt"

output_path <- "data/processed/ckd_analysis_2021_2023.rds"


# ------------------------------------------------------------
# 3. Read NHANES components
# ------------------------------------------------------------

demo <- read_xpt(demo_path)
albcr <- read_xpt(albcr_path)
biopro <- read_xpt(biopro_path)

cat("Raw NHANES component sizes:\n")
cat("DEMO_L:   ", nrow(demo), "\n")
cat("ALB_CR_L: ", nrow(albcr), "\n")
cat("BIOPRO_L: ", nrow(biopro), "\n\n")


# ------------------------------------------------------------
# 4. Restrict demographics to adults age >= 20
# ------------------------------------------------------------

demo_adult <- demo %>%
  filter(RIDAGEYR >= 20)

cat("Adults age >= 20:", nrow(demo_adult), "\n\n")


# ------------------------------------------------------------
# 5. Select variables needed for this project
# ------------------------------------------------------------

demo_keep <- demo_adult %>%
  select(
    SEQN,
    RIDAGEYR,
    RIAGENDR,
    RIDRETH3,
    DMDEDUC2,
    DMDMARTZ,
    INDFMPIR,
    WTMEC2YR,
    SDMVSTRA,
    SDMVPSU
  )

albcr_keep <- albcr %>%
  select(
    SEQN,
    URXUMA,
    URXUCR,
    URDACT
  )

biopro_keep <- biopro %>%
  select(
    SEQN,
    LBXSCR,
    LBXSUA,
    LBXSGL,
    LBXSBU,
    LBXSAL,
    LBXSCH,
    LBXSTR
  )


# ------------------------------------------------------------
# 6. Merge NHANES components
# ------------------------------------------------------------

ckd_data <- demo_keep %>%
  left_join(albcr_keep, by = "SEQN") %>%
  left_join(biopro_keep, by = "SEQN")

stopifnot(nrow(ckd_data) == nrow(demo_keep))


# ------------------------------------------------------------
# 7. Construct eGFR
#
# 2021 CKD-EPI race-free creatinine equation:
#
# eGFR =
#   142 *
#   min(Scr/kappa, 1)^alpha *
#   max(Scr/kappa, 1)^(-1.200) *
#   0.9938^Age *
#   1.012 [if female]
#
# where:
#   kappa = 0.7 for females, 0.9 for males
#   alpha = -0.241 for females, -0.302 for males
#
# NHANES RIAGENDR:
#   1 = Male
#   2 = Female
# Scr = serum creatinine (standardized)
# ------------------------------------------------------------

ckd_data <- ckd_data %>%
  mutate(
    kappa = case_when(
      RIAGENDR == 1 ~ 0.9,
      RIAGENDR == 2 ~ 0.7,
      TRUE ~ NA_real_
    ),
    
    alpha = case_when(
      RIAGENDR == 1 ~ -0.302,
      RIAGENDR == 2 ~ -0.241,
      TRUE ~ NA_real_
    ),
    
    sex_factor = case_when(
      RIAGENDR == 2 ~ 1.012,
      RIAGENDR == 1 ~ 1,
      TRUE ~ NA_real_
    ),
    
    scr_kappa = LBXSCR / kappa,
    
    eGFR = 142 *
      pmin(scr_kappa, 1)^alpha *
      pmax(scr_kappa, 1)^(-1.200) *
      (0.9938^RIDAGEYR) *
      sex_factor
  )


# ------------------------------------------------------------
# 8. Construct operational CKD classification
#
# CKD is classified as present when:
#
#   eGFR < 60 mL/min/1.73 m^2
#            OR
#   UACR >= 30 mg/g
#
# IMPORTANT:
# NHANES measurements are cross-sectional. Therefore, this
# variable identifies participants meeting laboratory criteria
# consistent with CKD at the examination; it does not establish
# persistence for >3 months required for clinical diagnosis.
#
# PS: CKD is defined as abnormalities of
# kidney structure or function, present
# for >3 months, with implications for
# health according to https://kdigo.org/wp-content/uploads/2024/07/07232024-KDIGO-CKD.pdf
#
# Also, Someone with eGFR < 60 is classified positive even if ACR 
# is missing, and likewise someone with UACR ≥ 30 is positive even 
# if eGFR is missing. But we do not classify someone CKD-negative 
# unless both criteria are observed and negative. ( I need to verify 
# this from professionals in CKD even tho According to the  ⁠KDIGO 2024 CKD guideline,
# reduced eGFR or elevated ACR can independently indicate kidney abnormalities. )
# ------------------------------------------------------------

ckd_data <- ckd_data %>%
  mutate(
    low_egfr = case_when(
      is.na(eGFR) ~ NA_integer_,
      eGFR < 60 ~ 1L,
      TRUE ~ 0L
    ),
    
    albuminuria = case_when(
      is.na(URDACT) ~ NA_integer_,
      URDACT >= 30 ~ 1L,
      TRUE ~ 0L
    ),
    
    CKD = case_when(
      low_egfr == 1 | albuminuria == 1 ~ 1L,
      
      low_egfr == 0 & albuminuria == 0 ~ 0L,
      
      TRUE ~ NA_integer_
    )
  )


# ------------------------------------------------------------
# 9. Add readable categorical variables
# ------------------------------------------------------------

ckd_data <- ckd_data %>%
  mutate(
    sex = factor(
      RIAGENDR,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    ckd_status = factor(
      CKD,
      levels = c(0, 1),
      labels = c("No CKD", "CKD")
    )
  )

# ------------------------------------------------------------
# 9A. Define eligibility for primary analysis
# For the first empirical analysis, I want to restrict the 
# primary analysis to participants with both kidney measurements available
# ------------------------------------------------------------

ckd_data <- ckd_data %>%
  mutate(
    complete_kidney_markers =
      !is.na(eGFR) & !is.na(URDACT),
    
    primary_analysis_eligible =
      complete_kidney_markers &
      !is.na(LBXSUA)
  )

stopifnot(
  sum(ckd_data$primary_analysis_eligible) == 5343
)

# ------------------------------------------------------------
# 10. Remove temporary eGFR calculation variables
# ------------------------------------------------------------

ckd_data <- ckd_data %>%
  select(
    -kappa,
    -alpha,
    -sex_factor,
    -scr_kappa
  )


# ------------------------------------------------------------
# 11. Basic validation
# ------------------------------------------------------------

cat("Validation summary\n")
cat("------------------\n")

cat("Adult analytic frame:", nrow(ckd_data), "\n")

cat(
  "Serum creatinine available:",
  sum(!is.na(ckd_data$LBXSCR)),
  "\n"
)

cat(
  "UACR available:",
  sum(!is.na(ckd_data$URDACT)),
  "\n"
)

cat(
  "Serum uric acid available:",
  sum(!is.na(ckd_data$LBXSUA)),
  "\n"
)

core_complete <- ckd_data %>%
  filter(
    !is.na(LBXSCR),
    !is.na(URDACT),
    !is.na(LBXSUA)
  )

cat(
  "Complete creatinine + UACR + uric acid:",
  nrow(core_complete),
  "\n"
)

cat(
  "CKD status available:",
  sum(!is.na(ckd_data$CKD)),
  "\n\n"
)


# ------------------------------------------------------------
# 12. Sanity checks
# ------------------------------------------------------------

stopifnot(nrow(ckd_data) == 7809)

stopifnot(
  nrow(core_complete) == 5343
)

stopifnot(
  all(
    ckd_data$CKD[!is.na(ckd_data$CKD)] %in% c(0, 1)
  )
)

stopifnot(
  all(
    ckd_data$eGFR[!is.na(ckd_data$eGFR)] > 0
  )
)


# ------------------------------------------------------------
# 13. Display preliminary unweighted summaries
#
# These are DATA CHECKS only.
# They must not be interpreted as population prevalence
# estimates. Survey-weighted estimates will be produced in
# 02_exploratory_analysis.R.
# ------------------------------------------------------------

cat("Unweighted CKD classification among participants with status:\n")

print(
  ckd_data %>%
    filter(!is.na(CKD)) %>%
    count(CKD) %>%
    mutate(percent = 100 * n / sum(n))
)

cat("\nUnweighted summary of eGFR:\n")
print(summary(ckd_data$eGFR))

cat("\nUnweighted summary of serum uric acid:\n")
print(summary(ckd_data$LBXSUA))


# ------------------------------------------------------------
# 14. Save processed dataset
# ------------------------------------------------------------

saveRDS(
  ckd_data,
  output_path
)

cat(
  "\nProcessed dataset saved to:",
  output_path,
  "\n"
)


# ------------------------------------------------------------
# End of 01_data_preparation.R
# ------------------------------------------------------------
