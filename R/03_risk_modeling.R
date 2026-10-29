# ============================================================
# CKD Statistical Risk Assessment Using NHANES
# 03_risk_modeling.R
#
# Purpose:
#   Evaluate the association between serum uric acid and
#   single-visit laboratory-defined CKD among U.S. adults
#   represented in NHANES 2021-2023.
#
# Primary outcome:
#   Laboratory-defined CKD based on:
#     eGFR < 60 mL/min/1.73 m^2 OR
#     UACR >= 30 mg/g
#
# Important:
#   NHANES provides single-visit laboratory measurements.
#   Therefore, this outcome does not establish the chronicity
#   required for a clinical diagnosis of CKD.
#
# Author: Denis Folitse
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(readr)
library(survey)


# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

data_path <- "data/processed/ckd_primary_analysis.rds"

table_dir <- "output/tables"
figure_dir <- "output/figures"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Read primary analytic dataset
# ------------------------------------------------------------

analysis_data <- readRDS(data_path)

cat(
  "Primary analytic dataset:",
  nrow(analysis_data),
  "\n"
)

stopifnot(nrow(analysis_data) == 5343)


# ------------------------------------------------------------
# 4. Prepare modeling variables
# ------------------------------------------------------------

analysis_data <- analysis_data %>%
  mutate(
    
    sex = factor(
      RIAGENDR,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    race_ethnicity = factor(
      RIDRETH3,
      levels = c(1, 2, 3, 4, 6, 7),
      labels = c(
        "Mexican American",
        "Other Hispanic",
        "Non-Hispanic White",
        "Non-Hispanic Black",
        "Non-Hispanic Asian",
        "Other/Multiracial"
      )
    )
  )


# ------------------------------------------------------------
# 5. Check model-variable completeness
# ------------------------------------------------------------

model_variables <- c(
  "CKD",
  "LBXSUA",
  "RIDAGEYR",
  "sex",
  "race_ethnicity",
  "WTPH2YR",
  "SDMVSTRA",
  "SDMVPSU"
)

model_missingness <- tibble(
  Variable = model_variables,
  
  N_missing = sapply(
    analysis_data[model_variables],
    function(x) sum(is.na(x))
  )
)

print(model_missingness)


# ------------------------------------------------------------
# 6. Define common modeling sample
#
# Models are fitted to the same observations so that changes
# between crude and adjusted estimates reflect covariate
# adjustment rather than changes in the analytic sample.
# ------------------------------------------------------------

model_data <- analysis_data %>%
  filter(
    complete.cases(
      CKD,
      LBXSUA,
      RIDAGEYR,
      sex,
      race_ethnicity,
      WTPH2YR,
      SDMVSTRA,
      SDMVPSU
    )
  )

cat(
  "\nCommon modeling sample:",
  nrow(model_data),
  "\n"
)

cat(
  "Excluded from modeling because of missing model variables:",
  nrow(analysis_data) - nrow(model_data),
  "\n"
)


# ------------------------------------------------------------
# 7. Construct NHANES survey design
# ------------------------------------------------------------

options(survey.lonely.psu = "adjust")

model_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTPH2YR,
  nest = TRUE,
  data = model_data
)


# ------------------------------------------------------------
# 8. Model 1: Crude association
#
# logit[P(CKD = 1)] =
#   beta_0 + beta_1(serum uric acid)
#
# beta_1 represents the change in log odds associated with a
# 1 mg/dL higher serum uric acid concentration.
# ------------------------------------------------------------

model_1 <- svyglm(
  CKD ~ LBXSUA,
  design = model_design,
  family = quasibinomial()
)

cat("\n---------------- MODEL 1 ----------------\n")
print(summary(model_1))


# ------------------------------------------------------------
# 9. Model 2: Demographically adjusted association
#
# Adjusts for:
#   age (continuous)
#   sex
#   race/ethnicity
#
# Age is initially modeled linearly. Nonlinear functional forms
# will be considered subsequently.
# ------------------------------------------------------------

model_2 <- svyglm(
  CKD ~ LBXSUA +
    RIDAGEYR +
    sex +
    race_ethnicity,
  design = model_design,
  family = quasibinomial()
)

cat("\n---------------- MODEL 2 ----------------\n")
print(summary(model_2))


# ------------------------------------------------------------
# 10. Extract odds ratios and 95% confidence intervals
#
# Survey-based coefficient estimates are on the log-odds scale.
#
# OR = exp(beta)
#
# Confidence intervals are obtained on the coefficient scale
# and exponentiated to produce OR confidence intervals.
# ------------------------------------------------------------

extract_or <- function(model, model_name) {
  
  beta <- coef(model)
  se <- sqrt(diag(vcov(model)))
  
  model_summary <- summary(model)
  coef_table <- coef(model_summary)
  
  tibble(
    Model = model_name,
    Term = names(beta),
    Beta = as.numeric(beta),
    SE = as.numeric(se),
    
    OR = exp(Beta),
    
    CI_lower = exp(
      Beta - qt(
        0.975,
        df = model$df.residual
      ) * SE
    ),
    
    CI_upper = exp(
      Beta + qt(
        0.975,
        df = model$df.residual
      ) * SE
    ),
    
    p_value =
      coef_table[, "Pr(>|t|)"]
  )
}


model_1_results <- extract_or(
  model_1,
  "Model 1: Crude"
)

model_2_results <- extract_or(
  model_2,
  "Model 2: Demographically adjusted"
)


model_results <- bind_rows(
  model_1_results,
  model_2_results
)

print(model_results, n = Inf)


write_csv(
  model_results,
  file.path(
    table_dir,
    "survey_logistic_model_results.csv"
  )
)


# ------------------------------------------------------------
# 11. Extract serum uric acid estimates for direct comparison
# ------------------------------------------------------------

uric_acid_model_comparison <- model_results %>%
  filter(Term == "LBXSUA") %>%
  select(
    Model,
    Beta,
    SE,
    OR,
    CI_lower,
    CI_upper,
    p_value
  )

print(uric_acid_model_comparison)

write_csv(
  uric_acid_model_comparison,
  file.path(
    table_dir,
    "uric_acid_model_comparison.csv"
  )
)


# ------------------------------------------------------------
# 12. Save fitted models locally
#
# Model objects are derived research outputs and are not
# required in the public GitHub repository.
# ------------------------------------------------------------

saveRDS(
  list(
    model_1 = model_1,
    model_2 = model_2
  ),
  "data/processed/survey_logistic_models.rds"
)


cat("\nSurvey-weighted logistic models completed.\n")

