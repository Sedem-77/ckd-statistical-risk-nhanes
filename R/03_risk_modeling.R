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
library(mgcv)
library(ggplot2)


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


# ------------------------------------------------------------
# 13. Weighted generalized additive model
#
# Purpose:
#   Examine whether the association between serum uric acid
#   and laboratory-defined CKD materially departs from a
#   linear functional form.
#
# Important:
#   The GAM uses NHANES phlebotomy sampling weights but does
#   not incorporate the complete complex survey variance
#   structure in the same way as svyglm().
#
#   Therefore:
#     - survey-weighted logistic regression remains the
#       primary inferential model;
#     - the GAM is used primarily to examine functional form.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Normalize sampling weights for GAM functional-form analysis
#
# mgcv treats supplied weights as prior weights. Raw NHANES
# sampling weights are therefore rescaled to have mean 1.
# This preserves their relative contribution while avoiding
# artificial inflation of the effective sample size.
# ------------------------------------------------------------

model_data <- model_data %>%
  mutate(
    gam_weight =
      WTPH2YR / mean(WTPH2YR, na.rm = TRUE)
  )

cat(
  "Sum of normalized GAM weights:",
  sum(model_data$gam_weight),
  "\n"
)

cat(
  "Mean normalized GAM weight:",
  mean(model_data$gam_weight),
  "\n"
)

gam_model <- gam(
  CKD ~
    s(LBXSUA, k = 5) +
    RIDAGEYR +
    sex +
    race_ethnicity,
  family = binomial(link = "logit"),
  weights = WTPH2YR,
  data = model_data,
  method = "REML"
)

cat("\n---------------- WEIGHTED GAM ----------------\n")
print(summary(gam_model))


# ------------------------------------------------------------
# 14. Inspect smooth-term complexity
#
# An estimated degrees of freedom (EDF) close to 1 suggests
# that the uric-acid relationship is approximately linear.
# Larger EDF indicates increasing evidence of curvature.
# ------------------------------------------------------------

gam_summary <- summary(gam_model)

smooth_results <- tibble(
  Term = rownames(gam_summary$s.table),
  EDF = gam_summary$s.table[, "edf"],
  Ref_df = gam_summary$s.table[, "Ref.df"],
  Chi_square = gam_summary$s.table[, "Chi.sq"],
  p_value = gam_summary$s.table[, "p-value"]
)

print(smooth_results)

write_csv(
  smooth_results,
  file.path(
    table_dir,
    "gam_smooth_results.csv"
  )
)


# ------------------------------------------------------------
# 15. Create prediction grid for serum uric acid
#
# Hold adjustment variables at representative values:
#   age = weighted mean age
#   sex = Female
#   race/ethnicity = Non-Hispanic White
#
# This curve is intended to visualize functional form.
# It is not a population prevalence curve.
# ------------------------------------------------------------

weighted_mean_age <- as.numeric(
  coef(
    svymean(
      ~RIDAGEYR,
      design = model_design,
      na.rm = TRUE
    )
  )
)

ua_limits <- quantile(
  model_data$LBXSUA,
  probs = c(0.01, 0.99),
  na.rm = TRUE
)

prediction_grid <- tibble(
  LBXSUA = seq(
    from = ua_limits[1],
    to = ua_limits[2],
    length.out = 200
  ),
  
  RIDAGEYR = weighted_mean_age,
  
  sex = factor(
    "Female",
    levels = levels(model_data$sex)
  ),
  
  race_ethnicity = factor(
    "Non-Hispanic White",
    levels = levels(model_data$race_ethnicity)
  )
)


# ------------------------------------------------------------
# 16. GAM predictions and confidence intervals
# ------------------------------------------------------------

gam_prediction <- predict(
  gam_model,
  newdata = prediction_grid,
  type = "link",
  se.fit = TRUE
)

prediction_grid <- prediction_grid %>%
  mutate(
    fit_link = as.numeric(gam_prediction$fit),
    se_link = as.numeric(gam_prediction$se.fit),
    
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link,
    
    predicted_probability = plogis(fit_link),
    lower_probability = plogis(lower_link),
    upper_probability = plogis(upper_link)
  )


# ------------------------------------------------------------
# 17. Adjusted uric-acid relationship figure
# ------------------------------------------------------------

ua_curve <- ggplot(
  prediction_grid,
  aes(
    x = LBXSUA,
    y = predicted_probability
  )
) +
  geom_ribbon(
    aes(
      ymin = lower_probability,
      ymax = upper_probability
    ),
    alpha = 0.20
  ) +
  geom_line(
    linewidth = 1
  ) +
  labs(
    title = "Adjusted Relationship Between Serum Uric Acid and Laboratory-Defined CKD",
    subtitle = paste0(
      "Weighted GAM; age held at ",
      round(weighted_mean_age, 1),
      " years, female, non-Hispanic White"
    ),
    x = "Serum uric acid (mg/dL)",
    y = "Predicted probability of laboratory-defined CKD",
    caption = paste(
      "NHANES 2021-2023.",
      "Outcome based on single-visit eGFR and UACR criteria;",
      "curve is for functional-form assessment and does not establish causality."
    )
  ) +
  theme_minimal(base_size = 12)

print(ua_curve)

ggsave(
  filename = file.path(
    figure_dir,
    "adjusted_uric_acid_ckd_curve.png"
  ),
  plot = ua_curve,
  width = 8,
  height = 5.5,
  dpi = 300
)


# ------------------------------------------------------------
# 18. Save prediction values
# ------------------------------------------------------------

write_csv(
  prediction_grid,
  file.path(
    table_dir,
    "adjusted_uric_acid_ckd_predictions.csv"
  )
)


# ------------------------------------------------------------
# 19. Simple model-performance assessment
#
# These are in-sample descriptive performance measures.
# They should NOT be described as externally validated
# predictive performance.
# ------------------------------------------------------------

model_data <- model_data %>%
  mutate(
    predicted_logistic = predict(
      model_2,
      type = "response"
    ),
    
    predicted_gam = predict(
      gam_model,
      type = "response"
    )
  )


# Brier score:
# mean squared difference between observed outcome and
# predicted probability.

weighted_brier <- function(y, p, w) {
  
  sum(
    w * (y - p)^2,
    na.rm = TRUE
  ) /
    sum(
      w,
      na.rm = TRUE
    )
}


brier_logistic <- weighted_brier(
  model_data$CKD,
  model_data$predicted_logistic,
  model_data$WTPH2YR
)

brier_gam <- weighted_brier(
  model_data$CKD,
  model_data$predicted_gam,
  model_data$WTPH2YR
)


performance_results <- tibble(
  Model = c(
    "Survey-weighted logistic model",
    "Weighted GAM"
  ),
  
  Brier_score = c(
    brier_logistic,
    brier_gam
  )
)

print(performance_results)

write_csv(
  performance_results,
  file.path(
    table_dir,
    "model_performance_results.csv"
  )
)


# ------------------------------------------------------------
# 20. Save GAM model
# ------------------------------------------------------------

saveRDS(
  gam_model,
  "data/processed/weighted_gam_model.rds"
)


cat("\nNonlinear modeling and preliminary model evaluation completed.\n")
