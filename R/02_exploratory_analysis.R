
# ============================================================
# CKD Statistical Risk Assessment Using NHANES
# exploratory_analysis.R
#
# Purpose:
#   Conduct exploratory analysis of laboratory-
#   defined CKD among U.S. adults in NHANES 2021-2023.
#
# Input:
#   data/processed/ckd_analysis_2021_2023.rds
#
# Outputs:
#   output/tables/
#   output/figures/
#
# Important:
#   CKD in this project is an operational, laboratory-defined
#   classification based on single-visit NHANES measurements.
#   It does not establish the KDIGO chronicity requirement.
#
# Author: Denis Folitse
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# if (!requireNamespace("survey", quietly = TRUE)) {
#   install.packages("survey")
# }
library(survey)

# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

data_path <- "data/processed/ckd_analysis_2021_2023.rds"

table_dir <- "output/tables"
figure_dir <- "output/figures"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Read processed dataset
# ------------------------------------------------------------

ckd_data <- readRDS(data_path)

cat("Processed adult dataset:", nrow(ckd_data), "\n")


# ------------------------------------------------------------
# 4. Verify required variables
# ------------------------------------------------------------

required_vars <- c(
  "SEQN",
  "RIDAGEYR",
  "RIAGENDR",
  "RIDRETH3",
  "WTMEC2YR",
  "WTPH2YR",
  "SDMVSTRA",
  "SDMVPSU",
  "LBXSCR",
  "LBXSUA",
  "URDACT",
  "eGFR",
  "low_egfr",
  "albuminuria",
  "CKD",
  "primary_analysis_eligible"
)

missing_vars <- setdiff(required_vars, names(ckd_data))

if (length(missing_vars) > 0) {
  stop(
    "Required variables missing: ",
    paste(missing_vars, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 5. Sample-flow table
# ------------------------------------------------------------

sample_flow <- tibble(
  Stage = c(
    "NHANES 2021-2023 participants",
    "Adults age >= 20 years",
    "Serum creatinine available",
    "Urine ACR available",
    "Serum uric acid available",
    "Both kidney markers available",
    "Primary analysis eligible"
  ),
  
  N = c(
    11933,
    nrow(ckd_data),
    sum(!is.na(ckd_data$LBXSCR)),
    sum(!is.na(ckd_data$URDACT)),
    sum(!is.na(ckd_data$LBXSUA)),
    sum(!is.na(ckd_data$eGFR) &
          !is.na(ckd_data$URDACT)),
    sum(ckd_data$primary_analysis_eligible)
  )
)

print(sample_flow)

write_csv(
  sample_flow,
  file.path(table_dir, "sample_flow.csv")
)


# ------------------------------------------------------------
# 6. Missingness table
#
# Missingness is calculated relative to the full adult
# demographic frame (age >= 20), not only participants
# eligible for the laboratory analysis.
# ------------------------------------------------------------

variables_for_missingness <- c(
  "RIDAGEYR",
  "RIAGENDR",
  "RIDRETH3",
  "DMDEDUC2",
  "INDFMPIR",
  "LBXSCR",
  "LBXSUA",
  "URDACT",
  "eGFR",
  "CKD",
  "WTMEC2YR",
  "WTPH2YR",
  "SDMVSTRA",
  "SDMVPSU"
)

missingness <- tibble(
  Variable = variables_for_missingness,
  
  N_total = nrow(ckd_data),
  
  N_nonmissing = sapply(
    ckd_data[variables_for_missingness],
    function(x) sum(!is.na(x))
  ),
  
  N_missing = sapply(
    ckd_data[variables_for_missingness],
    function(x) sum(is.na(x))
  )
) %>%
  
  mutate(
    Percent_missing =
      round(100 * N_missing / N_total, 2)
  )

print(missingness)

write_csv(
  missingness,
  file.path(table_dir, "missingness_summary.csv")
)


# ------------------------------------------------------------
# 7. Define primary analytic sample
#
# Primary analysis requires:
#   - eGFR available
#   - UACR available
#   - serum uric acid available
#
# Additional model-specific missingness will be handled and
# reported when adjustment variables are introduced.
# ------------------------------------------------------------

analysis_data <- ckd_data %>%
  filter(primary_analysis_eligible)

cat(
  "\nPrimary analysis sample:",
  nrow(analysis_data),
  "\n"
)

stopifnot(nrow(analysis_data) == 5343)


# ------------------------------------------------------------
# 8. Validate survey-design information
# ------------------------------------------------------------

cat(
  "Missing phlebotomy weights:",
  sum(is.na(analysis_data$WTPH2YR)),
  "\n"
)

cat(
  "Non-positive phlebotomy weights:",
  sum(
    analysis_data$WTPH2YR <= 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Missing strata:",
  sum(is.na(analysis_data$SDMVSTRA)),
  "\n"
)

cat(
  "Missing PSU:",
  sum(is.na(analysis_data$SDMVPSU)),
  "\n"
)

psu_structure <- analysis_data %>%
  distinct(SDMVSTRA, SDMVPSU) %>%
  count(SDMVSTRA, name = "Number_of_PSUs")

print(psu_structure)

cat(
  "Number of strata with only one observed PSU:",
  sum(psu_structure$Number_of_PSUs == 1),
  "\n"
)

# ------------------------------------------------------------
# 9. Construct survey design
#
# WTPH2YR:
#   2021-2023 phlebotomy weight retained because the primary
#   analysis requires serum analytes, including serum uric acid
#   and serum creatinine.
#
# SDMVSTRA:
#   masked variance stratum.
#
# SDMVPSU:
#   masked primary sampling unit.
#
# nest = TRUE:
#   identifies Primary Sampling Units (PSUs) as nested within strata.
# ------------------------------------------------------------

options(survey.lonely.psu = "adjust") # we don't need this since the number of strata with only one PSU is 0

nhanes_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTPH2YR, 
  nest = TRUE,
  data = analysis_data
)


# ------------------------------------------------------------
# 10. Unweighted outcome frequencies
# ------------------------------------------------------------

unweighted_outcomes <- analysis_data %>%
  summarise(
    N = n(),
    
    CKD_n = sum(CKD == 1),
    CKD_percent = 100 * mean(CKD == 1),
    
    Low_eGFR_n = sum(low_egfr == 1),
    Low_eGFR_percent = 100 * mean(low_egfr == 1),
    
    Albuminuria_n = sum(albuminuria == 1),
    Albuminuria_percent =
      100 * mean(albuminuria == 1)
  )

print(unweighted_outcomes)

write_csv(
  unweighted_outcomes,
  file.path(
    table_dir,
    "unweighted_outcome_summary.csv"
  )
)

# ------------------------------------------------------------
# 10A. Decompose laboratory-defined CKD classification
#
# Participants in the primary analytic sample have both eGFR
# and UACR available. This table shows whether classification
# arises from reduced eGFR, albuminuria, both, or neither.
# ------------------------------------------------------------

ckd_component_overlap <- analysis_data %>%
  mutate(
    CKD_component = case_when(
      low_egfr == 0 & albuminuria == 0 ~ "Neither criterion",
      low_egfr == 1 & albuminuria == 0 ~ "Reduced eGFR only",
      low_egfr == 0 & albuminuria == 1 ~ "Albuminuria only",
      low_egfr == 1 & albuminuria == 1 ~ "Both criteria"
    )
  ) %>%
  count(CKD_component) %>%
  mutate(
    Percent = 100 * n / sum(n)
  )

print(ckd_component_overlap)

write_csv(
  ckd_component_overlap,
  file.path(
    table_dir,
    "ckd_component_overlap.csv"
  )
)
# ------------------------------------------------------------
# 11. Survey-weighted prevalence
#
# ------------------------------------------------------------

estimate_prevalence <- function(formula, design, outcome_name) {
  
  est <- svyciprop(
    formula,
    design,
    method = "logit",
    na.rm = TRUE
  )
  
  ci <- confint(est)
  
  tibble(
    Outcome = outcome_name,
    Weighted_prevalence =
      100 * as.numeric(coef(est)),
    SE =
      100 * as.numeric(SE(est)),
    CI_lower =
      100 * ci[1],
    CI_upper =
      100 * ci[2]
  )
}


weighted_prevalence <- bind_rows(
  
  estimate_prevalence(
    ~I(CKD == 1),
    nhanes_design,
    "Laboratory-defined CKD"
  ),
  
  estimate_prevalence(
    ~I(low_egfr == 1),
    nhanes_design,
    "Reduced eGFR (<60)"
  ),
  
  estimate_prevalence(
    ~I(albuminuria == 1),
    nhanes_design,
    "Albuminuria (UACR >=30 mg/g)"
  )
)

print(weighted_prevalence)

write_csv(
  weighted_prevalence,
  file.path(
    table_dir,
    "survey_weighted_prevalence.csv"
  )
)

# ------------------------------------------------------------
# 12. Preliminary serum uric acid summary
# ------------------------------------------------------------

uric_acid_summary <- analysis_data %>%
  summarise(
    N = n(),
    Mean = mean(LBXSUA),
    SD = sd(LBXSUA),
    Median = median(LBXSUA),
    Q1 = quantile(LBXSUA, 0.25),
    Q3 = quantile(LBXSUA, 0.75),
    Min = min(LBXSUA),
    Max = max(LBXSUA)
  )

print(uric_acid_summary)

write_csv(
  uric_acid_summary,
  file.path(
    table_dir,
    "uric_acid_summary.csv"
  )
)


# ------------------------------------------------------------
# 13. Serum uric acid by laboratory-defined CKD
# ------------------------------------------------------------

uric_acid_by_ckd <- analysis_data %>%
  group_by(ckd_status) %>%
  summarise(
    N = n(),
    Mean = mean(LBXSUA),
    SD = sd(LBXSUA),
    Median = median(LBXSUA),
    Q1 = quantile(LBXSUA, 0.25),
    Q3 = quantile(LBXSUA, 0.75),
    .groups = "drop"
  )

print(uric_acid_by_ckd)

write_csv(
  uric_acid_by_ckd,
  file.path(
    table_dir,
    "uric_acid_by_ckd.csv"
  )
)


# ------------------------------------------------------------
# 14. Figure 1: Distribution of serum uric acid
# ------------------------------------------------------------

p1 <- ggplot(
  analysis_data,
  aes(x = LBXSUA)
) +
  geom_histogram(
    bins = 35,
    boundary = 0
  ) +
  labs(
    title = "Distribution of Serum Uric Acid",
    subtitle = "NHANES 2021-2023 primary analytic sample",
    x = "Serum uric acid (mg/dL)",
    y = "Number of participants"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    figure_dir,
    "serum_uric_acid_distribution.png"
  ),
  p1,
  width = 8,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 15. Figure 2: Serum uric acid by CKD classification
# ------------------------------------------------------------

p2 <- ggplot(
  analysis_data,
  aes(
    x = ckd_status,
    y = LBXSUA
  )
) +
  geom_boxplot(
    outlier.alpha = 0.25
  ) +
  labs(
    title =
      "Serum Uric Acid by Laboratory-Defined CKD Status",
    subtitle =
      "Single-visit laboratory classification; unweighted display",
    x = NULL,
    y = "Serum uric acid (mg/dL)"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    figure_dir,
    "serum_uric_acid_by_ckd.png"
  ),
  p2,
  width = 7,
  height = 5,
  dpi = 300
)


# ------------------------------------------------------------
# 16. Uric-acid quintiles
#
# Quintiles are used here only as an exploratory descriptive
# device. Continuous modeling will be preferred in the formal
# regression analysis.
# ------------------------------------------------------------

analysis_data <- analysis_data %>%
  mutate(
    uric_acid_quintile =
      ntile(LBXSUA, 5)
  )


quintile_summary <- analysis_data %>%
  group_by(uric_acid_quintile) %>%
  summarise(
    N = n(),
    UA_min = min(LBXSUA),
    UA_max = max(LBXSUA),
    UA_median = median(LBXSUA),
    CKD_n = sum(CKD == 1),
    CKD_percent_unweighted =
      100 * mean(CKD == 1),
    .groups = "drop"
  )

print(quintile_summary)

write_csv(
  quintile_summary,
  file.path(
    table_dir,
    "uric_acid_quintile_summary.csv"
  )
)


# ------------------------------------------------------------
# 17. Save analysis sample
#
# This is a derived local dataset and remains excluded from Git.
# ------------------------------------------------------------

saveRDS(
  analysis_data,
  "data/processed/ckd_primary_analysis.rds"
)


cat("\nExploratory analysis completed successfully.\n")
cat("Tables saved to:", table_dir, "\n")
cat("Figures saved to:", figure_dir, "\n")



# ============================================================
# 18. Descriptive characteristics of primary analytic sample
# ============================================================


# ------------------------------------------------------------
# 18A. Create analysis-friendly demographic variables
# ------------------------------------------------------------

analysis_data <- analysis_data %>%
  mutate(
    
    age_group = cut(
      RIDAGEYR,
      breaks = c(20, 40, 60, Inf),
      right = FALSE,
      labels = c("20-39", "40-59", "60+")
    ),
    
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
    ),
    
    ckd_group = factor(
      CKD,
      levels = c(0, 1),
      labels = c(
        "No laboratory-defined CKD",
        "Laboratory-defined CKD"
      )
    )
  )


# Reconstruct survey design after adding variables
nhanes_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTPH2YR,
  nest = TRUE,
  data = analysis_data
)


# ------------------------------------------------------------
# 18B. Unweighted sample counts by CKD status
# ------------------------------------------------------------

table1_counts <- analysis_data %>%
  count(ckd_group, name = "Unweighted_N")

print(table1_counts)

write_csv(
  table1_counts,
  file.path(
    table_dir,
    "table1_unweighted_counts.csv"
  )
)


# ------------------------------------------------------------
# 18C. Survey-weighted continuous characteristics
# ------------------------------------------------------------

continuous_vars <- c(
  "RIDAGEYR",
  "LBXSUA",
  "eGFR",
  "URDACT"
)

continuous_labels <- c(
  "Age (years)",
  "Serum uric acid (mg/dL)",
  "eGFR (mL/min/1.73 m^2)",
  "UACR (mg/g)"
)


continuous_results <- bind_rows(
  lapply(
    seq_along(continuous_vars),
    function(j) {
      
      variable <- continuous_vars[j]
      
      f <- as.formula(
        paste0("~", variable)
      )
      
      result <- svyby(
        f,
        ~ckd_group,
        nhanes_design,
        svymean,
        na.rm = TRUE,
        vartype = "se"
      )
      
      tibble(
        Characteristic = continuous_labels[j],
        CKD_group = result$ckd_group,
        Weighted_mean = result[[variable]],
        SE = result[[paste0("se.", variable)]]
      )
    }
  )
)

print(continuous_results)

write_csv(
  continuous_results,
  file.path(
    table_dir,
    "table1_continuous_characteristics.csv"
  )
)


# ------------------------------------------------------------
# 18D. Survey-weighted categorical characteristics
# ------------------------------------------------------------

weighted_category_table <- function(
    design,
    variable,
    label
) {
  
  f <- as.formula(
    paste0("~", variable)
  )
  
  result <- svyby(
    f,
    ~ckd_group,
    design,
    svymean,
    na.rm = TRUE
  )
  
  long <- result %>%
    select(
      ckd_group,
      starts_with(variable)
    ) %>%
    tidyr::pivot_longer(
      cols = -ckd_group,
      names_to = "Level",
      values_to = "Weighted_proportion"
    ) %>%
    mutate(
      Characteristic = label,
      Level = sub(
        paste0("^", variable),
        "",
        Level
      ),
      Weighted_percent =
        100 * Weighted_proportion
    ) %>%
    select(
      Characteristic,
      Level,
      ckd_group,
      Weighted_percent
    )
  
  long
}


categorical_results <- bind_rows(
  
  weighted_category_table(
    nhanes_design,
    "sex",
    "Sex"
  ),
  
  weighted_category_table(
    nhanes_design,
    "age_group",
    "Age group"
  ),
  
  weighted_category_table(
    nhanes_design,
    "race_ethnicity",
    "Race/ethnicity"
  )
)

print(categorical_results)

write_csv(
  categorical_results,
  file.path(
    table_dir,
    "table1_categorical_characteristics.csv"
  )
)


# ============================================================
# 19. Compare adult frame with primary analytic sample
# ============================================================


# ------------------------------------------------------------
# 19A. Create inclusion indicator and demographic variables
# ------------------------------------------------------------

sample_comparison <- ckd_data %>%
  mutate(
    
    included_primary =
      factor(
        if_else(
          primary_analysis_eligible,
          1L,
          0L
        ),
        levels = c(0, 1),
        labels = c(
          "Not included",
          "Primary analytic sample"
        )
      ),
    
    sex = factor(
      RIAGENDR,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    
    age_group = cut(
      RIDAGEYR,
      breaks = c(20, 40, 60, Inf),
      right = FALSE,
      labels = c(
        "20-39",
        "40-59",
        "60+"
      )
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
# 19B. Unweighted inclusion counts
# ------------------------------------------------------------

sample_inclusion_counts <- sample_comparison %>%
  count(included_primary) %>%
  mutate(
    Percent =
      100 * n / sum(n)
  )

print(sample_inclusion_counts)

write_csv(
  sample_inclusion_counts,
  file.path(
    table_dir,
    "analytic_sample_inclusion_counts.csv"
  )
)


# ------------------------------------------------------------
# 19C. Unweighted demographic comparison
#
# This comparison describes observed differences between adults
# retained versus not retained in the primary analytic sample.
# It is a diagnostic of analytic-sample construction rather than
# a population prevalence analysis.
# ------------------------------------------------------------

sample_demographic_comparison <- sample_comparison %>%
  group_by(included_primary) %>%
  summarise(
    
    N = n(),
    
    Mean_age =
      mean(RIDAGEYR, na.rm = TRUE),
    
    Female_percent =
      100 * mean(
        RIAGENDR == 2,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

print(sample_demographic_comparison)

write_csv(
  sample_demographic_comparison,
  file.path(
    table_dir,
    "analytic_sample_demographic_comparison.csv"
  )
)


# ------------------------------------------------------------
# 19D. Age-group distribution by inclusion status
# ------------------------------------------------------------

sample_age_distribution <- sample_comparison %>%
  count(
    included_primary,
    age_group
  ) %>%
  group_by(included_primary) %>%
  mutate(
    Percent =
      100 * n / sum(n)
  ) %>%
  ungroup()

print(sample_age_distribution)

write_csv(
  sample_age_distribution,
  file.path(
    table_dir,
    "analytic_sample_age_distribution.csv"
  )
)


# ------------------------------------------------------------
# 19E. Race/ethnicity distribution by inclusion status
# ------------------------------------------------------------

sample_race_distribution <- sample_comparison %>%
  count(
    included_primary,
    race_ethnicity
  ) %>%
  group_by(included_primary) %>%
  mutate(
    Percent =
      100 * n / sum(n)
  ) %>%
  ungroup()

print(sample_race_distribution)

write_csv(
  sample_race_distribution,
  file.path(
    table_dir,
    "analytic_sample_race_distribution.csv"
  )
)

# ------------------------------------------------------------
# End of 02_exploratory_analysis.R
# ------------------------------------------------------------

