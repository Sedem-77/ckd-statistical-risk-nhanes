# Targeted CKD Prediction Literature: Preliminary Gap Analysis

## Scope

This is a targeted, systematic-style landscape review rather than a
formal systematic review. Its purpose is to orient the CKD
risk-assessment project, identify representative model families and
datasets, and prevent unsupported novelty claims before the NHANES
analysis begins.

## What already exists

The literature does **not** support a claim that CKD lacks prediction
models. Traditional multivariable risk scores, logistic and Cox models,
the Kidney Failure Risk Equation (KFRE), machine-learning models, and
more recent interpretable ML approaches are already well represented.
CKD studies have addressed incident/prevalent CKD, progression to kidney
failure or renal replacement therapy, cardiovascular outcomes,
mortality, and disease-specific populations such as diabetes and
advanced CKD.

Recent NHANES work also means that **NHANES + uric acid + nonlinear
modeling or machine learning is not, by itself, a methodological novelty
claim**. Published studies have already used NHANES for CKD prediction
and have applied logistic regression, restricted cubic splines, SHAP,
CatBoost, and other ML methods.

## Recurring gaps supported by the literature

### 1. External validation and transportability

Models that perform well in their derivation samples do not necessarily
retain the same performance in new populations. External validation
remains limited in several reviews, and recent large validation work
continues to show population-dependent performance.

### 2. Calibration is often weaker or less emphasized than discrimination

Many CKD prediction papers emphasize AUC/C-statistics. Reviews and
external validations repeatedly show that a model can discriminate
reasonably well while substantially over- or under-predicting absolute
risk. Our project should therefore avoid presenting AUC alone as
sufficient evidence of model quality.

### 3. Population heterogeneity and subgroup performance

Generalizability across demographic and clinical subgroups remains an
issue. Performance may differ in people with diabetes, hypertension,
across ethnic groups, or across healthcare settings. Diverse-population
testing has been limited in parts of the ML literature.

### 4. Reproducibility and reporting

Some prognostic-model reviews report incomplete model equations,
heterogeneous performance metrics, and difficulty reproducing or
directly validating published models. A transparent public repository,
explicit variable definitions, reproducible code, and clear reporting
are therefore meaningful strengths for this project.

### 5. Clinical interpretability and usefulness

Machine learning is already common in CKD. The unresolved question is
not simply whether a more complex algorithm can increase apparent
accuracy. Clinical evaluation, interpretability, decision usefulness,
and validation in realistic settings remain less developed.

### 6. Continuous risk-factor relationships deserve careful treatment

Earlier reviews specifically noted dichotomization of continuous
predictors and untested linearity assumptions as modeling weaknesses.
Modern CKD studies increasingly use splines and interpretable nonlinear
tools, so nonlinearity itself is not novel. However, carefully
characterizing continuous risk relationships remains scientifically
useful, particularly when paired with transparent modeling and
appropriate validation.

## Implications for our first NHANES study

The first report should be framed as a **transparent, reproducible CKD
risk-assessment analysis**, not as a claim to have invented a new CKD
prediction method.

A defensible initial question is:

> Among U.S. adults represented in NHANES 2021--2023, how is serum uric
> acid associated with prevalent CKD after accounting for key
> demographic factors, does the association show evidence of
> nonlinearity, and how do transparent logistic and generalized additive
> models compare in discrimination and calibration?

The analysis should:

-   use the NHANES complex survey design where population-representative
    estimates are claimed;
-   clearly distinguish cross-sectional CKD classification from
    longitudinal CKD progression;
-   treat eGFR and UACR definitions transparently;
-   retain continuous predictors where scientifically appropriate rather
    than unnecessarily dichotomizing them;
-   compare a conventional logistic specification with an interpretable
    nonlinear model such as a GAM;
-   report calibration as well as discrimination;
-   document missingness and analytic exclusions;
-   avoid claiming that NHANES, uric acid, GAMs, splines, or machine
    learning are themselves novel.

## Working gap statement

> Although numerous statistical and machine-learning models for CKD
> identification and progression have been developed, recurring
> limitations remain in external validity, calibration, population
> transportability, reproducibility, and clinically interpretable
> characterization of risk relationships. This project begins by
> developing a transparent, survey-aware and reproducible assessment of
> CKD risk patterns in a contemporary U.S. population dataset, with
> explicit evaluation of continuous and potentially nonlinear risk
> relationships and model performance beyond discrimination alone.

This is a **working** gap statement. It should be refined as the
literature matrix expands and as the empirical analysis reveals what can
actually be supported.

## Selected sources

-   Tangri et al. (2013): https://pubmed.ncbi.nlm.nih.gov/23588748/
-   van Rijn et al. (2021): https://pubmed.ncbi.nlm.nih.gov/33051669/
-   Sanmarchi et al. (2023): https://pubmed.ncbi.nlm.nih.gov/36786976/
-   Echouffo-Tcheugui & Kengne (2012):
    https://pmc.ncbi.nlm.nih.gov/articles/PMC3502517/
-   2026 UK Biobank external validation:
    https://pmc.ncbi.nlm.nih.gov/articles/PMC13203005/
-   2026 NHANES elderly metabolic-syndrome study:
    https://pubmed.ncbi.nlm.nih.gov/41371294/
-   2024 NHANES abdominal-obesity study:
    https://pubmed.ncbi.nlm.nih.gov/39593076/
