# CKD Statistical Risk Assessment Using NHANES

This repository contains an ongoing research project examining chronic kidney disease (CKD) identification, risk factors, and risk assessment in U.S. adults using data from the National Health and Nutrition Examination Survey (NHANES).

The project emphasizes transparent and reproducible statistical analysis of CKD risk patterns, with particular attention to clinically relevant continuous risk factors, population heterogeneity, nonlinear relationships, and the evaluation of risk models beyond discrimination alone.

## Research Objectives

The project currently focuses on three interconnected objectives:

1. **CKD Risk Prediction Literature and Gap Assessment**  
   Review representative statistical and machine learning approaches used for CKD identification and risk prediction, with particular attention to external validation, calibration, generalizability, interpretability, and reproducibility.

2. **CKD Data and Variable Inventory**  
   Develop a structured inventory of demographic, clinical, laboratory, behavioral, socioeconomic, and other variables relevant to CKD research and identify their availability in NHANES.

3. **NHANES CKD Risk Analysis**  
   Conduct a reproducible empirical analysis of CKD among U.S. adults using NHANES, beginning with the relationship between serum uric acid and CKD status and comparing conventional and nonlinear statistical models.

## Initial Empirical Question

The first empirical analysis investigates:

> Among U.S. adults represented in NHANES 2021–2023, how is serum uric acid associated with prevalent CKD after accounting for key demographic factors, does the association show evidence of nonlinearity, and how do transparent logistic and generalized additive models compare in discrimination and calibration?

The analysis will distinguish cross-sectional CKD classification from longitudinal CKD progression and will account for the NHANES complex survey design when making population-level estimates.

## Repository Structure

```text
literature/       Targeted literature review and research gap assessment
data_inventory/   CKD variable and data-source inventory
R/                Reproducible data preparation and statistical analysis
output/figures/   Analysis figures
output/tables/    Analysis tables
report/           Research reports and summaries
```

## Current Status

This is an active research project. The current phase focuses on establishing the CKD literature landscape, developing the data inventory, and constructing the initial NHANES analysis.

## Data

The project uses publicly available NHANES data from the U.S. Centers for Disease Control and Prevention (CDC), National Center for Health Statistics.

Raw NHANES data files are not stored in this repository. Data-source documentation and instructions required to reproduce the analytical dataset will be provided as the project develops.

## Statistical Analysis

Planned analyses include descriptive and survey-weighted estimation, logistic regression, generalized additive models, assessment of nonlinear risk relationships, and evaluation of model discrimination and calibration.

Additional analyses will be introduced when scientifically justified by the research questions and available data.

## Reproducibility

Analysis scripts, variable definitions, data-processing steps, and research outputs are maintained in this repository to support transparent and reproducible research.

## Author

Denis Folitse
Ph.D. Candidate in Statistics
University of Massachusetts Amherst
