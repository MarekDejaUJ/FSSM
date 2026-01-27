# FSSM: Fuzzy Subjective Structure Model

Here is the full content for your `README.md` file. You can copy this code block and save it directly as `README.md` in your package root directory.

# FSSM: Fuzzy Subjective Structure Model

[![R Package](https://img.shields.io/badge/R%20Package-0.2.0-blue.svg)](https://github.com/marekdeja/fssm) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

## Overview

**FSSM** (Fuzzy Subjective Structure Model) is an R package implementing a hybrid causal procedure designed for information behaviour research and social sciences. It bridges the gap between **structural equation modeling (SEM)** and **necessary condition analysis (NCA)** by introducing fuzzy logic and propensity score weighting.

FSSM allows researchers to answer questions that standard regression misses: \* *"Is X necessary for Y, or just helpful?"* \* *"Does X trigger Y deterministically, or just increase its probability?"* \* *"How does the causal effect change as X increases (Dose-Response)?"*

### Key Features

-   **Dual Entry Routes:** Works with **PLS-SEM models** (via `seminr`) OR **Raw Data** (via entropy/MCDM weighting).
-   **Causal Balancing:** Implements **Generalized Propensity Scores (GPS)** for continuous treatments and **Binary IPW** for threshold-based analysis.
-   **Augmented cIPMA:** Extends the Importance-Performance Map with necessity bottlenecks and causal uplift.
-   **Fuzzy Logic:** Formalizes four distinct causal claims (Typical/Probabilistic Sufficiency & Necessity).

## The Four Causal Claims

FSSM moves beyond simple "significance" to test four specific logical structures:

| Claim | Estimand | Formula | Interpretation |
|-----------------|:---------------:|:---------------:|-----------------------|
| **Typical Sufficiency** | $\alpha$ | $\tilde{\alpha} \ge 0.85$ | *"If X is high, Y is **typically** high."* (Consistency) |
| **Typical Necessity** | $\epsilon$ | $\tilde{\epsilon} \le 0.05$ | *"If Y is high, X is **rarely** low."* (Exception Rate) |
| **Probabilistic Sufficiency** | $\Delta$ | $\tilde{\Delta} \ge 0.15$ | *"Increasing X **raises the probability** of high Y."* (Causal Uplift) |
| **Probabilistic Necessity** | $\beta$ | $\tilde{\beta} \le 0.20$ | *"If X is low, high Y is **unlikely**."* (Constraint Strength) |

------------------------------------------------------------------------

## Installation

``` r
# Install development version from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("MarekDejaUJ/FSSM")
```

------------------------------------------------------------------------

## Quick Start: The Two Routes

FSSM supports two workflows. Choose the one that fits your data availability.

### Route A: The PLS-SEM Route (Recommended)

Use this if you have a latent variable model defined in `seminr`. FSSM uses the model's outer weights and data.

``` r
library(fssm)
library(seminr)

# 1. Load data and define scales (Theoretical Min/Max)
data("tam_data")

# Define scales (theoretical min/max for each indicator)
my_scales <- list(
  PU_01 = c(1, 5), PU_02 = c(1, 5), PU_03 = c(1, 5),
  CO_01 = c(1, 5), CO_02 = c(1, 5), CO_03 = c(1, 5),
  EOU_01 = c(1, 5), EOU_02 = c(1, 5), EOU_03 = c(1, 5),
  EMV_01 = c(1, 5), EMV_02 = c(1, 5), EMV_03 = c(1, 5),
  AD_01 = c(1, 5), AD_02 = c(1, 5), AD_03 = c(1, 5),
  USE_01 = c(1, 7)
)

# 2. Estimate and Bootstrap PLS Model
tam_mm <- constructs(
  composite("Perceived_Usefulness", multi_items("PU_0", 1:3)),
  composite("Compatibility",        multi_items("CO_0", 1:3)),
  composite("Ease_of_Use",          multi_items("EOU_0", 1:3)),
  composite("Emotional_Value",      multi_items("EMV_0", 1:3)),
  composite("Adoption_Intention",   multi_items("AD_0", 1:3)),
  composite("Technology_Use",       single_item("USE_01"))
)

tam_sm <- relationships(
  paths(from = c("Compatibility", "Perceived_Usefulness", "Ease_of_Use", "Emotional_Value"),
        to = "Adoption_Intention"),
  paths(from = c("Adoption_Intention", "Compatibility", "Perceived_Usefulness",
                 "Ease_of_Use", "Emotional_Value"),
        to = "Technology_Use")
)

tam_pls <- estimate_pls(tam_data, tam_mm, tam_sm)
boot_tam <- bootstrap_model(tam_pls, nboot = 500)

# 3. Run FSSM with GPS Weighting
result_A <- fssm(
  input = boot_model,
  target_construct = "Technology_Use",
  data = tam_data,
  scales = my_scales,
  target_level = 85,
  confounders = c("Age", "Gender"), # Adjusts for confounding
  weighting = "gps"
)

# 4. View Results
print(result_A)
```

### Route B: The Raw Data / Entropy Route

Use this if you want to create composite scores directly from raw data using **Entropy Weighting** (Shannon's Entropy), which weights items based on their information content rather than loading strength.

``` r
# 1. Define syntax for raw data route (lavaan-style)
fssm_syntax <- "
  Perceived_Usefulness =~ PU_01 + PU_02 + PU_03;
  Compatibility =~ CO_01 + CO_02 + CO_03;
  Ease_of_Use =~ EOU_01 + EOU_02 + EOU_03;
  Emotional_Value =~ EMV_01 + EMV_02 + EMV_03;
  Adoption_Intention =~ AD_01 + AD_02 + AD_03;
  Technology_Use =~ USE_01
"

# 2. Extract scores using entropy weighting
fssm_obj_entropy <- fssm_extract(
  input = tam_data,
  syntax = fssm_syntax,
  scales = my_scales,  # Optional: use theoretical rescaling
  indicator_weighting = "entropy"
)

# 3. Run FSSM
result_B <- fssm(
  input = fssm_obj_entropy,
  target_construct = "Technology_Use",
  target_level = 85,
  confounders = c("Age", "Gender"),
  weighting = "gps"
)

print(result_B)
```

------------------------------------------------------------------------

## Visualization

FSSM includes a comprehensive plotting suite using `ggplot2`.

### 1. Summary Dashboard

The easiest way to see everything at once: IPMA, Heatmap, Dose-Response, and Bottlenecks.

``` r
plot_dashboard(result)
```

### 2. Augmented IPMA

Visualizes Importance vs. Performance, colored by Necessity status.

``` r
# Standard IPMA
plot(result, x_axis = "importance", y_axis = "performance", color_by = "necessity")

# Causal Strategy Map (Uplift vs. Consistency)
plot(result, x_axis = "delta", y_axis = "alpha", bubble_size = "performance")
```

### 3. Dose-Response Curves

Shows the causal "uplift": how the probability of the outcome changes as membership in X increases.

``` r
plot_dose_response(result)
```

------------------------------------------------------------------------

## Methodological Workflow

```         
┌─────────────────────────────────────────────────────────────────┐
│                       FSSM WORKFLOW                             │
├─────────────────────────────────────────────────────────────────┤
│  Route A: PLS-SEM Weights   │   Route B: Entropy/Equal Weights  │
│  (bootstrapped model)       │   (raw data extraction)           │
└───────────────┬─────────────┴──────────────┬────────────────────┘
                │                            │
                ▼                            ▼
┌─────────────────────────────────────────────────────────────────┐
│  Step 1: Rescaling & Standardization (0-100 Theoretical Range)  │
├─────────────────────────────────────────────────────────────────┤
│  Step 2: NCA Bottleneck Calculation (Ceiling Lines)             │
│          └─> Defines Thresholds T_X for fuzzy membership        │
├─────────────────────────────────────────────────────────────────┤
│  Step 3: Causal Balancing (Weighting)                           │
│          └─> GPS (Continuous) or Binary IPW based on Covariates │
├─────────────────────────────────────────────────────────────────┤
│  Step 4: Fuzzy Logic Causal Claims                              │
│          └─> Calculate α, ε, Δ, β using weighted aggregation    │
└───────────────────────────────┬─────────────────────────────────┘
                                │
                                ▼
                     ┌──────────────────────┐
                     │  Strategic Output    │
                     │  - Dashboard Plots   │
                     │  - Diagnostic Tables │
                     │  - Logic Checks      │
                     └──────────────────────┘
```

## Key Difference: GPS vs Binary IPW

| Aspect      | Binary IPW               | GPS (Recommended)          |
|-------------|--------------------------|----------------------------|
| Treatment   | Dichotomized:            | Continuous:                |
| Information | Lost (71 = 100 if both ) | Preserved (0.03 1.0)       |
| Output      | Two-point comparison     | Full dose-response curve   |
| Balance     | At threshold only        | Across entire distribution |

## Diagnostics

FSSM assumes that your measurement model and weighting model are valid. Check them easily:

``` r
# Check PLS assumptions (Reliability, HTMT, VIF) and Weighting (ESS)
check_assumptions(result)
```

## Citation

If you use FSSM in your research, please cite:

``` bibtex
@article{deja2026fssm,
  title={Fuzzy Subjective Structure Model: A Hybrid Causal Procedure for Information Behaviour Research},
  author={Deja, Marek},
  year={2026},
  note={R package version 0.2.0},
  url={[https://github.com/marekdeja/fssm](https://github.com/marekdeja/fssm)}
}
```

## Author

**Marek Deja** Assistant Professor, Jagiellonian University, Kraków, Poland

Email: [marek.deja\@uj.edu.pl](mailto:marek.deja@uj.edu.pl)

