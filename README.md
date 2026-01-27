# FSSM: Fuzzy Subjective Structure Model

[![R Package](https://img.shields.io/badge/R%20Package-0.1.0-blue.svg)](https://github.com/marekdeja/fssm)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**FSSM** (Fuzzy Subjective Structure Model) is an R package implementing a hybrid causal procedure for analyzing information behaviour in academic settings. It integrates:

- **PLS-SEM** sufficiency paths (what usually works)
- **NCA** necessity bottlenecks (what constrains action)
- **cIPMA** decision-oriented visualization
- **Entropy weighting** for subjective measurement
- **GPS/IPW weighting** for causal balance

The package formalizes four distinct, non-deterministic causal claims:

| Claim | Formula | Interpretation |
|-------|---------|----------------|
| **Typical Sufficiency** | α̃ ≥ 0.85 | "If X, then typically Y" |
| **Typical Necessity** | ε̃ ≤ 0.05 | "If not X, then typically not Y" |
| **Probabilistic Sufficiency** | Δ̃ ≥ 0.15 | "If X, then probably Y" (uplift) |
| **Probabilistic Necessity** | β̃ ≤ 0.20 | "If not X, then probably not Y" |

## Installation

```r
# Install from GitHub (when available)
# devtools::install_github("marekdeja/fssm")

# Or install from local source
install.packages("path/to/fssm", repos = NULL, type = "source")
```

## Dependencies

- `seminr` (>= 2.3.0) for PLS-SEM modeling
- `NCA` (>= 4.0.0) for Necessary Condition Analysis
- `ggplot2` (>= 3.4.0) for visualization
- `ggrepel` (>= 0.9.0) for label positioning

## Quick Start

```r
library(fssm)
library(seminr)

# 1. Run your PLS-SEM model with seminr (bootstrapped)
boot_model <- bootstrap_model(
  seminr_model = pls_model,
  nboot = 5000
)

# 2. Define theoretical scale bounds
scales <- list(
  "IL1" = c(1, 7),  # 7-point Likert
  "IC1" = c(1, 7),
  "IE1" = c(1, 7)
)

# 3. Run FSSM analysis
fssm_results <- fssm(
  model = boot_model,
  target_construct = "Empowerment",
  data = survey_data,
  scales = scales,
  target_level = 85,
  confounders = c("discipline", "gender", "seniority"),
  weighting = "gps"  # or "binary" for traditional IPW
)

# 4. View results
print(fssm_results)
summary(fssm_results)

# 5. Visualize
plot(fssm_results)
plot(fssm_results, claim = "alpha", y_axis = "beta")
plot_dose_response(fssm_results)
plot_claims_comparison(fssm_results)

# 6. Check diagnostics
check_assumptions.fssm(fssm_results)
```

## Key Functions

### Main Analysis
- `fssm()` - Run complete FSSM analysis
- `print.fssm()` - Display summary results
- `summary.fssm()` - Detailed statistics
- `extract_fssm()` - Extract specific components

### Membership Functions
- `calculate_fuzzy_membership()` - Compute μ_X^S, μ_X^N, μ_Y
- `summarize_membership()` - Membership distribution statistics
- `get_case_profiles()` - Case-level membership details

### Weighting
- `calculate_gps_weights()` - Generalized Propensity Score weights
- `calculate_binary_ipw_weights()` - Traditional binary IPW
- `check_gps_diagnostics()` - Weighting quality assessment

### Causal Claims
- `calculate_crisp_claims()` - Binary threshold formulations
- `calculate_fuzzy_claims()` - Continuous membership formulations
- `estimate_dose_response()` - Full sufficiency-outcome curve
- `bootstrap_claims()` - Confidence intervals via bootstrap
- `assess_claim_support()` - Decision support assessment

### Visualization
- `plot.fssm()` - Enhanced importance-performance map
- `plot_dose_response()` - Sufficiency-outcome relationship
- `plot_claims_comparison()` - Heatmap/bar of four claims
- `plot_membership()` - Membership distribution plots

### Diagnostics
- `check_assumptions.fssm()` - PLS-SEM and weighting diagnostics

## Methodological Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                     FSSM WORKFLOW                               │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Step 1: Entropy-weighted latent scores                         │
│          └─> X*_i, Y*_i ∈ [0, 100]                              │
│                                                                 │
│  Step 2: NCA bottleneck estimation (unweighted)                 │
│          └─> T_X = b_X(T_Y), bottleneck function                │
│                                                                 │
│  Step 3: Fuzzy membership (BEFORE weighting)                    │
│          └─> μ_X^S(i), μ_X^N(i), μ_Y(i) ∈ [0, 1]                │
│                                                                 │
│  Step 4: GPS weighting on continuous μ_X^S                      │
│          └─> Stabilized weights w_i                             │
│                                                                 │
│  Step 5: Dose-response estimation                               │
│          └─> E[μ_Y | μ_X^S = m] for m ∈ [0, 1]                  │
│                                                                 │
│  Step 6: Four causal claims (fuzzy formulations)                │
│          └─> α̃, ε̃, Δ̃, β̃ + violation mass + disc. power          │
│                                                                 │
│  Step 7: Crisp claims (optional comparison)                     │
│          └─> α̂, ε̂, Δ̂, β̂                                         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Key Difference: GPS vs Binary IPW

| Aspect | Binary IPW | GPS (Recommended) |
|--------|-----------|-------------------|
| Treatment | Dichotomized: T = 1{X ≥ T_X} | Continuous: μ_X^S ∈ [0,1] |
| Information | Lost (71 = 100 if both ≥ T_X) | Preserved (0.03 ≠ 1.0) |
| Output | Two-point comparison | Full dose-response curve |
| Balance | At threshold only | Across entire distribution |

## Interpretation Guide

### When to use which claim:

- **Typical Sufficiency (α̃)**: "Is the condition generally associated with success?"
  - High α̃: When X is present, Y typically occurs
  - Use for identifying enablers

- **Typical Necessity (ε̃)**: "Is the condition required for success?"
  - Low ε̃: High Y rarely occurs without sufficient X
  - Use for identifying gatekeepers/constraints

- **Probabilistic Sufficiency (Δ̃)**: "How much does X increase the probability of Y?"
  - High Δ̃: X meaningfully increases likelihood of Y
  - Use for effect size interpretation

- **Probabilistic Necessity (β̃)**: "What happens when X is absent?"
  - Low β̃: Without X, Y is unlikely
  - Use for risk assessment

### Construct Classification:

| α̃ | ε̃ | Classification |
|---|---|----------------|
| High | Low | **Strong Enabler & Gatekeeper** |
| High | High | Enabler (not necessary) |
| Low | Low | Gatekeeper (necessary but not sufficient) |
| Low | High | Weak or unclear relationship |

## Citation

If you use FSSM in your research, please cite:

```bibtex
@article{deja2026fssm,
  title={Fuzzy Subjective Structure Model: A Hybrid Causal Procedure for Information Behaviour Research},
  author={Deja, Marek},
  journal={Proceedings of ISIC 2026},
  year={2026}
}
```

## Further Development

The FSSM project is actively developed. Follow updates at:
- OSF: https://osf.io/b8maq/
- GitHub: https://github.com/marekdeja/fssm

## License

MIT License - see LICENSE file for details.

## Author

**Marek Deja**  
Associate Professor, Jagiellonian University, Kraków, Poland  
Email: marek.deja@uj.edu.pl
