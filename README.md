# Transferable Health Economic Evaluation Framework for Alzheimer's Disease-Modifying Therapies

**Interactive R Shiny decision-support platform for cost-effectiveness and budget impact analysis of lecanemab and donanemab in Ireland**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-Dashboard-orange.svg)](https://shiny.posit.co/)

---

## Overview

This platform implements a **multi-component health economic evaluation framework** for assessing Alzheimer's disease-modifying therapies (DMTs) from the HSE payer perspective, following HIQA guidelines. It integrates:

- **Cost-Effectiveness Analysis (CEA)** — 5-state Markov cohort model with probabilistic sensitivity analysis (up to 10,000 Monte Carlo iterations)
- **Budget Impact Analysis (BIA)** — Population cascade, multi-treatment market share modelling, 5-year horizon with cost decomposition
- **Threshold Price Analysis** — Identifies maximum cost-effective drug prices at country-specific WTP thresholds
- **Head-to-Head Comparison** — Lecanemab vs. donanemab at list prices with visual benchmarking

The framework is designed for transferability to additional small and medium healthcare systems.

---

## Model Specifications

| Feature | Specification |
|---|---|
| Model type | Markov cohort state-transition |
| Health states | MCI-AD, Mild AD, Moderate AD, Severe AD, Death |
| Cycle length | Monthly (360 cycles, 30-year horizon) |
| Half-cycle correction | Applied |
| Perspective | HSE payer (societal perspective optional) |
| Discount rate | 4% per annum (HIQA reference case) |
| WTP threshold | EUR 45,000/QALY (HIQA) |
| PSA | Monte Carlo simulation, configurable iterations |
| BIA approach | ISPOR Principles of Good Practice |

## Clinical Data Sources

| Treatment | Trial | Hazard Ratio | Duration |
|---|---|---|---|
| Lecanemab | CLARITY-AD (van Dyck et al., NEJM 2023) | 0.69 | 18 months |
| Donanemab | TRAILBLAZER-ALZ 2 (Sims et al., JAMA 2023) | 0.73 | 12 months |

## Base-Case Results (Ireland, List Prices)

| Outcome | Lecanemab vs. SoC | Donanemab vs. SoC |
|---|---|---|
| ICER (EUR/QALY) | EUR 296,506 | EUR 371,299 |
| Incremental Cost | EUR 47,110 | EUR 34,795 |
| Incremental QALYs | 0.159 | 0.094 |
| Times WTP Threshold | 6.6x | 8.3x |
| Cost-effective at list price? | No | No |

---

## Quick Start

### Prerequisites

- R >= 4.1
- RStudio (recommended)

### Installation

```bash
git clone https://github.com/sevincelifalkan-coder/ad-dmt-framework.git
cd ad-dmt-framework
```

### Run Locally

```r
install.packages(c("shiny", "ggplot2", "scales", "DT", "shinythemes", "plotly"))
shiny::runApp()
```

---

## Repository Structure

```
ad-dmt-framework/
|-- app.R                  # Main Shiny application (integrated CEA + BIA)
|-- R/
|   |-- markov_engine.R    # 5-state Markov model engine
|   |-- psa_engine.R       # Probabilistic sensitivity analysis
|   |-- bia_engine.R       # Budget impact analysis with population cascade
|   |-- tornado.R          # One-way deterministic sensitivity analysis
|   |-- parameters.R       # Default parameters (Ireland)
|-- data/
|   |-- ireland_params.csv # Country-specific parameter set
|-- docs/
|   |-- CHEERS_checklist.md
|   |-- technical_report.md
|-- LICENSE
|-- README.md
```

---

## Methodology

### Cost-Effectiveness Analysis

The CEA module implements a 5-state Markov cohort model evaluating lecanemab and donanemab versus standard of care. Treatment effects are modelled as hazard ratio reductions on disease progression transition probabilities during the active treatment period.

**Cost layers:**
1. Drug acquisition
2. Diagnostic workup (amyloid PET, MRI, APOE4 genotyping, neuropsychological assessment)
3. Infusion/administration
4. Safety monitoring (serial MRI for ARIA)
5. ARIA management
6. Disease-state-dependent healthcare utilisation
7. Informal caregiver costs (societal perspective, optional)

### Budget Impact Analysis

The BIA module follows ISPOR guidelines with a population cascade approach:

Total Population > Aged 65+ > AD Prevalence > Early-Stage > Diagnosed > Amyloid Confirmed > Eligible for DMT > Treated

Three uptake scenarios (conservative, base case, optimistic) are modelled with dynamic market share allocation between lecanemab and donanemab over 5 years.

---

## Cost Data Sources (Ireland)

| Parameter | Value | Source |
|---|---|---|
| MCI-AD annual cost | EUR 4,200 | HSE, community care estimates |
| Mild AD annual cost | EUR 12,800 | Connolly et al. 2014, HSE |
| Moderate AD annual cost | EUR 28,500 | HSE home care packages |
| Severe AD annual cost | EUR 52,000 | HSE Fair Deal nursing home scheme |
| Lecanemab annual drug cost | EUR 24,766 | US WAC converted |
| Donanemab total treatment cost | EUR 24,500 | US WAC converted |
| Diagnostic workup | EUR 3,200 | PET + MRI + APOE4 + neuropsych |
| WTP threshold | EUR 45,000/QALY | HIQA guidelines |

---

## Extending to Other Countries

The framework is designed for parameterisation across healthcare systems. To add a new country:

1. Create a parameter CSV in data/ following ireland_params.csv format
2. Specify: transition probabilities, utilities, health state costs, drug prices, WTP threshold, discount rate, population data
3. The Markov engine and BIA cascade accept country-specific inputs without code modification

---

## Citation

If you use this framework in your research, please cite:

Sen, S.E. (2026). A Transferable Health Economic Evaluation Framework for Alzheimer's Disease-Modifying Therapies Across Small and Medium Healthcare Systems: A Multi-Country Decision-Analytic Modelling Study with Interactive R Shiny Platform.

R Shiny platform: https://github.com/sevincelifalkan-coder/ad-dmt-framework

---

## Related Work

- Sen, S.E. (2026). Health Economic Evaluation of Lecanemab for Early Alzheimer's Disease in Ireland. Preprint.
- Sen, S.E. (2026). Cost-Effectiveness of Donanemab for Early Alzheimer's Disease in Ireland. Preprint.
- Sen, S.E. (2026). Budget Impact Analysis of Disease-Modifying Therapies for Alzheimer's Disease in Ireland. Preprint.

---

## Author

**Sevinc Elif Sen, MSc**
Independent Health Economist & HEOR Researcher
Specialisation: Decision-analytic modelling | HTA | Cost-effectiveness analysis | R Shiny

---

## License

This project is licensed under the MIT License. See LICENSE for details.

(c) 2026 Sevinc Elif Sen. All rights reserved.
