# AD-DMT Health Economic Evaluation Framework

Interactive R Shiny platform for cost-effectiveness and budget impact analysis of lecanemab & donanemab in Ireland.

[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-Dashboard-orange.svg)](https://shiny.posit.co/)

## What this does

Decision-support tool for evaluating Alzheimer's disease-modifying therapies (DMTs) from the HSE payer perspective, built around HIQA guidelines. The platform has four modules:

1. **Cost-Effectiveness Analysis** -- 5-state Markov cohort model, PSA with up to 10k iterations, tornado diagrams, CEAC
2. **Budget Impact Analysis** -- population cascade from Irish epidemiology, multi-treatment market share over 5 years, cost decomposition by category
3. **Threshold Price Analysis** -- finds the maximum drug price at which each DMT becomes cost-effective
4. **Head-to-Head Comparison** -- lecanemab vs donanemab side by side

## Model

| | |
|---|---|
| Structure | 5-state Markov (MCI-AD, Mild, Moderate, Severe, Death) |
| Cycle length | Monthly, 30-year horizon |
| Perspective | HSE payer (societal optional) |
| Discount rate | 4% (HIQA) |
| WTP | EUR 45,000/QALY |
| Lecanemab | HR 0.69, 18 months (CLARITY-AD) |
| Donanemab | HR 0.73, 12 months (TRAILBLAZER-ALZ 2) |

## Results at list prices (Ireland)

| | Lecanemab | Donanemab |
|---|---|---|
| ICER | EUR 296,506/QALY | EUR 371,299/QALY |
| Inc. cost | EUR 47,110 | EUR 34,795 |
| Inc. QALYs | 0.159 | 0.094 |
| vs threshold | 6.6x | 8.3x |

Neither treatment is cost-effective at current list prices. Price reductions of ~38-72% would be needed depending on the country context.

## Running it

```r
# install dependencies
install.packages(c("shiny", "ggplot2", "scales", "DT", "shinythemes", "plotly"))

# run
shiny::runApp()
```

## Files

- `app.R` -- main Shiny application (CEA + BIA + comparison + price analysis)
- `R/` -- modular engine components (markov, PSA, BIA, tornado, parameters)
- `data/ireland_params.csv` -- all model inputs with sources
- `docs/` -- CHEERS checklist, technical validation notes

## Data sources

- Clinical: CLARITY-AD (NEJM 2023), TRAILBLAZER-ALZ 2 (JAMA 2023)
- Costs: HSE Fair Deal scheme, Connolly et al. 2014
- Utilities: Landeiro et al. 2020 (EQ-5D)
- Epidemiology: CSO Ireland, National Dementia Office

## Adding other countries

The framework takes country-specific inputs without code changes. Create a CSV following `data/ireland_params.csv` format with local costs, utilities, drug prices, and WTP threshold.

## Citation

Sen SE (2026). A Transferable Health Economic Evaluation Framework for Alzheimer's Disease-Modifying Therapies Across Small and Medium Healthcare Systems. *Submitted.*

## Related preprints

- Sen SE (2026). Health Economic Evaluation of Lecanemab for Early Alzheimer's Disease in Ireland.
- Sen SE (2026). Cost-Effectiveness of Donanemab for Early Alzheimer's Disease in Ireland.
- Sen SE (2026). Budget Impact Analysis of DMTs for Alzheimer's Disease in Ireland.

---

**Sevinc Elif Sen** -- Independent Health Economist | Decision-analytic modelling, HTA, R Shiny

MIT License
