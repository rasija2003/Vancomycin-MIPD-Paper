# Vancomycin-MIPD-Paper
Monte Carlo simulation study for Vancomycin                    precision dosing in sepsis — Pharm.D. research project
## Repository Files

| File | Description | Status |
|------|-------------|--------|
| 01_base_pk_model.R | Base one-compartment PK model | ✅ Complete |
| 02_covariate_model.R | Renal function covariate model | ✅ Complete |
| 03_monte_carlo_simulation.R | 3000 patient Monte Carlo simulation | ✅ Complete |
| 04_dose_optimisation.R | 23 regimen dose optimisation | ✅ Complete |
| pta_results.csv | PTA results by stratum | ✅ Complete |
| all_results.csv | Full simulation dataset | ✅ Complete |
| optimisation_r2.csv | Dose optimisation results | ✅ Complete |
| manuscript_draft.docx | Paper manuscript draft | 🔄 In progress |

## Key Findings So Far
- Standard vancomycin 1500mg q12h achieves PTA of only 
  1.9% (ARC), 23.2% (Normal), and 8.9% (AKI)
- No fixed dosing regimen achieved PTA ≥ 90% across 
  23 regimens and 7,000 virtual patients
- Findings support model-informed precision dosing (MIPD)
