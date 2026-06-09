# ================================================
# PROJECT: Vancomycin MIPD Monte Carlo Study
# FILE:    02_covariate_model.R
# AUTHOR:  Sivakumar [your surname]
# DATE:    June 2026
# PURPOSE: Add renal function covariate to PK model
#          Simulate ARC vs Normal vs AKI patients
#          Generates Figure 2 of manuscript
# REF:     Matzke 1984, Rybak 2020 ASHP/IDSA/SIDP
# ================================================

library(mrgsolve)
library(ggplot2)
library(dplyr)

# ── STEP 1: Define model with CrCl covariate ──────────────────

van_covariate_model <- '
$PARAM
  TVCL = 4.5,
  V    = 50.0,
  CrCl = 90

$CMT CENT

$MAIN
  double CL = TVCL * pow(CrCl/90.0, 0.8);

$ODE
  dxdt_CENT = -(CL/V) * CENT;

$TABLE
  double CP = CENT/V;

$CAPTURE CP CL CrCl
'

# ── STEP 2: Build the model ───────────────────────────────────

mod2 <- mread("van_covariate", tempdir(), van_covariate_model)

# ── STEP 3: Define dosing regimen ─────────────────────────────
# Same dose for all three patients — 1500mg IV q12h x 4 doses

dose <- ev(
  amt  = 1500,
  rate = 1500,
  ii   = 12,
  addl = 3
)

# ── STEP 4: Simulate three patient types ──────────────────────

# Patient 1 — ARC (Augmented Renal Clearance)
# Young septic patient, hyperdynamic circulation
sim_ARC <- mod2 %>%
  param(CrCl = 150) %>%
  ev(dose) %>%
  mrgsim(end = 48, delta = 0.1) %>%
  as.data.frame() %>%
  mutate(Patient = "ARC (CrCl = 150 mL/min)")

# Patient 2 — Normal renal function
sim_Normal <- mod2 %>%
  param(CrCl = 90) %>%
  ev(dose) %>%
  mrgsim(end = 48, delta = 0.1) %>%
  as.data.frame() %>%
  mutate(Patient = "Normal (CrCl = 90 mL/min)")

# Patient 3 — AKI (Acute Kidney Injury)
# Critically ill, organ dysfunction
sim_AKI <- mod2 %>%
  param(CrCl = 20) %>%
  ev(dose) %>%
  mrgsim(end = 48, delta = 0.1) %>%
  as.data.frame() %>%
  mutate(Patient = "AKI (CrCl = 20 mL/min)")

# ── STEP 5: Combine all three simulations ─────────────────────

sim_combined <- bind_rows(sim_ARC, sim_Normal, sim_AKI)

# Set the order for the legend
sim_combined$Patient <- factor(sim_combined$Patient,
                               levels = c(
                                 "ARC (CrCl = 150 mL/min)",
                                 "Normal (CrCl = 90 mL/min)",
                                 "AKI (CrCl = 20 mL/min)"
                               )
)

# ── STEP 6: Plot Figure 2 ──────────────────────────────────────

ggplot(sim_combined, aes(x = time, y = CP,
                         color = Patient,
                         linetype = Patient)) +
  geom_line(linewidth = 1.2) +
  
  # Therapeutic window
  geom_hline(yintercept = 10, linetype = "dashed",
             color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 20, linetype = "dashed",
             color = "darkred", linewidth = 0.8) +
  
  # Shaded therapeutic window
  annotate("rect",
           xmin = 0, xmax = 48,
           ymin = 10, ymax = 20,
           alpha = 0.08, fill = "green") +
  
  # Labels for reference lines
  annotate("text", x = 46, y = 9,
           label = "Min target: 10 mg/L",
           color = "red", size = 3.2, hjust = 1) +
  annotate("text", x = 46, y = 21,
           label = "Toxicity: 20 mg/L",
           color = "darkred", size = 3.2, hjust = 1) +
  
  # Colours for three patient types
  scale_color_manual(values = c(
    "ARC (CrCl = 150 mL/min)"    = "#2196F3",  # Blue
    "Normal (CrCl = 90 mL/min)"  = "#1B3A6B",  # Navy
    "AKI (CrCl = 20 mL/min)"     = "#C62828"   # Red
  )) +
  
  scale_linetype_manual(values = c(
    "ARC (CrCl = 150 mL/min)"    = "solid",
    "Normal (CrCl = 90 mL/min)"  = "solid",
    "AKI (CrCl = 20 mL/min)"     = "solid"
  )) +
  
  labs(
    title    = "Vancomycin Concentration — Impact of Renal Function in Sepsis",
    subtitle = "1500mg IV q12h | Same dose, three patient types",
    x        = "Time (hours)",
    y        = "Vancomycin Concentration (mg/L)",
    color    = "Patient Type",
    linetype = "Patient Type"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold"),
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(color = "gray40")
  ) +
  
  scale_x_continuous(breaks = seq(0, 48, by = 6)) +
  scale_y_continuous(limits = c(0, 80))

# ── STEP 7: Print clearance values for each patient ───────────

cat("=== Calculated Clearance Values ===\n")
cat("ARC    (CrCl=150): CL =",
    round(4.5 * (150/90)^0.8, 2), "L/hr\n")
cat("Normal (CrCl=90):  CL =",
    round(4.5 * (90/90)^0.8,  2), "L/hr\n")
cat("AKI    (CrCl=20):  CL =",
    round(4.5 * (20/90)^0.8,  2), "L/hr\n")