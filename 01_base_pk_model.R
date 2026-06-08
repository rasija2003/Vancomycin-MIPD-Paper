# ================================================
# PROJECT: Vancomycin MIPD Monte Carlo Study
# FILE:    01_base_pk_model.R
# AUTHOR:  Sivakumar Kalyana sundaram
# DATE:    June 2026
# PURPOSE: Base one-compartment vancomycin PK model
#          Simulates typical adult with normal renal
#          function — 1500mg IV q12h over 48 hours
#          Generates Figure 1 of manuscript
# REF:     Rybak et al. 2020 ASHP/IDSA/SIDP
#          Vancomycin Guidelines
# ================================================

library(mrgsolve)
library(ggplot2)

library(mrgsolve)
library(ggplot2)

# Step 1 — Define the PK model
van_model <- '
$PARAM
  CL = 4.5,
  V  = 50.0

$CMT CENT

$ODE
  dxdt_CENT = -(CL/V) * CENT;

$TABLE
  double CP = CENT/V;

$CAPTURE CP
'

# Step 2 — Build the model
mod <- mread("vancomycin", tempdir(), van_model)

# Step 3 — Design the dosing regimen
# 1500mg IV every 12 hours, 1-hour infusion, 4 doses total
dose <- ev(
  amt  = 1500,
  rate = 1500,
  ii   = 12,
  addl = 3
)

# Step 4 — Simulate one typical patient
sim_result <- mod %>%
  ev(dose) %>%
  mrgsim(end = 48, delta = 0.1) %>%
  as.data.frame()

# Step 5 — Plot
ggplot(sim_result, aes(x = time, y = CP)) +
  geom_line(color = "#1B3A6B", linewidth = 1.2) +
  geom_hline(yintercept = 10, linetype = "dashed",
             color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 20, linetype = "dashed",
             color = "darkred", linewidth = 0.8) +
  annotate("text", x = 45, y = 10.8,
           label = "Trough target: 10 mg/L",
           color = "red", size = 3.5) +
  annotate("text", x = 45, y = 20.8,
           label = "Toxicity threshold: 20 mg/L",
           color = "darkred", size = 3.5) +
  labs(
    title    = "Vancomycin Plasma Concentration — Typical Adult Patient",
    subtitle = "1500mg IV q12h | CL = 4.5 L/hr | Vd = 50L | Normal Renal Function",
    x        = "Time (hours)",
    y        = "Vancomycin Concentration (mg/L)"
  ) +
  theme_bw(base_size = 13) +
  scale_x_continuous(breaks = seq(0, 48, by = 6)) +
  scale_y_continuous(limits = c(0, 45))
