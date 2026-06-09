# ================================================
# PROJECT: Vancomycin MIPD Monte Carlo Study
# FILE:    03_monte_carlo_simulation.R
# AUTHOR : Raksana Sivakumar
# DATE:    June 2026
# PURPOSE: Monte Carlo simulation — 1000 virtual patients
#          Calculate AUC and Probability of Target
#          Attainment (PTA) by renal function stratum
# REF:     Rybak 2020, Matzke 1984
# ================================================

library(mrgsolve)
library(ggplot2)
library(dplyr)

set.seed(42)  # Makes results reproducible — same every time you run

# ── STEP 1: Define number of patients ────────────────────────

n_patients <- 1000

# ── STEP 2: Define population parameters ─────────────────────
# Inter-individual variability (IIV) from published literature
# Expressed as coefficient of variation (CV%)
# CL CV = 28%, V CV = 23% (Rybak 2020 / Winter 2010)

TVCL <- 4.5   # Typical clearance L/hr at CrCl=90
TVV  <- 50.0  # Typical volume L
CV_CL <- 0.28  # 28% coefficient of variation for CL
CV_V  <- 0.23  # 23% coefficient of variation for V

# Convert CV to log-normal variance (omega squared)
omega_CL <- sqrt(log(1 + CV_CL^2))
omega_V  <- sqrt(log(1 + CV_V^2))

# ── STEP 3: Define the PK model with variability ─────────────

van_mc_model <- '
$PARAM
  TVCL = 4.5,
  TVV  = 50.0,
  CrCl = 90,
  ETA1 = 0,
  ETA2 = 0

$CMT CENT

$MAIN
  double CL = TVCL * pow(CrCl/90.0, 0.8) * exp(ETA1);
  double V  = TVV * exp(ETA2);

$ODE
  dxdt_CENT = -(CL/V) * CENT;

$TABLE
  double CP  = CENT/V;

$CAPTURE CP CL V CrCl
'

mod_mc <- mread("van_mc", tempdir(), van_mc_model)

# ── STEP 4: Define dosing regimen ─────────────────────────────
# 1500mg IV q12h x 4 doses (48 hours)

dose <- ev(
  amt  = 1500,
  rate = 1500,
  ii   = 12,
  addl = 3
)

# ── STEP 5: Monte Carlo simulation function ───────────────────
# Simulates n patients for a given CrCl range

simulate_population <- function(n, crcl_mean, crcl_sd,
                                stratum_label) {
  
  cat("Simulating", n, "patients —", stratum_label, "...\n")
  
  # Generate individual patient parameters
  # CrCl varies within each stratum (realistic variability)
  CrCl_values <- rnorm(n,
                       mean = crcl_mean,
                       sd   = crcl_sd)
  
  # Keep CrCl within realistic bounds
  CrCl_values <- pmax(CrCl_values, 5)
  CrCl_values <- pmin(CrCl_values, 250)
  
  # Generate individual variability for CL and V
  ETA1_values <- rnorm(n, mean = 0, sd = omega_CL)
  ETA2_values <- rnorm(n, mean = 0, sd = omega_V)
  
  # Store results
  results <- data.frame()
  
  for (i in 1:n) {
    
    # Simulate one patient
    sim <- mod_mc %>%
      param(
        CrCl = CrCl_values[i],
        ETA1 = ETA1_values[i],
        ETA2 = ETA2_values[i]
      ) %>%
      ev(dose) %>%
      mrgsim(end = 48, delta = 0.25) %>%
      as.data.frame()
    
    # Calculate AUC at steady state (last dosing interval)
    # Using trapezoidal rule on hours 36-48 (dose 4 interval)
    ss_data <- sim %>%
      filter(time >= 36 & time <= 48)
    
    # Trapezoidal AUC calculation
    AUC_ss <- sum(diff(ss_data$time) *
                    (head(ss_data$CP, -1) +
                       tail(ss_data$CP, -1)) / 2)
    
    # Trough at steady state (just before dose 4)
    trough_ss <- sim %>%
      filter(time >= 35.9 & time <= 36.1) %>%
      summarise(trough = min(CP)) %>%
      pull(trough)
    
    # Peak at steady state (dose 4)
    peak_ss <- sim %>%
      filter(time >= 36 & time <= 37) %>%
      summarise(peak = max(CP)) %>%
      pull(peak)
    
    # Get individual CL and V
    ind_CL <- unique(sim$CL)[1]
    ind_V  <- unique(sim$V)[1]
    
    results <- bind_rows(results, data.frame(
      Patient  = i,
      Stratum  = stratum_label,
      CrCl     = round(CrCl_values[i], 1),
      CL       = round(ind_CL, 2),
      V        = round(ind_V, 1),
      AUC_24   = round(AUC_ss, 1),
      Trough   = round(trough_ss, 1),
      Peak     = round(peak_ss, 1)
    ))
  }
  
  return(results)
}

# ── STEP 6: Run simulation for all three strata ───────────────

# ARC: CrCl mean=150, SD=20
results_ARC <- simulate_population(
  n            = n_patients,
  crcl_mean    = 150,
  crcl_sd      = 20,
  stratum_label = "ARC (CrCl ~150)"
)

# Normal: CrCl mean=90, SD=15
results_Normal <- simulate_population(
  n            = n_patients,
  crcl_mean    = 90,
  crcl_sd      = 15,
  stratum_label = "Normal (CrCl ~90)"
)

# AKI: CrCl mean=20, SD=5
results_AKI <- simulate_population(
  n            = n_patients,
  crcl_mean    = 20,
  crcl_sd      = 5,
  stratum_label = "AKI (CrCl ~20)"
)

# Combine all results
all_results <- bind_rows(results_ARC,
                         results_Normal,
                         results_AKI)

# ── STEP 7: Calculate PTA ─────────────────────────────────────
# Target: AUC 400-600 mg.h/L (ASHP/IDSA 2020)

pta_results <- all_results %>%
  group_by(Stratum) %>%
  summarise(
    N               = n(),
    Mean_CrCl       = round(mean(CrCl), 1),
    Mean_CL         = round(mean(CL), 2),
    Mean_AUC        = round(mean(AUC_24), 1),
    SD_AUC          = round(sd(AUC_24), 1),
    Mean_Trough     = round(mean(Trough), 1),
    Mean_Peak       = round(mean(Peak), 1),
    PTA_target      = round(mean(AUC_24 >= 400 &
                                   AUC_24 <= 600) * 100, 1),
    PTA_above_600   = round(mean(AUC_24 > 600) * 100, 1),
    PTA_below_400   = round(mean(AUC_24 < 400) * 100, 1)
  )

cat("\n========================================\n")
cat("PROBABILITY OF TARGET ATTAINMENT RESULTS\n")
cat("Dose: 1500mg IV q12h\n")
cat("Target: AUC 400-600 mg.h/L (ASHP/IDSA 2020)\n")
cat("========================================\n")
print(pta_results)

# ── STEP 8: Plot Figure 3 — AUC distributions ─────────────────

all_results$Stratum <- factor(all_results$Stratum,
                              levels = c("ARC (CrCl ~150)",
                                         "Normal (CrCl ~90)",
                                         "AKI (CrCl ~20)"))

ggplot(all_results, aes(x = AUC_24, fill = Stratum)) +
  geom_histogram(binwidth = 30, alpha = 0.7,
                 color = "white", position = "identity") +
  facet_wrap(~Stratum, ncol = 1) +
  
  # Target window
  geom_vline(xintercept = 400, linetype = "dashed",
             color = "red", linewidth = 0.9) +
  geom_vline(xintercept = 600, linetype = "dashed",
             color = "darkred", linewidth = 0.9) +
  
  # Shaded target zone
  annotate("rect",
           xmin = 400, xmax = 600,
           ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "green") +
  
  scale_fill_manual(values = c(
    "ARC (CrCl ~150)"    = "#2196F3",
    "Normal (CrCl ~90)"  = "#1B3A6B",
    "AKI (CrCl ~20)"     = "#C62828"
  )) +
  
  labs(
    title    = "Vancomycin AUC Distribution — Monte Carlo Simulation",
    subtitle = "1000 virtual patients per stratum | 1500mg IV q12h | Target: AUC 400-600 mg·h/L",
    x        = "AUC₀₋₂₄ (mg·h/L)",
    y        = "Number of Patients",
    fill     = "Renal Stratum"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 11),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40")
  ) +
  
  scale_x_continuous(
    limits = c(0, 2000),
    breaks = seq(0, 2000, by = 200)
  )