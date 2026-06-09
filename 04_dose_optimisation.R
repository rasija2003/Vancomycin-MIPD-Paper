# ================================================
# PROJECT: Vancomycin MIPD Monte Carlo Study
# FILE:    04_dose_optimisation.R
# AUTHOR:  Raksana Sivakumar
# DATE:    June 2026
# PURPOSE: Systematic dose optimisation
#          Find regimen achieving PTA >= 90%
#          per renal function stratum
#          Generates Table 3 and Figure 4
# REF:     Rybak 2020 ASHP/IDSA/SIDP
# ================================================

library(mrgsolve)
library(ggplot2)
library(dplyr)

set.seed(42)

# ── STEP 1: Load saved data if available ─────────────────────
# If pta_results not in memory, reload from CSV

if (!exists("all_results")) {
  cat("Loading saved simulation data...\n")
  all_results  <- read.csv("data/all_results.csv")
  pta_results  <- read.csv("data/pta_results.csv")
  cat("Data loaded successfully\n")
}

# ── STEP 2: Define PK model ───────────────────────────────────

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
  double CP = CENT/V;

$CAPTURE CP CL V CrCl
'

mod_opt <- mread("van_opt", tempdir(), van_mc_model)

# ── STEP 3: Define variability parameters ─────────────────────

n_patients <- 500  # 500 per regimen is sufficient for optimisation
CV_CL      <- 0.28
CV_V       <- 0.23
omega_CL   <- sqrt(log(1 + CV_CL^2))
omega_V    <- sqrt(log(1 + CV_V^2))

# ── STEP 4: Core simulation function ─────────────────────────

simulate_regimen <- function(n, crcl_mean, crcl_sd,
                             dose_amt, dose_interval,
                             stratum_label, regimen_label) {
  
  # Generate patient parameters
  CrCl_values <- rnorm(n, mean = crcl_mean, sd = crcl_sd)
  CrCl_values <- pmax(CrCl_values, 5)
  CrCl_values <- pmin(CrCl_values, 250)
  ETA1_values <- rnorm(n, mean = 0, sd = omega_CL)
  ETA2_values <- rnorm(n, mean = 0, sd = omega_V)
  
  # Calculate number of doses for 48 hours
  n_additional <- floor(48 / dose_interval) - 1
  
  # Dosing event
  dose <- ev(
    amt  = dose_amt,
    rate = dose_amt,  # 1-hour infusion
    ii   = dose_interval,
    addl = n_additional
  )
  
  # Simulation end time
  end_time <- 48
  
  # Steady state window — last dosing interval
  ss_start <- end_time - dose_interval
  
  results <- data.frame()
  
  for (i in 1:n) {
    
    sim <- mod_opt %>%
      param(
        CrCl = CrCl_values[i],
        ETA1 = ETA1_values[i],
        ETA2 = ETA2_values[i]
      ) %>%
      ev(dose) %>%
      mrgsim(end = end_time, delta = 0.25) %>%
      as.data.frame()
    
    # AUC at steady state
    ss_data <- sim %>%
      filter(time >= ss_start & time <= end_time)
    
    AUC_ss <- sum(diff(ss_data$time) *
                    (head(ss_data$CP, -1) +
                       tail(ss_data$CP, -1)) / 2)
    
    # Trough at steady state
    trough_ss <- sim %>%
      filter(time >= (ss_start - 0.1) &
               time <= (ss_start + 0.1)) %>%
      summarise(trough = min(CP)) %>%
      pull(trough)
    
    results <- bind_rows(results, data.frame(
      Stratum  = stratum_label,
      Regimen  = regimen_label,
      Dose     = dose_amt,
      Interval = dose_interval,
      CrCl     = round(CrCl_values[i], 1),
      AUC_24   = round(AUC_ss, 1),
      Trough   = round(trough_ss, 1)
    ))
  }
  return(results)
}

# ── STEP 5: Define all regimens to test ───────────────────────

cat("Starting dose optimisation...\n")
cat("This will take 5-8 minutes — please wait\n\n")

# ── ARC REGIMENS ──────────────────────────────────────────────
cat("Testing ARC regimens...\n")

arc_1 <- simulate_regimen(n_patients, 150, 20,
                          1500, 12, "ARC", "1500mg q12h")
arc_2 <- simulate_regimen(n_patients, 150, 20,
                          2000, 12, "ARC", "2000mg q12h")
arc_3 <- simulate_regimen(n_patients, 150, 20,
                          2500, 12, "ARC", "2500mg q12h")
arc_4 <- simulate_regimen(n_patients, 150, 20,
                          1500, 8,  "ARC", "1500mg q8h")
arc_5 <- simulate_regimen(n_patients, 150, 20,
                          2000, 8,  "ARC", "2000mg q8h")

# ── NORMAL REGIMENS ───────────────────────────────────────────
cat("Testing Normal regimens...\n")

nor_1 <- simulate_regimen(n_patients, 90, 15,
                          1500, 12, "Normal", "1500mg q12h")
nor_2 <- simulate_regimen(n_patients, 90, 15,
                          2000, 12, "Normal", "2000mg q12h")
nor_3 <- simulate_regimen(n_patients, 90, 15,
                          2500, 12, "Normal", "2500mg q12h")
nor_4 <- simulate_regimen(n_patients, 90, 15,
                          1500, 8,  "Normal", "1500mg q8h")

# ── AKI REGIMENS ──────────────────────────────────────────────
cat("Testing AKI regimens...\n")

aki_1 <- simulate_regimen(n_patients, 20, 5,
                          1500, 12, "AKI", "1500mg q12h")
aki_2 <- simulate_regimen(n_patients, 20, 5,
                          750,  12, "AKI", "750mg q12h")
aki_3 <- simulate_regimen(n_patients, 20, 5,
                          500,  12, "AKI", "500mg q12h")
aki_4 <- simulate_regimen(n_patients, 20, 5,
                          750,  24, "AKI", "750mg q24h")
aki_5 <- simulate_regimen(n_patients, 20, 5,
                          500,  24, "AKI", "500mg q24h")

# ── STEP 6: Combine all results ───────────────────────────────

all_regimens <- bind_rows(
  arc_1, arc_2, arc_3, arc_4, arc_5,
  nor_1, nor_2, nor_3, nor_4,
  aki_1, aki_2, aki_3, aki_4, aki_5
)

# ── STEP 7: Calculate PTA for each regimen ────────────────────

optimisation_results <- all_regimens %>%
  group_by(Stratum, Regimen, Dose, Interval) %>%
  summarise(
    N            = n(),
    Mean_AUC     = round(mean(AUC_24), 1),
    SD_AUC       = round(sd(AUC_24), 1),
    Mean_Trough  = round(mean(Trough), 1),
    PTA_target   = round(mean(AUC_24 >= 400 &
                                AUC_24 <= 600) * 100, 1),
    PTA_above600 = round(mean(AUC_24 > 600) * 100, 1),
    PTA_below400 = round(mean(AUC_24 < 400) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(Stratum, Interval, Dose)

# ── STEP 8: Print results ─────────────────────────────────────

cat("\n=============================================\n")
cat("DOSE OPTIMISATION RESULTS\n")
cat("Target: PTA >= 90% within AUC 400-600\n")
cat("=============================================\n\n")

cat("--- ARC PATIENTS (CrCl ~150 mL/min) ---\n")
print(optimisation_results %>%
        filter(Stratum == "ARC") %>%
        select(Regimen, Mean_AUC, PTA_target,
               PTA_above600, PTA_below400))

cat("\n--- NORMAL PATIENTS (CrCl ~90 mL/min) ---\n")
print(optimisation_results %>%
        filter(Stratum == "Normal") %>%
        select(Regimen, Mean_AUC, PTA_target,
               PTA_above600, PTA_below400))

cat("\n--- AKI PATIENTS (CrCl ~20 mL/min) ---\n")
print(optimisation_results %>%
        filter(Stratum == "AKI") %>%
        select(Regimen, Mean_AUC, PTA_target,
               PTA_above600, PTA_below400))

# ── STEP 9: Identify optimal regimens ────────────────────────

cat("\n=============================================\n")
cat("OPTIMAL REGIMENS (PTA >= 90%)\n")
cat("=============================================\n")

optimal <- optimisation_results %>%
  filter(PTA_target >= 90) %>%
  select(Stratum, Regimen, Mean_AUC,
         PTA_target, PTA_above600, PTA_below400)

print(optimal)

# ── STEP 10: Figure 4 — PTA comparison bar chart ─────────────

# Set factor order for display
optimisation_results$Stratum <- factor(
  optimisation_results$Stratum,
  levels = c("ARC", "Normal", "AKI")
)

ggplot(optimisation_results,
       aes(x    = reorder(Regimen, PTA_target),
           y    = PTA_target,
           fill = Stratum)) +
  geom_bar(stat = "identity", width = 0.7) +
  
  # 90% threshold line
  geom_hline(yintercept = 90,
             linetype   = "dashed",
             color      = "red",
             linewidth  = 1) +
  
  annotate("text", x = 1, y = 92,
           label = "PTA 90% threshold",
           color = "red", size = 3.5,
           hjust = 0) +
  
  facet_wrap(~Stratum, scales = "free_x") +
  
  scale_fill_manual(values = c(
    "ARC"    = "#2196F3",
    "Normal" = "#1B3A6B",
    "AKI"    = "#C62828"
  )) +
  
  labs(
    title    = "Vancomycin Dose Optimisation — PTA by Regimen",
    subtitle = "Target: PTA ≥ 90% within AUC 400–600 mg·h/L",
    x        = "Dosing Regimen",
    y        = "Probability of Target Attainment (%)",
    fill     = "Renal Stratum"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 35,
                                    hjust = 1,
                                    size  = 9),
    strip.text       = element_text(face = "bold"),
    legend.position  = "none",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40")
  ) +
  
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 10))

# ── STEP 11: Save optimisation results ───────────────────────

write.csv(optimisation_results,
          "data/optimisation_results.csv",
          row.names = FALSE)

cat("\nOptimisation results saved to data/\n")