###############################################################################
# ADEMP — M3: STRUCTURAL EQUATION MODEL (SEM)
# ---------------------------------------------------------------------------
# Purpose:
#   Estimate bidirectional causal effects between X and Y using SEM
#   applied to longitudinal panel data.
#
# Method:
#   Fit a cross-lagged path model with lagged predictors:
#       X_{t+1} ~ φ_x * X_t + β_yx * Y_t + γ_cx * c_t
#       Y_{t+1} ~ φ_y * Y_t + β_xy * X_t + γ_cy * c_t
#
#   SEM is fit via lavaan, jointly estimating all parameters.
#   Primary estimands: β_xy (X → Y) and β_yx (Y → X).
#
#   Performance metrics compare estimated βs to true ACEs from
#   the DGM (bias, variance, MSE).
###############################################################################

###############################################################################
# ADEMP — M3: STRUCTURAL EQUATION MODEL (SEM)
###############################################################################

library(lavaan)
library(dplyr)
library(purrr)
library(glue)

# ---------------------------------------------------------------------------
# Helper: True causal estimands by scenario
# ---------------------------------------------------------------------------
true_estimands <- function(scenario) {
  if (scenario == 1) {
    tibble(ACE_XY = 0, ACE_YX = 0)
  } else if (scenario == 2) {
    tibble(ACE_XY = 0.8, ACE_YX = 0)
  } else if (scenario == 3) {
    tibble(ACE_XY = 0.8, ACE_YX = 0.5)
  } else {
    stop("Scenario must be 1, 2, or 3")
  }
}

# ---------------------------------------------------------------------------
# 1. Fit SEM for a single scenario
# ---------------------------------------------------------------------------
fit_sem_scenario <- function(data, scenario, T_obs = 4) {

  # Convert to wide format for lavaan
  df_wide <- data %>%
    select(id, time, X, Y, c) %>%
    pivot_wider(names_from = time, values_from = c(X, Y, c), names_sep = "_")

  # Cross-lagged SEM specification
  sem_model <- "
    # Autoregressive paths
    X_2 ~ phi_x*X_1 + beta_yx*Y_1 + gamma_cx*c_2
    Y_2 ~ phi_y*Y_1 + beta_xy*X_1 + gamma_cy*c_2

    X_3 ~ phi_x*X_2 + beta_yx*Y_2 + gamma_cx*c_3
    Y_3 ~ phi_y*Y_2 + beta_xy*X_2 + gamma_cy*c_3

    X_4 ~ phi_x*X_3 + beta_yx*Y_3 + gamma_cx*c_4
    Y_4 ~ phi_y*Y_3 + beta_xy*X_3 + gamma_cy*c_4

    # Latent variances and covariances
    X_1 ~~ Y_1
    X_1 ~~ X_1
    Y_1 ~~ Y_1
    c_1 ~~ c_1
  "

  # Fit SEM
  fit <- sem(sem_model, data = df_wide, fixed.x = FALSE, missing = "ML")

  # Extract parameter estimates with CIs
  est <- parameterEstimates(fit, standardized = TRUE, ci = TRUE) %>%
    filter(label %in% c("beta_xy", "beta_yx", "phi_x", "phi_y")) %>%
    select(label, est, se, z, pvalue, ci.lower, ci.upper)

  # Append true causal effects
  truth <- true_estimands(scenario)

  tibble(
    scenario = scenario,
    beta_xy_hat  = est$est[est$label == "beta_xy"],
    beta_xy_se   = est$se[est$label == "beta_xy"],
    beta_xy_low  = est$ci.lower[est$label == "beta_xy"],
    beta_xy_high = est$ci.upper[est$label == "beta_xy"],
    beta_xy_true = truth$ACE_XY,

    beta_yx_hat  = est$est[est$label == "beta_yx"],
    beta_yx_se   = est$se[est$label == "beta_yx"],
    beta_yx_low  = est$ci.lower[est$label == "beta_yx"],
    beta_yx_high = est$ci.upper[est$label == "beta_yx"],
    beta_yx_true = truth$ACE_YX
  )
}

# ---------------------------------------------------------------------------
# 2. Fit SEM across all scenarios
# ---------------------------------------------------------------------------
set.seed(20251025)
sem_results <- map_dfr(1:3, function(s) {
  data_sim <- simulate_panel(scenario = s, N = 1000, T_obs = 4)
  fit_sem_scenario(data_sim, scenario = s)
})

# ---------------------------------------------------------------------------
# 3. Compute performance metrics
# ---------------------------------------------------------------------------
sem_performance <- sem_results %>%
  summarise(
    bias_XY = mean(beta_xy_hat - beta_xy_true),
    bias_YX = mean(beta_yx_hat - beta_yx_true),
    var_XY  = var(beta_xy_hat),
    var_YX  = var(beta_yx_hat),
    mse_XY  = mean((beta_xy_hat - beta_xy_true)^2),
    mse_YX  = mean((beta_yx_hat - beta_yx_true)^2)
  )

# ---------------------------------------------------------------------------
# 4. Output
# ---------------------------------------------------------------------------

# Add scenario labels and confidence intervals
sem_output <- sem_results %>%
  mutate(
    scenario_label = recode(as.character(scenario),
      "1" = "S1: No causal effect",
      "2" = "S2: X → Y",
      "3" = "S3: X ↔ Y")
  ) %>%
  select(scenario_label,
         beta_xy_hat, beta_xy_se, beta_xy_low, beta_xy_high, beta_xy_true,
         beta_yx_hat, beta_yx_se, beta_yx_low, beta_yx_high, beta_yx_true)

# Display results
print(sem_output)

# Display aggregate performance
sem_performance
