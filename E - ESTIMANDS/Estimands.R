###############################################################################
# ADEMP — E: ESTIMANDS
# -----------------------------------------------------------------------------
# Goal:
#   Quantify the true (known) Average Causal Effects (ACE) embedded in the
#   Data-Generating Mechanism (DGM).
#
# Context:
#   From the DGM:
#     X_{t+1} = φ_X * X_t + β_YX * Y_t + γ_cX * c_t + ε_X
#     Y_{t+1} = φ_Y * Y_t + β_XY * X_t + γ_cY * c_t + ε_Y
#
#   => ACE_{X→Y}(t) = β_XY
#      ACE_{Y→X}(t) = β_YX
#
#   These coefficients are constant over time (time-invariant) but
#   confounder c_t varies over time, generating time-varying confounding.
#
# Purpose:
#   Provide reference values ("truth") to compare against estimator outputs
#   (G-Formula, DYNAMITE, SEM) in the Performance section.
###############################################################################

library(dplyr)
library(tibble)
library(purrr)
library(glue)

# ---------------------------------------------------------------------
# True structural coefficients per scenario
# ---------------------------------------------------------------------
true_estimands <- function(scenario, beta_xy = NULL, beta_yx = NULL) {
  bxy <- if (is.null(beta_xy)) ifelse(scenario >= 2, 0.8, 0) else beta_xy
  byx <- if (is.null(beta_yx)) ifelse(scenario == 3, 0.5, 0) else beta_yx
  list(ACE_XY = bxy, ACE_YX = byx)
}

# ---------------------------------------------------------------------
# Analytic ACE validation — exact, no simulation
# ---------------------------------------------------------------------
validate_estimands_analytic <- function(scenario, T_obs = 4, delta = 1) {
  par <- true_estimands(scenario)
  tibble(
    scenario = scenario,
    time_intervened = 0:(T_obs - 1),
    delta = delta,
    ACE_XY_true = par$ACE_XY,
    ACE_YX_true = par$ACE_YX,
    ACE_XY_estimated = par$ACE_XY * delta,
    ACE_YX_estimated = par$ACE_YX * delta,
    method = "analytic"
  )
}

# ---------------------------------------------------------------------
# Monte Carlo validation — checks analytic ACEs numerically
# ---------------------------------------------------------------------
validate_estimands_mc <- function(
  scenario = 3, N = 100000, T_obs = 4, t = 2, delta = 1,
  phi_x = 0.5, phi_y = 0.5, gamma_cx = 0.3, gamma_cy = 0.4,
  sigma_x = 1, sigma_y = 1, sigma_c = 0.5
) {
  par <- true_estimands(scenario)
  bxy <- par$ACE_XY; byx <- par$ACE_YX

  base <- simulate_panel(
    scenario = scenario, N = N, T_obs = T_obs,
    phi_x = phi_x, phi_y = phi_y,
    gamma_cx = gamma_cx, gamma_cy = gamma_cy,
    sigma_x = sigma_x, sigma_y = sigma_y, sigma_c = sigma_c
  )

  d <- base %>%
    filter(time %in% c(t, t + 1)) %>%
    select(id, time, X, Y, c) %>%
    tidyr::pivot_wider(names_from = time, values_from = c(X, Y, c), names_sep = "_")

  Xt <- glue("X_{t}")
  Yt <- glue("Y_{t}")
  Ct <- glue("c_{t}")  # use confounder at time t

  Y_next0 <- phi_y * d[[Yt]] + bxy * d[[Xt]] + gamma_cy * d[[Ct]]
  Y_next1 <- phi_y * d[[Yt]] + bxy * (d[[Xt]] + delta) + gamma_cy * d[[Ct]]

  X_next0 <- phi_x * d[[Xt]] + byx * d[[Yt]] + gamma_cx * d[[Ct]]
  X_next1 <- phi_x * d[[Xt]] + byx * (d[[Yt]] + delta) + gamma_cx * d[[Ct]]

  tibble(
    scenario = scenario,
    time_intervened = t,
    delta = delta,
    ACE_XY_true = bxy,
    ACE_YX_true = byx,
    ACE_XY_estimated = mean(Y_next1 - Y_next0),
    ACE_YX_estimated = mean(X_next1 - X_next0),
    method = "montecarlo"
  )
}

# ---------------------------------------------------------------------
# Generate analytic reference table
# ---------------------------------------------------------------------
results_estimands_full <- map_dfr(1:3, ~validate_estimands_analytic(.x, T_obs = 4, delta = 1))

# ---------------------------------------------------------------------
# Optional Monte Carlo sanity check
# ---------------------------------------------------------------------
mc_checks <- map_dfr(1:3, ~bind_rows(
  validate_estimands_mc(scenario = .x, t = 1, N = 20000),
  validate_estimands_mc(scenario = .x, t = 3, N = 20000)
))

# ---------------------------------------------------------------------
# Combine and print
# ---------------------------------------------------------------------
final_estimands <- results_estimands_full %>%
  mutate(
    scenario_label = recode(as.character(scenario),
      "1" = "S1: No causal effect",
      "2" = "S2: X → Y",
      "3" = "S3: X ↔ Y")
  ) %>%
  select(scenario_label, time_intervened, delta,
         ACE_XY_true, ACE_XY_estimated,
         ACE_YX_true, ACE_YX_estimated, method)

final_estimands
mc_checks