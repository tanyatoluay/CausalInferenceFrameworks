###############################################################################
# MONTE CARLO DRIVER — STRATOS SEM SIMULATION 
###############################################################################

library(tidyverse)
library(lavaan)

# --- Safe SEM wrapper --------------------------------------------------------
fit_sem_safe <- function(data, scenario, T_obs = 4) {
  out <- try(fit_sem_scenario(data, scenario = scenario, T_obs = T_obs), silent = TRUE)
  if (inherits(out, "try-error") || is.null(out)) {
    tibble(
      scenario = scenario,
      beta_xy_hat = NA_real_, beta_xy_se = NA_real_,
      beta_xy_low = NA_real_, beta_xy_high = NA_real_, beta_xy_true = true_estimands(scenario)$ACE_XY,
      beta_yx_hat = NA_real_, beta_yx_se = NA_real_,
      beta_yx_low = NA_real_, beta_yx_high = NA_real_, beta_yx_true = true_estimands(scenario)$ACE_YX,
      converged = FALSE
    )
  } else {
    out %>%
      mutate(converged = !any(is.na(c(beta_xy_hat, beta_yx_hat))))
  }
}

# --- One-scenario Monte Carlo ------------------------------------------------
mc_one_scenario <- function(
  scenario,
  n_sims = 500,
  N = 1000,
  T_obs = 4,
  phi_x = 0.5, phi_y = 0.5,
  beta_xy = NULL, beta_yx = NULL,
  gamma_cx = 0.3, gamma_cy = 0.4,
  sigma_x = 1, sigma_y = 1, sigma_c = 0.5,
  seed = 20251030
) {
  if (is.null(beta_xy)) beta_xy <- ifelse(scenario >= 2, 0.8, 0)
  if (is.null(beta_yx)) beta_yx <- ifelse(scenario == 3, 0.5, 0)

  results <- vector("list", n_sims)

  for (rep_id in seq_len(n_sims)) {
    set.seed(seed + 1000 * scenario + rep_id)
    dat <- simulate_panel(
      scenario = scenario, N = N, T_obs = T_obs,
      phi_x = phi_x, phi_y = phi_y,
      beta_xy = beta_xy, beta_yx = beta_yx,
      gamma_cx = gamma_cx, gamma_cy = gamma_cy,
      sigma_x = sigma_x, sigma_y = sigma_y, sigma_c = sigma_c
    )

    fit_res <- fit_sem_safe(dat, scenario, T_obs)
    results[[rep_id]] <- fit_res %>%
      mutate(rep = rep_id) %>%
      mutate(
        cover_xy = between(beta_xy_true, beta_xy_low, beta_xy_high),
        cover_yx = between(beta_yx_true, beta_yx_low, beta_yx_high),
        rej_xy = !(0 >= beta_xy_low & 0 <= beta_xy_high),
        rej_yx = !(0 >= beta_yx_low & 0 <= beta_yx_high)
      )
  }

  sims <- bind_rows(results)

  sims %>%
    mutate(scenario_label = recode(as.character(scenario),
      "1" = "S1: No causal effect",
      "2" = "S2: X → Y",
      "3" = "S3: X ↔ Y"
    ))
}

# --- Summarise performance ---------------------------------------------------
summarise_mc <- function(sims) {
  sims %>%
    group_by(scenario, scenario_label) %>%
    summarise(
      n_fit      = sum(converged, na.rm = TRUE),
      fail_rate  = mean(!converged, na.rm = TRUE),
      bias_XY    = mean(beta_xy_hat - beta_xy_true, na.rm = TRUE),
      sd_XY      = sd(beta_xy_hat, na.rm = TRUE),
      se_XY      = mean(beta_xy_se, na.rm = TRUE),
      rmse_XY    = sqrt(mean((beta_xy_hat - beta_xy_true)^2, na.rm = TRUE)),
      cover95_XY = mean(cover_xy, na.rm = TRUE),
      power_XY   = mean(rej_xy, na.rm = TRUE),
      bias_YX    = mean(beta_yx_hat - beta_yx_true, na.rm = TRUE),
      sd_YX      = sd(beta_yx_hat, na.rm = TRUE),
      se_YX      = mean(beta_yx_se, na.rm = TRUE),
      rmse_YX    = sqrt(mean((beta_yx_hat - beta_yx_true)^2, na.rm = TRUE)),
      cover95_YX = mean(cover_yx, na.rm = TRUE),
      power_YX   = mean(rej_yx, na.rm = TRUE),
      .groups = "drop"
    )
}

# --- Run all scenarios -------------------------------------------------------
run_mc_sem <- function(
  scenarios = 1:3,
  n_sims = 500,
  N = 1000,
  T_obs = 4,
  seed = 20251030
) {
  sims <- map_dfr(
    scenarios,
    ~ mc_one_scenario(scenario = .x, n_sims = n_sims, N = N, T_obs = T_obs, seed = seed)
  )
  list(sims = sims, summary = summarise_mc(sims))
}

# --- Example execution (RMarkdown chunk) ------------------------------------
set.seed(123)
mc_out <- run_mc_sem(n_sims = 10)

cat("\n### Monte Carlo Summary Table\n")
print(mc_out$summary, digits = 3)

cat("\n### First 10 Simulation Results\n")
print(head(mc_out$sims, 10), digits = 3)