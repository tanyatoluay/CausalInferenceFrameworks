###############################################################################
# MONTE CARLO SIMULATION FOR DYNAMITE METHOD
# ADEMP: DGM → ESTIMAND → METHOD → PERFORMANCE
###############################################################################

library(tidyverse)
library(data.table)
library(purrr)
library(dynamite)

# ---------------------------------------------------------------------------
# 1. True causal effects (estimands)
# ---------------------------------------------------------------------------
true_estimands <- function(scenario) {
  tibble(
    scenario = scenario,
    ACE_XY = ifelse(scenario >= 2, 0.8, 0),
    ACE_YX = ifelse(scenario == 3, 0.5, 0)
  )
}

# ---------------------------------------------------------------------------
# 2. One Dynamite run
# ---------------------------------------------------------------------------
run_dyn_once <- function(scenario_id, N = 1000, T_obs = 4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data_sim <- simulate_panel(scenario = scenario_id, N = N, T_obs = T_obs)
  res <- fit_dynamite_scenario(data_sim, scenario = scenario_id, seed = seed)

  res %>%
    select(
      scenario,
      beta_xy_hat, beta_xy_q05, beta_xy_q95,
      beta_yx_hat, beta_yx_q05, beta_yx_q95,
      beta_xy_true, beta_yx_true
    )
}

# ---------------------------------------------------------------------------
# 3. Monte Carlo repetitions
# ---------------------------------------------------------------------------
set.seed(20251029)
nsim <- 10  # increase to 500+ for real runs

dyn_sims <- map_dfr(1:3, function(s) {
  reps <- replicate(
    nsim,
    run_dyn_once(s, seed = 20251029 + s * 1000 + sample.int(1e5, 1)),
    simplify = FALSE
  )
  bind_rows(reps)
})

# ---------------------------------------------------------------------------
# 4. Monte Carlo SE computation
# ---------------------------------------------------------------------------
mc_se <- function(estimates, theta, ci_low = NULL, ci_high = NULL) {
  nsim <- length(estimates)
  bias <- mean(estimates) - theta
  bias_mcse <- sqrt(var(estimates) / nsim)
  empSE <- sd(estimates)
  empSE_mcse <- empSE * sqrt(1 / (2 * (nsim - 1)))
  MSE <- mean((estimates - theta)^2)
  MSE_mcse <- sqrt(sum(((estimates - theta)^2 - MSE)^2) / (nsim * (nsim - 1)))
  if (!is.null(ci_low) && !is.null(ci_high)) {
    cover <- mean(ci_low <= theta & ci_high >= theta)
    cover_mcse <- sqrt(cover * (1 - cover) / nsim)
  } else cover <- cover_mcse <- NA
  data.frame(bias, bias_mcse, empSE, empSE_mcse, MSE, MSE_mcse, cover, cover_mcse)
}

# ---------------------------------------------------------------------------
# 5. Summarize performance per scenario
# ---------------------------------------------------------------------------
setDT(dyn_sims)

dyn_perf <- dyn_sims[, {
  theta_xy <- unique(beta_xy_true)
  theta_yx <- unique(beta_yx_true)
  xy <- mc_se(beta_xy_hat, theta_xy, beta_xy_q05, beta_xy_q95)
  yx <- mc_se(beta_yx_hat, theta_yx, beta_yx_q05, beta_yx_q95)
  cbind(method = c("XY", "YX"), rbind(xy, yx))
}, by = scenario]

# ---------------------------------------------------------------------------
# 6. Output
# ---------------------------------------------------------------------------
print(dyn_perf, digits = 4)
