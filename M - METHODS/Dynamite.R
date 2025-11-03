#######################################################################
# ADEMP — M: BAYESIAN DYNAMITE MODEL
# -----------------------------------------------------------------------
# Purpose:
#   Estimate lagged causal effects (X → Y and Y → X) using Bayesian
#   dynamic structural equation models implemented via the 'dynamite' package.
#
# Description:
#   - Fits autoregressive and cross-lagged Gaussian models for X and Y.
#   - Uses hierarchical Bayesian inference to estimate time-varying effects.
#   - Quantifies uncertainty via posterior distributions.
#   - Summarizes posterior means and 95% credible intervals.
#   - Evaluates bias, variance, and MSE against known ACEs from the
#     data-generating mechanisms (DGMs).
#######################################################################


# -----------------------------------------------------------------------
# 1. Helper function: Fit Bayesian dynamite model for one scenario
# -----------------------------------------------------------------------

library("dynamite")

fit_dynamite_scenario <- function(data, scenario, T_obs = 4, seed = 123) {

  set.seed(seed)

  # Define model using 'obs()' for each response variable
  # Each equation specifies lagged predictors and a Gaussian outcome
  model_formula <-
    obs(Y ~ lag(Y) + lag(X) + c, family = "gaussian") +
    obs(X ~ lag(X) + lag(Y) + c, family = "gaussian")

  # Fit Bayesian dynamic SEM
  dyn_fit <- dynamite(
    dformula = model_formula,
    data = data,
    time = "time",
    group = "id",
    chains = 2,
    iter = 1500,
    warmup = 750,
    control = list(
      max_treedepth = 12,
      adapt_delta = 0.95    # increase to reduce divergent transitions
    ),
    refresh = 0,
    verbose = FALSE
  )

  # Extract posterior summaries for regression coefficients
  post_sum <- as.data.frame(dyn_fit, types = "beta", summary = TRUE)

  # Return posterior means, SDs, and quantiles for cross-lagged paths
  tibble(
    scenario = scenario,
    beta_xy_hat = post_sum$mean[post_sum$parameter == "beta_Y_X_lag1"],
    beta_yx_hat = post_sum$mean[post_sum$parameter == "beta_X_Y_lag1"],
    beta_xy_sd  = post_sum$sd[post_sum$parameter == "beta_Y_X_lag1"],
    beta_yx_sd  = post_sum$sd[post_sum$parameter == "beta_X_Y_lag1"],
    beta_xy_q05 = post_sum$q5[post_sum$parameter == "beta_Y_X_lag1"],
    beta_xy_q95 = post_sum$q95[post_sum$parameter == "beta_Y_X_lag1"],
    beta_yx_q05 = post_sum$q5[post_sum$parameter == "beta_X_Y_lag1"],
    beta_yx_q95 = post_sum$q95[post_sum$parameter == "beta_X_Y_lag1"]
  ) %>%
    mutate(
      beta_xy_true = true_estimands(scenario)$ACE_XY,
      beta_yx_true = true_estimands(scenario)$ACE_YX
    )
}


# -----------------------------------------------------------------------
# 2. Run Bayesian Dynamite model across all scenarios (1–3)
# -----------------------------------------------------------------------
# Each scenario represents a distinct causal structure:
#   1 = No causal effects
#   2 = X → Y effect only
#   3 = Bidirectional effects (X ↔ Y)
# Output: Combined posterior summaries for all scenarios
# -----------------------------------------------------------------------

set.seed(20251025)
dyn_results <- map_dfr(1:3, function(s) {
  data_sim <- simulate_panel(scenario = s, N = 1000, T_obs = 4)
  fit_dynamite_scenario(data_sim, scenario = s, seed = 20251025 + s)
})


# -----------------------------------------------------------------------
# 3. Compute performance metrics from posterior means
# -----------------------------------------------------------------------
# Calculates:
#   - Bias (mean deviation from true ACE)
#   - Variance (empirical variance of posterior means)
#   - Mean Squared Error (bias² + variance)
# -----------------------------------------------------------------------

dyn_performance <- dyn_results %>%
  summarise(
    bias_XY = mean(beta_xy_hat - beta_xy_true),
    bias_YX = mean(beta_yx_hat - beta_yx_true),
    var_XY  = var(beta_xy_hat),
    var_YX  = var(beta_yx_hat),
    mse_XY  = mean((beta_xy_hat - beta_xy_true)^2),
    mse_YX  = mean((beta_yx_hat - beta_yx_true)^2)
  )


# -----------------------------------------------------------------------
# 4. Prepare formatted output tables
# -----------------------------------------------------------------------
# dyn_output : Posterior summaries by scenario and direction
# dyn_performance : Aggregate bias/variance/MSE metrics
# -----------------------------------------------------------------------

dyn_output <- dyn_results %>%
  mutate(
    scenario_label = recode(as.character(scenario),
      "1" = "S1: No causal effect",
      "2" = "S2: X → Y",
      "3" = "S3: X ↔ Y")
  ) %>%
  select(
    scenario_label,
    beta_xy_hat, beta_xy_sd, beta_xy_q05, beta_xy_q95, beta_xy_true,
    beta_yx_hat, beta_yx_sd, beta_yx_q05, beta_yx_q95, beta_yx_true
  )

print(dyn_output)
print(dyn_performance)


###############################################################################
# VISUALS FOR DYNAMITE MONTE CARLO RESULTS
# ---------------------------------------------------------------------------
# Inputs required:
#   - dyn_results: posterior results from fit_dynamite_scenario()
#   - metrics_flat: summary performance metrics (bias, coverage, etc.)
#   - meta_dir: output directory for figures
#
# Outputs:
#   - Fig1_EstimateDistributions.png
#   - Fig2_SECalibration.png
#   - Fig3_Coverage.png
#   - Fig4_ZipPlot.png
#   - Fig5_BiasMCSE.png
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Ensure parameter factor order for plotting
metrics_flat$parameter <- factor(metrics_flat$parameter,
                                 levels = c("β_Y←X", "β_X←Y"))

# ---------------------------------------------------------------------------
# 1) Posterior distributions of parameter estimates
# ---------------------------------------------------------------------------
p_dist <- dyn_results |>
  pivot_longer(c(beta_xy_hat, beta_yx_hat),
               names_to = "parameter", values_to = "estimate") |>
  mutate(parameter = recode(parameter,
                            beta_xy_hat = "β_Y←X",
                            beta_yx_hat = "β_X←Y")) |>
  ggplot(aes(x = estimate, fill = parameter)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~scenario, ncol = 1, scales = "free_y") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 13) +
  labs(title = "Posterior estimate distributions by scenario",
       x = "Estimate", y = "Density", fill = "Parameter")

ggsave(file.path(meta_dir, "Fig1_EstimateDistributions.png"),
       p_dist, width = 7, height = 8)

# ---------------------------------------------------------------------------
# 2) SE calibration: model-estimated SE vs empirical SE
# ---------------------------------------------------------------------------
p_sec <- metrics_flat |>
  ggplot(aes(x = factor(scenario), y = ratio,
             shape = parameter, color = parameter)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ratio - 1.96 * empSE_mcse / empSE,
                    ymax = ratio + 1.96 * empSE_mcse / empSE),
                width = 0.15) +
  theme_minimal(base_size = 13) +
  labs(title = "SE calibration (Model SE / Empirical SE)",
       x = "Scenario", y = "Ratio (ideal = 1)",
       shape = "Parameter", color = "Parameter")

ggsave(file.path(meta_dir, "Fig2_SECalibration.png"),
       p_sec, width = 6.5, height = 4.5)

# ---------------------------------------------------------------------------
# 3) Coverage: proportion of 95% credible intervals covering the truth
# ---------------------------------------------------------------------------
p_cov <- metrics_flat |>
  ggplot(aes(x = factor(scenario), y = cover,
             shape = parameter, color = parameter)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = pmax(0, cover - 1.96 * cover_mcse),
                    ymax = pmin(1, cover + 1.96 * cover_mcse)),
                width = 0.15) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  theme_minimal(base_size = 13) +
  labs(title = "95% coverage with Monte Carlo SE",
       x = "Scenario", y = "Coverage", shape = "Parameter",
       color = "Parameter")

ggsave(file.path(meta_dir, "Fig3_Coverage.png"),
       p_cov, width = 6.5, height = 4.5)

# ---------------------------------------------------------------------------
# 4) Zip plot: ranked 95% credible intervals vs true values
# ---------------------------------------------------------------------------
zip_df <- dyn_results |>
  select(scenario,
         xy_low = beta_xy_q05, xy_high = beta_xy_q95, xy_true = beta_xy_true,
         yx_low = beta_yx_q05, yx_high = beta_yx_q95, yx_true = beta_yx_true) |>
  pivot_longer(-scenario, names_to = c("param", ".value"),
               names_pattern = "(..)_?(low|high|true)?") |>
  group_by(scenario, param) |>
  arrange(scenario, param, low, .by_group = TRUE) |>
  mutate(rank = row_number())

p_zip <- ggplot(zip_df, aes(x = rank, ymin = low, ymax = high)) +
  geom_linerange(alpha = 0.6) +
  facet_grid(param ~ scenario, scales = "free_y",
             labeller = labeller(param = c(xy = "β_Y←X", yx = "β_X←Y"))) +
  geom_hline(data = distinct(zip_df, scenario, param, truth),
             aes(yintercept = truth), linetype = "dashed") +
  theme_minimal(base_size = 11) +
  labs(title = "Zip plot: ranked 95% credible intervals",
       x = "Rank (replication)", y = "95% CI")

ggsave(file.path(meta_dir, "Fig4_ZipPlot.png"),
       p_zip, width = 7, height = 6)

# ---------------------------------------------------------------------------
# 5) Optional: Bias estimates with Monte Carlo SE
# ---------------------------------------------------------------------------
p_bias <- metrics_flat |>
  ggplot(aes(x = factor(scenario), y = bias,
             shape = parameter, color = parameter)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = bias - 1.96 * bias_mcse,
                    ymax = bias + 1.96 * bias_mcse),
                width = 0.15) +
  theme_minimal(base_size = 13) +
  labs(title = "Bias with Monte Carlo SE",
       x = "Scenario", y = "Bias", shape = "Parameter",
       color = "Parameter")

ggsave(file.path(meta_dir, "Fig5_BiasMCSE.png"),
       p_bias, width = 6.5, height = 4.5)

