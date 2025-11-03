###############################################################################
# ADEMP SIMULATION STUDY: DATA-GENERATING MECHANISM (DGM)
# -----------------------------------------------------------------------------
# Panel data for 3 continous variables X, Y and confounder C is created for
# - 3 scenarios:
#     1. No causal effect (confounding effect only)
#     2. X → Y causal effect
#     3. Bidirectional cross-lag effect X ↔ Y
#
# - Each subject has baseline + 4 follow-ups = 5 total time points
# - Confounder c_t is *time-varying* and autocorrelated (AR(1)), mimic of aging
# - Both exposure X_t and outcome Y_t depend on lagged confounder c_{t−1}
###############################################################################

simulate_panel <- function(
  scenario = 1,       # 1 = no effect, 2 = X→Y, 3 = X↔Y
  N = 1000,           # number of individuals
  T_obs = 4,          # number of follow-ups after baseline (→ 5 total time points)
  # Autoregressive parameters
  phi_x = 0.5,        # AR(1) parameter for X_t
  phi_y = 0.5,        # AR(1) parameter for Y_t
  rho_c = 0.7,        # AR(1) parameter for confounder c_t
  drift_c = 1.0,      # mean increase in confounder per follow-up (aging trend)
  # Causal path coefficients
  beta_xy = NULL,     # effect of X_{t−1} on Y_t
  beta_yx = NULL,     # effect of Y_{t−1} on X_t (only active in scenario 3)
  gamma_cx = 0.3,     # effect of c_{t−1} on X_t
  gamma_cy = 0.4,     # effect of c_{t−1} on Y_t
  # Noise terms
  sigma_x = 1,        # residual SD for X_t
  sigma_y = 1,        # residual SD for Y_t
  sigma_c = 0.5,      # residual SD for c_t
  # Baseline distribution for confounder
  mu_c0 = 30,         # mean of c_0
  sd_c0 = 5           # SD of c_0
) {

  # --- True causal parameters by scenario ------------------------------------
  if (is.null(beta_xy)) beta_xy <- ifelse(scenario >= 2, 0.8, 0)
  if (is.null(beta_yx)) beta_yx <- ifelse(scenario == 3, 0.5, 0)

  # --- Time points (baseline + follow-ups) -----------------------------------
  times <- 0:T_obs
  n_time <- length(times)

  # --- Initialize storage matrices -------------------------------------------
  C <- X <- Y <- matrix(NA_real_, nrow = N, ncol = n_time)

  # --- Baseline values -------------------------------------------------------
  C[, 1] <- rnorm(N, mu_c0, sd_c0)                     # baseline confounder
  X[, 1] <- gamma_cx * C[, 1] + rnorm(N, 0, sigma_x)   # baseline exposure
  Y[, 1] <- gamma_cy * C[, 1] + rnorm(N, 0, sigma_y)   # baseline outcome

  # --- Temporal evolution -----------------------------------------------------
for (t in 2:n_time) {
  # target means at t-1 and t
  mu_tm1 <- mu_c0 + drift_c * (t - 2)
  mu_t   <- mu_c0 + drift_c * (t - 1)

  # Aging-like process: everyone increases ~ +1 per wave with small noise
  C[, t] <- C[, t - 1] + drift_c + rnorm(N, 0, sigma_c)

  # exposure and outcome with lagged confounder
  X[, t] <- phi_x * X[, t - 1] + beta_yx * Y[, t - 1] +
            gamma_cx * C[, t - 1] + rnorm(N, 0, sigma_x)

  Y[, t] <- phi_y * Y[, t - 1] + beta_xy * X[, t - 1] +
            gamma_cy * C[, t - 1] + rnorm(N, 0, sigma_y)
}


  # --- Return in long format -------------------------------------------------
  data.frame(
    id   = rep(1:N, each = n_time),
    time = rep(times, times = N),
    c    = as.vector(t(C)),
    X    = as.vector(t(X)),
    Y    = as.vector(t(Y)),
    scenario = paste0("Scenario_", scenario)
  )
}