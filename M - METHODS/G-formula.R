################################################################################
# ADEMP — M: g-FORMULA TESTING PIPELINE
#
# Purpose:
#   Run and evaluate the parametric g-formula under multiple simulated
#   causal scenarios and static treatment interventions.
#
#
# Dependencies:
#   gfoRmula, data.table, ggplot2
################################################################################

library(gfoRmula)
library(data.table)
library(tidyverse)

################################################################################
# 1. MASTER FUNCTION: RUN G-FORMULA UNDER DIFFERENT INTERVENTIONS
################################################################################
# Function: run_gformula_scenario()
#
# Arguments:
#   scenario_id   : integer, selects data-generating mechanism
#   N             : number of subjects
#   T_obs         : number of time points per subject
#   nsamples      : number of bootstrap samples for SE estimation
#   treat_levels  : numeric vector of static exposure levels to test
#   reverse       : logical, if TRUE flips causal direction (Y→X)
#
# Returns:
#   data.table with estimated g-formula results for all interventions,
#   labeled by scenario and causal direction.
#
# Workflow:
#   1. Simulate longitudinal data
#   2. Define causal direction (X→Y or Y→X)
#   3. Specify covariate and outcome models
#   4. Define static treatment interventions
#   5. Run parametric g-formula via gfoRmula::gformula()
#   6. Assemble and return results with readable labels
################################################################################

run_gformula_scenario <- function(scenario_id,
                                  N = 1000,
                                  T_obs = 4,
                                  nsamples = 200,
                                  treat_levels = c(0, 0.25, 0.5, 0.75, 1),
                                  reverse = FALSE) {
  # --- Step 1: Simulate longitudinal panel data
  dat <- simulate_panel(scenario = scenario_id, N = N, T_obs = T_obs)
  setDT(dat)
  K <- max(dat$time)

  # --- Step 2: Define causal direction
  # Forward (X→Y): treat X as exposure, Y as outcome
  # Reverse (Y→X): treat Y as exposure, X as outcome
  if (!reverse) {
    dat[time < K, Y := NA_real_]
    exposure <- "X"
    outcome_name <- "Y"
  } else {
    dat[time < K, X := NA_real_]
    exposure <- "Y"
    outcome_name <- "X"
  }

  # Model setup
  outcome_type <- "continuous_eof"
  covnames <- c("c", exposure)
  covtypes <- c("normal", "normal")
  id <- "id"
  time_name <- "time"
  time_points <- length(unique(dat[[time_name]]))

  # Define covariate histories
  histvars  <- list(c(exposure), c("c"))
  histories <- c(lagged, lagged)

  # --- Step 3: Covariate models
  covparams <- list(
    covmodels = c(
      c ~ lag1_c + time,
      as.formula(paste0(exposure, " ~ lag1_", exposure, " + lag1_c + c + time"))
    )
  )

  # --- Step 4: Static intervention definitions
  intvars_static <- rep(list(exposure), length(treat_levels))
  interventions_static <- lapply(treat_levels, function(val) {
    list(c(static, rep(val, time_points)))
  })
  int_descript_static <- paste("Static treat =", treat_levels)

  intvars       <- intvars_static
  interventions <- interventions_static
  int_descript  <- int_descript_static

  # --- Step 5: Outcome model
  ymodel <- as.formula(
    paste0(outcome_name, " ~ ", exposure, " + lag1_", exposure, " + c + lag1_c")
  )

  # --- Step 6: Run parametric g-formula
  message("Running scenario ", scenario_id,
          if (reverse) " (Y→X direction)" else " (X→Y direction)", " ...")

  res <- gformula(
    obs_data       = dat,
    id             = id,
    time_name      = time_name,
    time_points    = time_points,
    covnames       = covnames,
    covtypes       = covtypes,
    outcome_name   = outcome_name,
    outcome_type   = outcome_type,
    covparams      = covparams,
    ymodel         = ymodel,
    intvars        = intvars,
    interventions  = interventions,
    int_descript   = int_descript,
    histories      = histories,
    histvars       = histvars,
    nsimul         = 10000,
    nsamples       = nsamples,
    seed           = 1234
  )

  # --- Step 7: Format output
  out <- data.table(
    scenario   = scenario_id,
    direction  = ifelse(reverse, "Y→X", "X→Y"),
    res$result
  )

  # Align descriptive labels with numeric intervention codes
  int_labels_full <- c(NA, int_descript)
  out[, int_descript := int_labels_full[`Interv.` + 1]]
  out[`Interv.` == 0, int_descript := "Natural course"]

  return(out)
}


################################################################################
# 2. RUN ALL SCENARIOS AND BOTH DIRECTIONS
#
# Scenarios:
#   1 = No causal effect
#   2 = X→Y effect only
#   3 = Bidirectional effect (X↔Y)
#
# Combines results for both causal directions.
################################################################################

set.seed(20251031)

results_all <- rbindlist(c(
  lapply(1:3, run_gformula_scenario),                  # X→Y
  lapply(1:3, function(s) run_gformula_scenario(s, reverse = TRUE))  # Y→X
))


################################################################################
# 3. PERFORMANCE EVALUATION — ESTIMATE BIAS, VARIANCE, AND MSE
################################################################################

# Filter out natural course
gf_eval <- results_all[`Interv.` > 0, .(
  scenario,
  direction,
  int_descript,
  MD = `Mean difference`,
  MD_SE = `MD SE`,
  MD_low = `MD lower 95% CI`,
  MD_high = `MD upper 95% CI`
)]

# True ACEs by scenario and direction
gf_eval[, ACE_true := fifelse(direction == "X→Y" & scenario >= 2, 0.8,
                       fifelse(direction == "Y→X" & scenario == 3, 0.5, 0))]

# Derive performance metrics
gf_eval[, bias := MD - ACE_true]
gf_eval[, var := MD_SE^2]
gf_eval[, mse := bias^2 + var]

# Clean summary table
gf_summary <- gf_eval[, .(
  scenario,
  direction,
  int_descript,
  ACE_true,
  est_MD   = round(MD, 3),
  est_low  = round(MD_low, 3),
  est_high = round(MD_high, 3),
  bias     = round(bias, 3),
  var      = signif(var, 3),
  mse      = signif(mse, 3)
)]

print(gf_summary)


################################################################################
# 4. VISUALIZATION — STATIC TREATMENTS WITH 95% CIs
################################################################################

library(ggplot2)
library(data.table)

# Ordering for consistent display
results_all[, int_descript := factor(
  int_descript,
  levels = paste("Static treat =", c(0, 0.25, 0.5, 0.75, 1))
)]
results_all[, direction := factor(direction, levels = c("X→Y", "Y→X"))]

# Color palette (colorblind-friendly)
direction_colors <- c("X→Y" = "#E41A1C", "Y→X" = "#377EB8")

# Facet labels by scenario
scenario_labels <- c(
  `1` = "S1: No causal effect",
  `2` = "S2: X→Y causal effect",
  `3` = "S3: X↔Y causal effect"
)

# Plot estimated dose–response curves
p <- ggplot(results_all[`Interv.` > 0],
            aes(x = int_descript,
                y = `g-form mean`,
                color = direction,
                group = direction,
                shape = direction)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = `Mean lower 95% CI`,
                    ymax = `Mean upper 95% CI`),
                width = 0.1, alpha = 0.7) +
  facet_wrap(~ scenario, nrow = 1, labeller = labeller(scenario = scenario_labels)) +
  scale_color_manual(values = direction_colors) +
  scale_shape_manual(values = c(16, 17)) +
  labs(
    x = "Static Treatment Level",
    y = "Estimated g-formula Mean (95% CI)",
    color = "Causal Direction",
    shape = "Causal Direction",
    title = "Dose–Response Across Scenarios and Causal Directions"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p)
