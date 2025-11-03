###############################################################################
# ADEMP — D: DATASET GENERATION PIPELINE
# -----------------------------------------------------------------------------
# PURPOSE:
#   This section pregenerates all synthetic panel datasets used throughout the
#   Monte Carlo study. Each dataset corresponds to one scenario (1–3) and one
#   repetition (1..nsim). The files are stored to disk so that every method
#   can operate on identical input data for fair comparison.
#
# DESIGN:
#   • Three causal structures ("scenarios"):
#       1. No causal effect (confounding only)
#       2. X → Y causal effect
#       3. Bidirectional feedback X ↔ Y
#   • Five time points per subject (baseline + 4 follow-ups)
#   • N = 1000 subjects per dataset
#   • nsim = 2000 Monte Carlo repetitions per scenario
#   • Randomness controlled via a deterministic RNG log (rng_log.csv)
#
# OUTPUTS:
#   - /content/data/dgm{scenario}_rep{rep}.rds   → single dataset per replicate
#   - Each file contains a long-format data.frame with columns:
#       id, time, c (confounder), X (exposure), Y (outcome), scenario
#
# REPRODUCIBILITY:
#   Each (scenario, repetition) pair has a unique seed taken from rng_log.
#   This ensures full reproducibility and cross-method comparability.
#
# RUNTIME:
#   With N = 1000 and nsim = 2000 per scenario, expect ≈15–25 minutes total on
#   modern CPU (depends on host). Use smaller nsim for test runs.
###############################################################################

message("Starting dataset generation. nsim = ", nsim)

for (s in 1:3) {
  set.seed(master_seeds[s])  # scenario-level reproducibility

  for (r in 1:nsim) {
    # Retrieve the deterministic seed for this repetition
    seed_r <- rng_log %>%
      filter(scenario == s, rep_id == r) %>%
      pull(seed)
    set.seed(seed_r)

    # Simulate dataset (baseline + 4 follow-ups = 5 total time points)
    dat <- simulate_panel(scenario = s, N = 1000, T_obs = 4)

    # Save to disk in consistent structure
    saveRDS(dat, file = file.path(save_dir, sprintf("dgm%i_rep%i.rds", s, r)))

    # Lightweight progress indicator
    if (r %% 100 == 0) {
      message("Scenario ", s, " — ", r, " / ", nsim, " datasets saved.")
    }
  }
}

message("Dataset generation completed. Files stored in ", save_dir)