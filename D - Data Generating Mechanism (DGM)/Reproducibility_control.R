###############################################################################
# GLOBAL REPRODUCIBILITY CONTROL
# -----------------------------------------------------------------------------
# Creates a central RNG table and a workspace folder for simulated data.
# All simulation steps (SEM, DYNAMITE, etc.) will reference this same folder.
###############################################################################

library(tidyverse)
library(data.table)

# ---- 1. Define directories ---------------------------------------------------
# Adjust here if needed for your own working environment
save_dir <- "/content/data"     # folder to store simulated datasets
meta_dir <- "/content/meta"     # folder to store logs (e.g., RNG table)

dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(meta_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Global seeds --------------------------------------------------------
master_seeds <- c(20251001, 20251002, 20251003)  # one per scenario
nsim <- 2000                                     # number of Monte Carlo repetitions

# ---- 3. Build RNG log table -------------------------------------------------
# Each repetition within a scenario gets a unique deterministic seed
rng_log <- tibble(
  scenario = rep(1:3, each = nsim),
  rep_id   = rep(1:nsim, times = 3),
  seed     = rep(master_seeds, each = nsim) + 10000 * rep(1:3, each = nsim) + 1:nsim
)

# ---- 4. Save RNG log --------------------------------------------------------
write_csv(rng_log, file.path(meta_dir, "rng_log.csv"))

# ---- 5. Confirm -------------------------------------------------------------
message("RNG table created at: ", file.path(meta_dir, "rng_log.csv"))
message("Data directory: ", save_dir)