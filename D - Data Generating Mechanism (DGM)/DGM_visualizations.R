###############################################################################
# VISUALS (ADEMP — DATA-GENERATING MECHANISM)
# -----------------------------------------------------------------------------
# Purpose:
#   Produce high-quality figures illustrating the data-generating mechanism
#   and temporal behavior of simulated variables (c, X, Y) across time points.
#
###############################################################################

library(tidyverse)

# -----------------------------------------------------------------------------
# Shared ggplot theme for publication
# -----------------------------------------------------------------------------
theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11)
  )

scen_labs <- c(
  "Scenario 1: No causal effect",
  "Scenario 2: X → Y causal effect",
  "Scenario 3: Bidirectional feedback (X ↔ Y)"
)

# -----------------------------------------------------------------------------
# DAG schematic (schematic, clean, color-coded)
# -----------------------------------------------------------------------------
plot_dag_scenario_pub <- function(scenario_id, outfile) {
  png(outfile, width = 1200, height = 650, res = 160)
  par(mar = c(2,2,2,2))
  plot.new(); plot.window(xlim = c(0, 5), ylim = c(0, 6))

  node <- function(x, y, label, col) {
    symbols(x, y, circles = 0.28, inches = FALSE, add = TRUE, bg = col, fg = NA)
    text(x, y, label, cex = 0.9, font = 2)
  }
  arrow <- function(x0, y0, x1, y1, lty = 1, lwd = 2, col = "grey20") {
    arrows(x0, y0, x1, y1, length = 0.08, lty = lty, lwd = lwd, col = col)
  }

  times <- 1:4  # baseline + 3 follow-ups (4 time points total)
  for (t in times) {
    node(t - 0.5, 5, paste0("c", t), "#F6C142")
    node(t, 3, paste0("X", t), "#5AB4E6")
    node(t, 1, paste0("Y", t), "#00A37A")
  }
  # temporal confounder evolution
  for (t in 2:4) arrow((t-1)-0.5+0.25, 5, (t-0.5)-0.25, 5, lty = 2, col="grey40")
  # confounder → X, Y
  for (t in times) {
    arrow(t-0.5, 4.75, t, 3.25)
    arrow(t-0.5, 4.75, t, 1.25)
  }
  # AR(1) within X and Y
  for (t in 2:4) {
    arrow(t-1+0.25, 3, t-0.25, 3, lty = 2, col="grey40")
    arrow(t-1+0.25, 1, t-0.25, 1, lty = 2, col="grey40")
  }
  # scenario-specific causal arrows
  if (scenario_id == 2) {
    for (t in 1:3) arrow(t+0.2, 2.9, t+1-0.2, 1.1, col = "#1F78B4")
    title("Scenario 2 — X → Y causal effect", font.main = 2)
  } else if (scenario_id == 3) {
    for (t in 2:4) {
      arrow(t-1+0.2, 2.9, t-0.2, 1.1, col = "#1F78B4")
      arrow(t-1+0.2, 1.1, t-0.2, 2.9, col = "#E45756")
    }
    title("Scenario 3 — Bidirectional feedback X ↔ Y", font.main = 2)
  } else {
    title("Scenario 1 — Confounding only", font.main = 2)
  }

  dev.off()
  message("Saved DAG schematic: ", outfile)
}

# -----------------------------------------------------------------------------
# Generate and save figures
# -----------------------------------------------------------------------------
set.seed(master_seeds[1])
vis_s1 <- simulate_panel(1, N = 500, T_obs = 4)
vis_s2 <- simulate_panel(2, N = 500, T_obs = 4)
vis_s3 <- simulate_panel(3, N = 500, T_obs = 4)

plot_dag_scenario_pub(1, file.path(meta_dir, "FigDAG_S1.png"))
plot_dag_scenario_pub(2, file.path(meta_dir, "FigDAG_S2.png"))
plot_dag_scenario_pub(3, file.path(meta_dir, "FigDAG_S3.png"))

# -----------------------------------------------------------------------------
# Sanity check
# -----------------------------------------------------------------------------
library(glue)

check_means_pretty <- function(dat, label) {
  msg <- dat |>
    group_by(time) |>
    summarise(
      mean_c = glue("{round(mean(c),2)} ({round(mean(c)-1.96*sd(c)/sqrt(n()),2)}–{round(mean(c)+1.96*sd(c)/sqrt(n()),2)})"),
      mean_X = glue("{round(mean(X),2)} ({round(mean(X)-1.96*sd(X)/sqrt(n()),2)}–{round(mean(X)+1.96*sd(X)/sqrt(n()),2)})"),
      mean_Y = glue("{round(mean(Y),2)} ({round(mean(Y)-1.96*sd(Y)/sqrt(n()),2)}–{round(mean(Y)+1.96*sd(Y)/sqrt(n()),2)})"),
      .groups = "drop"
    )
  cat("\n", label, "\n")
  print(msg)
}

# Example usage
set.seed(1)
tmp1 <- simulate_panel(1, N = 10000, T_obs = 4)
tmp2 <- simulate_panel(2, N = 10000, T_obs = 4)
tmp3 <- simulate_panel(3, N = 10000, T_obs = 4)

check_means_pretty(tmp1, "S1")
check_means_pretty(tmp2, "S2")
check_means_pretty(tmp3, "S3")

# =============================================================================
# VISUALIZATION SUITE FOR PANEL DATA SIMULATIONS
# =============================================================================
# Purpose: Create figures that reveal longitudinal structure,
#          causal relationships, and confounding in simulated panel data
# =============================================================================

# -----------------------------------------------------------------------------
# Shared ggplot theme
# -----------------------------------------------------------------------------
theme_pub <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Color palette (consistent across all plots)
col_c <- "#E69F00"  # Confounder
col_X <- "#0072B2"  # Exposure
col_Y <- "#009E73"  # Outcome

scen_labs <- c(
  "Scenario 1: No Causal Effect (Confounding Only)",
  "Scenario 2: Unidirectional X → Y",
  "Scenario 3: Bidirectional Feedback X ↔ Y"
)

# =============================================================================
# 1. SPAGHETTI PLOTS: Individual Trajectories
# =============================================================================
#' Plot Individual Trajectories (Spaghetti Plot)
#'
#' Shows trajectories for a random sample of individuals to reveal:
#' - Within-person persistence (autoregressive structure)
#' - Between-person variability
#' - Monotonic trends (e.g., confounder drift)
#'
#' @param data Data frame from simulate_panel() in long format
#' @param scenario_id Integer 1-3 indicating scenario
#' @param n_lines Number of individuals to plot (default 30)
#' @param alpha_lines Transparency for individual lines (default 0.4)
#' @return ggplot object (printed and returned invisibly)
plot_trajectories <- function(data, scenario_id, n_lines = 30,
                              alpha_lines = 0.4) {

  # Sample random individuals
  sampled_ids <- sample(unique(data$id), min(n_lines, length(unique(data$id))))

  plot_data <- data %>%
    filter(id %in% sampled_ids) %>%
    pivot_longer(cols = c(c, X, Y), names_to = "variable", values_to = "value") %>%
    mutate(
      variable = factor(variable,
                       levels = c("c", "X", "Y"),
                       labels = c("Confounder (c)", "Exposure (X)", "Outcome (Y)"))
    )

  # Calculate population means for reference
  pop_means <- data %>%
    pivot_longer(cols = c(c, X, Y), names_to = "variable", values_to = "value") %>%
    mutate(
      variable = factor(variable,
                       levels = c("c", "X", "Y"),
                       labels = c("Confounder (c)", "Exposure (X)", "Outcome (Y)"))
    ) %>%
    group_by(time, variable) %>%
    summarise(mean_val = mean(value), .groups = "drop")

  p <- ggplot() +
    # Individual trajectories (spaghetti)
    geom_line(data = plot_data,
              aes(x = time, y = value, group = id, color = variable),
              alpha = alpha_lines, linewidth = 0.5) +
    # Population mean (bold line)
    geom_line(data = pop_means,
              aes(x = time, y = mean_val, color = variable),
              linewidth = 1.5) +
    facet_wrap(~ variable, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("Confounder (c)" = col_c,
                                   "Exposure (X)" = col_X,
                                   "Outcome (Y)" = col_Y)) +
    labs(
      title = scen_labs[scenario_id],
      subtitle = paste("Individual trajectories (n =", n_lines,
                      ") with population mean (bold line)"),
      x = "Time",
      y = "Value"
    ) +
    theme_pub +
    theme(legend.position = "none")  # Legend redundant with facets

  print(p)
  invisible(p)
}


# =============================================================================
# 2. MEAN TRAJECTORIES WITH UNCERTAINTY
# =============================================================================
#' Plot Population Mean Trajectories with Confidence Bands
#'
#' Shows mean trajectories with ±1 SD ribbons. Useful for:
#' - Showing population-level trends (e.g., confounder drift)
#' - Comparing means across scenarios
#' - Publication-ready summary statistics
#'
#' @param data Data frame from simulate_panel() in long format
#' @param scenario_id Integer 1-3 indicating scenario
#' @return ggplot object (printed and returned invisibly)
plot_mean_trajectories <- function(data, scenario_id) {

  summary_data <- data %>%
    pivot_longer(cols = c(c, X, Y), names_to = "variable", values_to = "value") %>%
    mutate(
      variable = factor(variable,
                       levels = c("c", "X", "Y"),
                       labels = c("Confounder (c)", "Exposure (X)", "Outcome (Y)"))
    ) %>%
    group_by(time, variable) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      .groups = "drop"
    )

  p <- ggplot(summary_data, aes(x = time, y = mean, color = variable, fill = variable)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3, shape = 21, color = "black") +
    scale_color_manual(values = c("Confounder (c)" = col_c,
                                   "Exposure (X)" = col_X,
                                   "Outcome (Y)" = col_Y)) +
    scale_fill_manual(values = c("Confounder (c)" = col_c,
                                  "Exposure (X)" = col_X,
                                  "Outcome (Y)" = col_Y)) +
    labs(
      title = scen_labs[scenario_id],
      subtitle = "Population mean trajectories ± SD",
      x = "Time",
      y = "Mean Value"
    ) +
    theme_pub

  print(p)
  invisible(p)
}


# =============================================================================
# 3. SCATTERPLOT MATRIX (Baseline vs Final)
# =============================================================================
#' Plot Scatterplot Matrix for Baseline and Final Timepoints
#'
#' Shows pairwise relationships between c, X, Y at baseline (t=0) and
#' final follow-up. Reveals:
#' - Cross-sectional confounding structure
#' - Stability of relationships over time
#' - Bivariate distributions
#'
#' @param data Data frame from simulate_panel() in long format
#' @param scenario_id Integer 1-3 indicating scenario
#' @return ggplot object (printed and returned invisibly)
plot_scatter_matrix <- function(data, scenario_id) {

  max_time <- max(data$time)

  # Extract baseline and final timepoints
  scatter_data <- data %>%
    filter(time %in% c(0, max_time)) %>%
    mutate(timepoint = ifelse(time == 0, "Baseline (t=0)",
                              paste0("Final (t=", max_time, ")"))) %>%
    select(id, timepoint, c, X, Y) %>%
    pivot_longer(cols = c(c, X, Y), names_to = "variable", values_to = "value")

  # Create wide format with separate columns for c, X, Y
  wide_scatter <- scatter_data %>%
    pivot_wider(names_from = variable, values_from = value)

  # Create all pairwise combinations manually
  var_names <- c("c", "X", "Y")
  pairs_list <- list()

  for (vx in var_names) {
    for (vy in var_names) {
      if (vx != vy) {
        pairs_list[[length(pairs_list) + 1]] <- wide_scatter %>%
          select(id, timepoint, all_of(c(vx, vy))) %>%
          rename(val_x = all_of(vx), val_y = all_of(vy)) %>%
          mutate(var_x = vx, var_y = vy)
      }
    }
  }

  pairs_data <- bind_rows(pairs_list) %>%
    mutate(
      var_x = factor(var_x, levels = c("c", "X", "Y")),
      var_y = factor(var_y, levels = c("c", "X", "Y"))
    )

  p <- ggplot(pairs_data, aes(x = val_x, y = val_y)) +
    geom_point(alpha = 0.3, size = 1, color = "grey30") +
    geom_smooth(method = "lm", se = TRUE, color = "#D55E00", linewidth = 1) +
    facet_grid(var_y ~ var_x + timepoint, scales = "free",
               labeller = labeller(
                 var_x = c("c" = "Confounder", "X" = "Exposure", "Y" = "Outcome"),
                 var_y = c("c" = "Confounder", "X" = "Exposure", "Y" = "Outcome")
               )) +
    labs(
      title = scen_labs[scenario_id],
      subtitle = "Pairwise relationships at baseline and final follow-up",
      x = NULL,
      y = NULL
    ) +
    theme_pub +
    theme(strip.text = element_text(size = 9))

  print(p)
  invisible(p)
}


# =============================================================================
# 4. COMPOSITE FIGURE: All Visualizations for One Scenario
# =============================================================================
#' Create Comprehensive Multi-Panel Figure
#'
#' Combines trajectory plot and mean trajectories into a single
#' publication-ready figure using patchwork.
#'
#' @param data Data frame from simulate_panel() in long format
#' @param scenario_id Integer 1-3 indicating scenario
#'
#' @return Combined patchwork object (printed and returned invisibly)
plot_scenario_comprehensive <- function(data, scenario_id) {

  # Generate individual plots (without printing)
  p1 <- plot_trajectories(data, scenario_id, n_lines = 30)
  suppressMessages(print(p1))  # Print with trajectory plot

  p2 <- plot_mean_trajectories(data, scenario_id)
  suppressMessages(print(p2))  # Print with mean plot

  invisible(list(trajectories = p1, means = p2))
}


# =============================================================================
# USAGE WITH MASTER SEEDS
# =============================================================================
if (TRUE) {
  # Assuming simulate_panel() function is loaded
  # Define master seeds (one per scenario)
  master_seeds <- c(20251001, 20251002, 20251003)

  # Generate data for all scenarios
  set.seed(master_seeds[1])
  data_s1 <- simulate_panel(scenario = 1, N = 500, T_obs = 4)

  set.seed(master_seeds[2])
  data_s2 <- simulate_panel(scenario = 2, N = 500, T_obs = 4)

  set.seed(master_seeds[3])
  data_s3 <- simulate_panel(scenario = 3, N = 500, T_obs = 4)

  # Generate all plots for each scenario (auto-printed)
  plot_trajectories(data_s1, scenario_id = 1)
  plot_mean_trajectories(data_s1, scenario_id = 1)
  plot_scatter_matrix(data_s1, scenario_id = 1)

  plot_trajectories(data_s2, scenario_id = 2)
  plot_mean_trajectories(data_s2, scenario_id = 2)
  plot_scatter_matrix(data_s2, scenario_id = 2)

  plot_trajectories(data_s3, scenario_id = 3)
  plot_mean_trajectories(data_s3, scenario_id = 3)
  plot_scatter_matrix(data_s3, scenario_id = 3)

}

