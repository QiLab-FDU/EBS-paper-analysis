#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

script_dir <- local({
  for (i in rev(seq_len(sys.nframe()))) {
    if (identical(sys.function(i), base::source)) {
      source_frame <- sys.frame(i)
      path <- source_frame$ofile
      if (is.character(path) && length(path) == 1L) {
        if (!is.null(source_frame$owd)) path <- basename(path)
        return(dirname(normalizePath(path, mustWork = TRUE)))
      }
    }
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(file_arg)) {
    stop("Open the script and use Source, or run it with Rscript.", call. = FALSE)
  }
  path <- sub("^--file=", "", file_arg[[1L]])
  path <- gsub("~+~", " ", path, fixed = TRUE)  # Decode spaces in --file.
  dirname(normalizePath(path, mustWork = TRUE))
})
data_dir <- file.path(script_dir, "data")

results <- read.csv(
  file.path(data_dir, "ratio_metrics_results_improved.csv"),
  check.names = FALSE
)

overall_summary <- read.csv(
  file.path(data_dir, "overall_summary_improved.csv"),
  check.names = FALSE
)

# Load correlations for inspection.
correlation_summary <- read.csv(
  file.path(data_dir, "correlation_summary_improved.csv"),
  check.names = FALSE
)

required_result_columns <- c(
  "metric", "percent_rank", "observed", "obs_ci_low", "obs_ci_high",
  "random_mean", "random_ci_low", "random_ci_high"
)

if (!all(required_result_columns %in% names(results))) {
  stop(
    "Missing columns in ratio_metrics_results_improved.csv: ",
    paste(
      setdiff(required_result_columns, names(results)),
      collapse = ", "
    )
  )
}

metric_order <- c(
  "wt43_ratio",
  "spo11_ratio",
  "spo11_ebs_ratio",
  "spo11_wt43_ratio",
  "spo11_ebs_wt43_ratio"
)

metric_titles <- c(
  wt43_ratio = "EBS peaks overlapping WT43",
  spo11_ratio = "EBS peaks overlapping SPO11",
  spo11_ebs_ratio = "EBS peaks overlapping SPO11-EBS",
  spo11_wt43_ratio = "SPO11 overlap among EBS-WT43 peaks",
  spo11_ebs_wt43_ratio = "SPO11-EBS overlap among EBS-WT43 peaks"
)

observed_color <- "#E64B35"
random_color <- "#2C7FB8"

make_plot <- function(metric_name) {
  z <- results[results$metric == metric_name, , drop = FALSE]
  overall <- overall_summary[
    overall_summary$metric == metric_name, , drop = FALSE
  ]

  if (nrow(z) == 0L || nrow(overall) != 1L) {
    stop("Missing or duplicated data for metric: ", metric_name)
  }

  y_max <- max(z$obs_ci_high, z$random_ci_high, na.rm = TRUE)
  y_upper <- min(1, max(0.05, y_max * 1.08))

  ggplot(z, aes(x = percent_rank)) +
    geom_hline(
      yintercept = overall$observed_overall,
      color = observed_color,
      linewidth = 0.45,
      linetype = "dashed",
      alpha = 0.65
    ) +
    geom_hline(
      yintercept = overall$random_mean_overall,
      color = random_color,
      linewidth = 0.45,
      linetype = "dashed",
      alpha = 0.70
    ) +
    geom_errorbar(
      aes(
        ymin = obs_ci_low,
        ymax = obs_ci_high,
        color = "Observed"
      ),
      width = 0,
      linewidth = 0.45,
      alpha = 0.90
    ) +
    geom_point(
      aes(y = observed, color = "Observed"),
      size = 1.45
    ) +
    geom_errorbar(
      aes(
        ymin = random_ci_low,
        ymax = random_ci_high,
        color = "Randomized control"
      ),
      width = 0,
      linewidth = 0.42,
      alpha = 0.85
    ) +
    geom_point(
      aes(y = random_mean, color = "Randomized control"),
      size = 1.30
    ) +
    scale_color_manual(values = c(
      "Observed" = observed_color,
      "Randomized control" = random_color
    )) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 25),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 1.05),
      expand = expansion(mult = c(0.025, 0.04))
    ) +
    labs(
      title = unname(metric_titles[[metric_name]]),
      x = expression("EBS peak percentile rank (-log"[10] * " q-value)"),
      y = metric_name,
      color = NULL
    ) +
    theme_classic(base_size = 10.5, base_family = "sans") +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0),
      axis.title = element_text(size = 9.5, color = "black"),
      axis.text = element_text(size = 8.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.line = element_line(linewidth = 0.40, color = "black"),
      panel.border = element_rect(
        color = "black", fill = NA, linewidth = 0.45
      ),
      legend.position = "top",
      legend.text = element_text(size = 8.3),
      legend.key.width = grid::unit(0.85, "cm"),
      plot.margin = margin(6, 8, 6, 6)
    ) +
    guides(color = guide_legend(
      override.aes = list(linewidth = 0.7, size = 2)
    ))
}

p_wt43_ratio <- make_plot("wt43_ratio")
p_spo11_ratio <- make_plot("spo11_ratio")
p_spo11_ebs_ratio <- make_plot("spo11_ebs_ratio")
p_spo11_wt43_ratio <- make_plot("spo11_wt43_ratio")
p_spo11_ebs_wt43_ratio <- make_plot("spo11_ebs_wt43_ratio")

plots <- list(
  wt43_ratio = p_wt43_ratio,
  spo11_ratio = p_spo11_ratio,
  spo11_ebs_ratio = p_spo11_ebs_ratio,
  spo11_wt43_ratio = p_spo11_wt43_ratio,
  spo11_ebs_wt43_ratio = p_spo11_ebs_wt43_ratio
)

if (interactive()) {
  print(p_wt43_ratio)
  print(p_spo11_ratio)
  print(p_spo11_ebs_ratio)
  print(p_spo11_wt43_ratio)
  print(p_spo11_ebs_wt43_ratio)
  message(
    "Existing results loaded. Five separate plots were sent to ",
    "the RStudio Plots pane; no randomization was rerun."
  )
} else {
  message(
    "Existing results loaded. Non-interactive mode: plot objects were created ",
    "but not printed, so no Rplots.pdf file was produced."
  )
}

