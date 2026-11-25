#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(GenomicRanges)
  library(IRanges)
})

set.seed(20260928)

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
input_dir <- file.path(dirname(script_dir), "peak_annotation")
data_dir <- file.path(script_dir, "data")
report_dir <- file.path(script_dir, "reports")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

n_group <- 50L
n_boot <- 2000L
n_sim <- 1000L

input_files <- c(
  EBS = "EBS-FLAG_peak_tair12_annotated.xlsx",
  WT43 = "WT-H3K4me3_peak_tair12_annotated.xlsx",
  SPO11 = "SPO11-1-FLAG_peak_tair12_annotated.xlsx",
  SPO11_EBS = "ebs-SPO11-1-FLAG_peak_tair12_annotated.xlsx"
)

read_unique_peaks <- function(filename, keep_signal = FALSE) {
  x <- read_excel(file.path(input_dir, filename), sheet = "all_annotated") %>%
    filter(type == "euchromatin") %>%
    transmute(
      chr = as.character(seqnames),
      start = as.integer(start),
      end = as.integer(end),
      signal = if (keep_signal) as.numeric(V9) else NA_real_,
      geneId = as.character(geneId)
    ) %>%
    filter(!is.na(chr), !is.na(start), !is.na(end), start >= 1L, end >= start) %>%
    arrange(chr, start, end)

  x %>% distinct(chr, start, end, .keep_all = TRUE)
}

ebs <- read_unique_peaks(input_files[["EBS"]], keep_signal = TRUE) %>%
  filter(is.finite(signal)) %>%
  arrange(signal, chr, start, end) %>%
  mutate(
    signal_group = ntile(row_number(), n_group),
    percent_rank = 100 * signal_group / n_group
  )

wt43 <- read_unique_peaks(input_files[["WT43"]])
spo11 <- read_unique_peaks(input_files[["SPO11"]])
spo11_ebs <- read_unique_peaks(input_files[["SPO11_EBS"]])

to_gr <- function(x) {
  GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end))
}

ebs_gr <- to_gr(ebs)
wt43_gr <- to_gr(wt43)
spo11_gr <- to_gr(spo11)
spo11_ebs_gr <- to_gr(spo11_ebs)

# Use the observed peak extent on each chromosome.
all_peaks <- bind_rows(
  ebs %>% select(chr, start, end),
  wt43 %>% select(chr, start, end),
  spo11 %>% select(chr, start, end),
  spo11_ebs %>% select(chr, start, end)
)

chrom_bounds <- all_peaks %>%
  group_by(chr) %>%
  summarise(lo = min(start), hi = max(end), .groups = "drop")

bound_lo <- setNames(chrom_bounds$lo, chrom_bounds$chr)
bound_hi <- setNames(chrom_bounds$hi, chrom_bounds$chr)

# Split wrapped peaks to preserve their total width.
circular_shift_gr <- function(x) {
  out_chr <- character(0)
  out_start <- integer(0)
  out_end <- integer(0)

  for (ch in unique(x$chr)) {
    z <- x[x$chr == ch, , drop = FALSE]
    lo <- unname(bound_lo[[ch]])
    hi <- unname(bound_hi[[ch]])
    domain_width <- hi - lo + 1L
    offset <- sample.int(domain_width, 1L) - 1L
    widths <- z$end - z$start + 1L
    shifted_start <- lo + ((z$start - lo + offset) %% domain_width)
    shifted_end <- shifted_start + widths - 1L

    nowrap <- shifted_end <= hi
    if (any(nowrap)) {
      out_chr <- c(out_chr, rep(ch, sum(nowrap)))
      out_start <- c(out_start, shifted_start[nowrap])
      out_end <- c(out_end, shifted_end[nowrap])
    }

    if (any(!nowrap)) {
      wrap_n <- sum(!nowrap)
      overflow <- shifted_end[!nowrap] - hi
      out_chr <- c(out_chr, rep(ch, 2L * wrap_n))
      out_start <- c(out_start, shifted_start[!nowrap], rep(lo, wrap_n))
      out_end <- c(out_end, rep(hi, wrap_n), lo + overflow - 1L)
    }
  }

  GRanges(seqnames = out_chr, ranges = IRanges(start = out_start, end = out_end))
}

overlap_flag <- function(query, subject) countOverlaps(query, subject, ignore.strand = TRUE) > 0L

obs_wt43 <- overlap_flag(ebs_gr, wt43_gr)
obs_spo11 <- overlap_flag(ebs_gr, spo11_gr)
obs_spo11_ebs <- overlap_flag(ebs_gr, spo11_ebs_gr)

metric_def <- list(
  wt43_ratio = list(num = obs_wt43, den = rep(TRUE, nrow(ebs))),
  spo11_ratio = list(num = obs_spo11, den = rep(TRUE, nrow(ebs))),
  spo11_ebs_ratio = list(num = obs_spo11_ebs, den = rep(TRUE, nrow(ebs))),
  spo11_wt43_ratio = list(num = obs_wt43 & obs_spo11, den = obs_wt43),
  spo11_ebs_wt43_ratio = list(num = obs_wt43 & obs_spo11_ebs, den = obs_wt43)
)

ratio_or_na <- function(num, den) {
  n_den <- sum(den)
  if (n_den == 0L) return(NA_real_)
  sum(num & den) / n_den
}

bootstrap_ci <- function(num, den, groups, n_rep = n_boot) {
  ans <- matrix(NA_real_, nrow = n_group, ncol = 2L)
  for (g in seq_len(n_group)) {
    ids <- which(groups == g)
    vals <- replicate(n_rep, {
      s <- sample(ids, length(ids), replace = TRUE)
      ratio_or_na(num[s], den[s])
    })
    if (any(is.finite(vals))) {
      ans[g, ] <- quantile(vals, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    }
  }
  ans
}

message("Calculating observed values and bootstrap intervals...")
observed_tables <- lapply(names(metric_def), function(metric) {
  d <- metric_def[[metric]]
  ci <- bootstrap_ci(d$num, d$den, ebs$signal_group)
  tibble(
    metric = metric,
    signal_group = seq_len(n_group),
    percent_rank = 100 * seq_len(n_group) / n_group,
    mean_signal = vapply(seq_len(n_group), function(g) mean(ebs$signal[ebs$signal_group == g]), numeric(1)),
    n_total = vapply(seq_len(n_group), function(g) sum(ebs$signal_group == g), integer(1)),
    n_denominator = vapply(seq_len(n_group), function(g) sum(d$den[ebs$signal_group == g]), integer(1)),
    observed = vapply(seq_len(n_group), function(g) {
      ids <- ebs$signal_group == g
      ratio_or_na(d$num[ids], d$den[ids])
    }, numeric(1)),
    obs_ci_low = ci[, 1],
    obs_ci_high = ci[, 2]
  )
})
names(observed_tables) <- names(metric_def)

random_mats <- lapply(names(metric_def), function(x) matrix(NA_real_, nrow = n_group, ncol = n_sim))
names(random_mats) <- names(metric_def)

message("Calculating chromosome-wise circular-shift controls...")
for (i in seq_len(n_sim)) {
  rand_wt43 <- overlap_flag(ebs_gr, circular_shift_gr(wt43))
  rand_spo11 <- overlap_flag(ebs_gr, circular_shift_gr(spo11))
  rand_spo11_ebs <- overlap_flag(ebs_gr, circular_shift_gr(spo11_ebs))

  for (g in seq_len(n_group)) {
    ids <- ebs$signal_group == g
    random_mats$wt43_ratio[g, i] <- mean(rand_wt43[ids])
    random_mats$spo11_ratio[g, i] <- mean(rand_spo11[ids])
    random_mats$spo11_ebs_ratio[g, i] <- mean(rand_spo11_ebs[ids])

    # Keep the observed EBS–WT43 denominator fixed.
    den_ids <- ids & obs_wt43
    random_mats$spo11_wt43_ratio[g, i] <- if (any(den_ids)) mean(rand_spo11[den_ids]) else NA_real_
    random_mats$spo11_ebs_wt43_ratio[g, i] <- if (any(den_ids)) mean(rand_spo11_ebs[den_ids]) else NA_real_
  }

  if (i %% 100L == 0L) message("  permutation ", i, "/", n_sim)
}

summarise_random <- function(m) {
  tibble(
    random_mean = apply(m, 1L, mean, na.rm = TRUE),
    random_ci_low = apply(m, 1L, quantile, probs = 0.025, na.rm = TRUE, names = FALSE),
    random_ci_high = apply(m, 1L, quantile, probs = 0.975, na.rm = TRUE, names = FALSE),
    random_sd = apply(m, 1L, sd, na.rm = TRUE)
  )
}

results <- bind_rows(lapply(names(metric_def), function(metric) {
  bind_cols(observed_tables[[metric]], summarise_random(random_mats[[metric]]))
}))

write.csv(results, file.path(data_dir, "ratio_metrics_results_improved.csv"), row.names = FALSE)

safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || sd(y[ok]) == 0) return(c(estimate = NA_real_, p = NA_real_))
  z <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  c(estimate = unname(z$estimate), p = z$p.value)
}

correlation_summary <- bind_rows(lapply(names(metric_def), function(metric) {
  z <- results %>% filter(.data$metric == .env$metric)
  obs_s <- safe_cor(z$percent_rank, z$observed, "spearman")
  obs_p <- safe_cor(z$percent_rank, z$observed, "pearson")
  rnd_s <- safe_cor(z$percent_rank, z$random_mean, "spearman")
  rnd_p <- safe_cor(z$percent_rank, z$random_mean, "pearson")
  tibble(
    metric = metric,
    obs_spearman_rho = obs_s[["estimate"]], obs_spearman_p = obs_s[["p"]],
    obs_pearson_r = obs_p[["estimate"]], obs_pearson_p = obs_p[["p"]],
    random_spearman_rho = rnd_s[["estimate"]], random_spearman_p = rnd_s[["p"]],
    random_pearson_r = rnd_p[["estimate"]], random_pearson_p = rnd_p[["p"]]
  )
}))
write.csv(correlation_summary, file.path(data_dir, "correlation_summary_improved.csv"), row.names = FALSE)

overall_summary <- bind_rows(lapply(names(metric_def), function(metric) {
  d <- metric_def[[metric]]
  rand_overall <- apply(random_mats[[metric]], 2L, function(v) {
    if (metric %in% c("wt43_ratio", "spo11_ratio", "spo11_ebs_ratio")) {
      weighted.mean(v, observed_tables[[metric]]$n_total, na.rm = TRUE)
    } else {
      weighted.mean(v, observed_tables[[metric]]$n_denominator, na.rm = TRUE)
    }
  })
  tibble(
    metric = metric,
    observed_overall = ratio_or_na(d$num, d$den),
    random_mean_overall = mean(rand_overall, na.rm = TRUE),
    random_ci_low_overall = quantile(rand_overall, 0.025, na.rm = TRUE, names = FALSE),
    random_ci_high_overall = quantile(rand_overall, 0.975, na.rm = TRUE, names = FALSE),
    fold_over_random = observed_overall / random_mean_overall,
    excess_over_random = observed_overall - random_mean_overall
  )
}))
write.csv(overall_summary, file.path(data_dir, "overall_summary_improved.csv"), row.names = FALSE)

old_geneid_summary <- tibble(
  metric = names(metric_def),
  old_geneid_ratio = c(
    mean(ebs$geneId %in% unique(wt43$geneId)),
    mean(ebs$geneId %in% unique(spo11$geneId)),
    mean(ebs$geneId %in% unique(spo11_ebs$geneId)),
    sum((ebs$geneId %in% unique(wt43$geneId)) & (ebs$geneId %in% unique(spo11$geneId))) / sum(ebs$geneId %in% unique(wt43$geneId)),
    sum((ebs$geneId %in% unique(wt43$geneId)) & (ebs$geneId %in% unique(spo11_ebs$geneId))) / sum(ebs$geneId %in% unique(wt43$geneId))
  )
) %>%
  left_join(overall_summary %>% select(metric, direct_peak_ratio = observed_overall), by = "metric") %>%
  mutate(geneid_minus_direct = old_geneid_ratio - direct_peak_ratio)
write.csv(old_geneid_summary, file.path(data_dir, "old_geneid_vs_direct_peak_overlap.csv"), row.names = FALSE)

fmt <- function(x, digits = 4L) formatC(x, digits = digits, format = "f")
summary_lines <- c(
  "Fig. S7c-style peak-overlap analysis report",
  "============================================",
  "",
  paste0("Run date: ", Sys.Date()),
  paste0("Random seed: 20260928"),
  paste0("Signal bins: ", n_group),
  paste0("Observed bootstrap replicates per bin: ", n_boot),
  paste0("Circular-shift permutations: ", n_sim),
  "",
  "1. Problem found in the original script",
  "---------------------------------------",
  "The original observed statistic was defined by shared annotated geneId values,",
  "whereas its randomized statistic moved peaks in genomic space and then assigned",
  "genes. These are different statistical objects. A peak can be assigned to a target",
  "gene without directly overlapping a target peak, and long/dense gene regions raise",
  "the randomized baseline. The original approach therefore overstates both observed",
  "and background overlap relative to the interval-overlap quantity shown in Fig. S7c.",
  "",
  "2. Revised processing",
  "---------------------",
  "- Read the all_annotated sheet from each supplied XLSX file.",
  "- Retained rows labelled type == euchromatin, matching the original script.",
  "- Removed any duplicate genomic intervals; the peak interval is the unit of analysis.",
  "- Ranked EBS peaks by V9 (-log10 q-value) and divided them into 50 nearly equal bins.",
  "- Calculated all overlaps with direct chromosome/start/end interval intersections.",
  "- Calculated 95% confidence intervals for observed values with 2,000 bootstrap samples.",
  "- Constructed the null by shifting every target peak set around its own chromosome",
  "  using one random offset per chromosome. Wrapped intervals were split at the boundary.",
  "  This preserves peak number, peak widths, chromosome assignment and peak clustering.",
  "- Used the randomized mean and the 2.5th/97.5th percentiles of 1,000 permutations.",
  "- For conditional ratios, retained the real EBS-WT43 denominator and randomized only",
  "  SPO11 or SPO11_EBS. This tests extra co-localization conditional on WT43 overlap.",
  "",
  "3. Metric definitions",
  "---------------------",
  "wt43_ratio = EBS peaks overlapping at least one WT43 peak / all EBS peaks.",
  "spo11_ratio = EBS peaks overlapping at least one SPO11 peak / all EBS peaks.",
  "spo11_ebs_ratio = EBS peaks overlapping at least one SPO11_EBS peak / all EBS peaks.",
  "spo11_wt43_ratio = EBS peaks overlapping both WT43 and SPO11 / EBS peaks overlapping WT43.",
  "spo11_ebs_wt43_ratio = EBS peaks overlapping both WT43 and SPO11_EBS / EBS peaks overlapping WT43.",
  "",
  "4. Overall results",
  "------------------"
)

for (i in seq_len(nrow(overall_summary))) {
  z <- overall_summary[i, ]
  summary_lines <- c(
    summary_lines,
    paste0(
      z$metric, ": observed=", fmt(z$observed_overall),
      "; random mean=", fmt(z$random_mean_overall),
      " (95% interval ", fmt(z$random_ci_low_overall), "-", fmt(z$random_ci_high_overall), ")",
      "; fold/random=", fmt(z$fold_over_random, 2L),
      "; observed-random=", fmt(z$excess_over_random)
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "Interpretation of whether the random control should be near zero:",
  "A valid random control is not forced to zero. It equals the chance overlap implied by",
  "target-peak density, widths and chromosome distribution. Sparse SPO11_EBS peaks should",
  "produce a control close to zero; the much denser WT43 peak set is expected to have a",
  "higher non-zero chance-overlap baseline. Subtracting or rescaling the null purely to make",
  "it visually zero would bias the comparison. The excess_over_random column is provided if",
  "a zero-centered effect size is needed.",
  "",
  "5. Changes from the original R script",
  "-------------------------------------",
  "- Replaced geneId membership with direct peak-coordinate overlap.",
  "- Replaced randomizeRegions over gene ranges with chromosome-wise circular shifts.",
  "- Increased random simulations from 50 to 1,000.",
  "- Replaced mean +/- SD random bars with empirical 95% permutation intervals.",
  "- Enabled all five requested metrics.",
  "- Made signal bins robust with ntile, avoiding duplicate-quantile break failures.",
  "- Used fixed seeds and wrote complete numeric outputs for reproducibility.",
  "",
  "6. Limitations",
  "--------------",
  "No chromosome-size/blacklist/euchromatin-mask file was supplied. Circular-shift bounds",
  "were therefore inferred from the combined coordinate extent of the supplied peak files.",
  "This is preferable to inventing a genome build, but a validated callable-region mask can",
  "be substituted later if available. The random baseline may then change slightly.",
  "",
  "7. Output files",
  "---------------",
  "- ratio_metrics_results_improved.csv (per-bin values and intervals).",
  "- overall_summary_improved.csv (overall observed/random comparisons).",
  "- correlation_summary_improved.csv (rank correlations).",
  "- old_geneid_vs_direct_peak_overlap.csv (diagnostic comparison).",
  "- analysis_report.txt (analysis summary)."
)

writeLines(summary_lines, file.path(report_dir, "analysis_report.txt"), useBytes = TRUE)

message("Analysis complete. Outputs written to: ", script_dir)
