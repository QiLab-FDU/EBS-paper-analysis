library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(ggh4x)
library(corrplot)

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

# SPO11-1-FLAG
spo11_flag_df <- read.table(
  file.path(script_dir, "SPO11-1-FLAG_clean_log2_nosm.bedgraph"),
  header = FALSE
)

colnames(spo11_flag_df) <- c(
  "chr_id", "start", "end", "log2_signal"
)

chr_info <- data.frame(
  chr_id = paste0("Chr", 1:5),
  length_bp = c(
    32637894,
    27782540,
    26149336,
    25768865,
    30142610
  )
)

bin_size_bp   <- 1000 * 1000
slide_step_bp <- 50 * 1000

# Retain truncated windows at chromosome ends.
bin_df <- lapply(seq_len(nrow(chr_info)), function(i) {
  
  chr <- chr_info$chr_id[i]
  chr_len <- chr_info$length_bp[i]
  
  starts <- seq(
    from = 0,
    to = chr_len - 1,
    by = slide_step_bp
  )
  
  data.frame(
    chr_id = chr,
    bin_start = starts,
    bin_end = pmin(starts + bin_size_bp, chr_len)
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    window_len = bin_end - bin_start,
    is_full_window = window_len == bin_size_bp,
    mid_Mb = (bin_start + bin_end) / 2 / 1e6
  )

overlap_df <- spo11_flag_df %>%
  inner_join(bin_df, by = "chr_id") %>%
  filter(start < bin_end, end > bin_start) %>%
  mutate(
    overlap_bp = pmin(end, bin_end) -
      pmax(start, bin_start)
  )

# Weight signals by the number of bases overlapping each window.
bin_signal <- overlap_df %>%
  group_by(chr_id, bin_start, bin_end) %>%
  summarise(
    weighted_mean = sum(
      log2_signal * overlap_bp,
      na.rm = TRUE
    ) / sum(overlap_bp),
    .groups = "drop"
  )

spo11_final <- bin_df %>%
  left_join(
    bin_signal,
    by = c("chr_id", "bin_start", "bin_end")
  ) %>%
  arrange(chr_id, bin_start) %>%
  group_by(chr_id) %>%
  mutate(
    mean_log2 = if_else(
      is_full_window,
      weighted_mean,
      NA_real_
    )
  ) %>%
  # Extend the last complete-window value to the chromosome end.
  fill(mean_log2, .direction = "down") %>%
  mutate(
    # Compute Z-score parameters from complete windows only.
    chr_mean = mean(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    chr_sd = sd(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    z_score = if_else(
      is.na(chr_sd) | chr_sd == 0,
      0,
      (mean_log2 - chr_mean) / chr_sd
    )
  ) %>%
  ungroup() %>%
  mutate(group = "SPO11-1-FLAG") %>%
  select(
    chr_id,
    bin_start,
    bin_end,
    mid_Mb,
    window_len,
    is_full_window,
    mean_log2,
    z_score,
    group
  )


# EBS-FLAG
ebs_flag_df = read.table(file.path(script_dir, "EBS-FLAG_clean_log2_nosm.bedgraph"), header = F)

colnames(ebs_flag_df) <- c(
  "chr_id", "start", "end", "log2_signal"
)

chr_info <- data.frame(
  chr_id = paste0("Chr", 1:5),
  length_bp = c(
    32637894,
    27782540,
    26149336,
    25768865,
    30142610
  )
)

bin_size_bp   <- 1000 * 1000
slide_step_bp <- 50 * 1000

bin_df <- lapply(seq_len(nrow(chr_info)), function(i) {
  
  chr <- chr_info$chr_id[i]
  chr_len <- chr_info$length_bp[i]
  
  starts <- seq(
    from = 0,
    to = chr_len - 1,
    by = slide_step_bp
  )
  
  data.frame(
    chr_id = chr,
    bin_start = starts,
    bin_end = pmin(starts + bin_size_bp, chr_len)
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    window_len = bin_end - bin_start,
    is_full_window = window_len == bin_size_bp,
    mid_Mb = (bin_start + bin_end) / 2 / 1e6
  )

overlap_df <- ebs_flag_df %>%
  inner_join(bin_df, by = "chr_id") %>%
  filter(start < bin_end, end > bin_start) %>%
  mutate(
    overlap_bp = pmin(end, bin_end) -
      pmax(start, bin_start)
  )

bin_signal <- overlap_df %>%
  group_by(chr_id, bin_start, bin_end) %>%
  summarise(
    weighted_mean = sum(
      log2_signal * overlap_bp,
      na.rm = TRUE
    ) / sum(overlap_bp),
    .groups = "drop"
  )

ebs_final <- bin_df %>%
  left_join(
    bin_signal,
    by = c("chr_id", "bin_start", "bin_end")
  ) %>%
  arrange(chr_id, bin_start) %>%
  group_by(chr_id) %>%
  mutate(
    mean_log2 = if_else(
      is_full_window,
      weighted_mean,
      NA_real_
    )
  ) %>%
  fill(mean_log2, .direction = "down") %>%
  mutate(
    chr_mean = mean(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    chr_sd = sd(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    z_score = if_else(
      is.na(chr_sd) | chr_sd == 0,
      0,
      (mean_log2 - chr_mean) / chr_sd
    )
  ) %>%
  ungroup() %>%
  mutate(group = "EBS-FLAG") %>%
  select(
    chr_id,
    bin_start,
    bin_end,
    mid_Mb,
    window_len,
    is_full_window,
    mean_log2,
    z_score,
    group
  )


# ebs-SPO11-1-FLAG
spo11_ebs_flag_34_df <- read.table(file.path(script_dir, "ebs-SPO11-1-FLAG_clean_log2_nosm.bedgraph"), header = F)

colnames(spo11_ebs_flag_34_df) <- c(
  "chr_id", "start", "end", "log2_signal"
)

chr_info <- data.frame(
  chr_id = paste0("Chr", 1:5),
  length_bp = c(
    32637894,
    27782540,
    26149336,
    25768865,
    30142610
  )
)

bin_size_bp   <- 1000 * 1000
slide_step_bp <- 50 * 1000

bin_df <- lapply(seq_len(nrow(chr_info)), function(i) {
  
  chr <- chr_info$chr_id[i]
  chr_len <- chr_info$length_bp[i]
  
  starts <- seq(
    from = 0,
    to = chr_len - 1,
    by = slide_step_bp
  )
  
  data.frame(
    chr_id = chr,
    bin_start = starts,
    bin_end = pmin(starts + bin_size_bp, chr_len)
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    window_len = bin_end - bin_start,
    is_full_window = window_len == bin_size_bp,
    mid_Mb = (bin_start + bin_end) / 2 / 1e6
  )

overlap_df <- spo11_ebs_flag_34_df %>%
  inner_join(bin_df, by = "chr_id") %>%
  filter(start < bin_end, end > bin_start) %>%
  mutate(
    overlap_bp = pmin(end, bin_end) -
      pmax(start, bin_start)
  )

bin_signal <- overlap_df %>%
  group_by(chr_id, bin_start, bin_end) %>%
  summarise(
    weighted_mean = sum(
      log2_signal * overlap_bp,
      na.rm = TRUE
    ) / sum(overlap_bp),
    .groups = "drop"
  )

ebs_spo11_final <- bin_df %>%
  left_join(
    bin_signal,
    by = c("chr_id", "bin_start", "bin_end")
  ) %>%
  arrange(chr_id, bin_start) %>%
  group_by(chr_id) %>%
  mutate(
    mean_log2 = if_else(
      is_full_window,
      weighted_mean,
      NA_real_
    )
  ) %>%
  fill(mean_log2, .direction = "down") %>%
  mutate(
    chr_mean = mean(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    chr_sd = sd(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    z_score = if_else(
      is.na(chr_sd) | chr_sd == 0,
      0,
      (mean_log2 - chr_mean) / chr_sd
    )
  ) %>%
  ungroup() %>%
  mutate(group = "ebs-SPO11-1-FLAG") %>%
  select(
    chr_id,
    bin_start,
    bin_end,
    mid_Mb,
    window_len,
    is_full_window,
    mean_log2,
    z_score,
    group
  )


# WT-H3K4me3
h43_wt_df = read.table(file.path(script_dir, "WT-H3K4me3_clean_log2_nosm.bedgraph"), header = F)
colnames(h43_wt_df) <- c(
  "chr_id", "start", "end", "log2_signal"
)

chr_info <- data.frame(
  chr_id = paste0("Chr", 1:5),
  length_bp = c(
    32637894,
    27782540,
    26149336,
    25768865,
    30142610
  )
)

bin_size_bp   <- 1000 * 1000
slide_step_bp <- 50 * 1000

bin_df <- lapply(seq_len(nrow(chr_info)), function(i) {
  
  chr <- chr_info$chr_id[i]
  chr_len <- chr_info$length_bp[i]
  
  starts <- seq(
    from = 0,
    to = chr_len - 1,
    by = slide_step_bp
  )
  
  data.frame(
    chr_id = chr,
    bin_start = starts,
    bin_end = pmin(starts + bin_size_bp, chr_len)
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    window_len = bin_end - bin_start,
    is_full_window = window_len == bin_size_bp,
    mid_Mb = (bin_start + bin_end) / 2 / 1e6
  )

overlap_df <- h43_wt_df %>%
  inner_join(bin_df, by = "chr_id") %>%
  filter(start < bin_end, end > bin_start) %>%
  mutate(
    overlap_bp = pmin(end, bin_end) -
      pmax(start, bin_start)
  )

bin_signal <- overlap_df %>%
  group_by(chr_id, bin_start, bin_end) %>%
  summarise(
    weighted_mean = sum(
      log2_signal * overlap_bp,
      na.rm = TRUE
    ) / sum(overlap_bp),
    .groups = "drop"
  )

h43_wt_final <- bin_df %>%
  left_join(
    bin_signal,
    by = c("chr_id", "bin_start", "bin_end")
  ) %>%
  arrange(chr_id, bin_start) %>%
  group_by(chr_id) %>%
  mutate(
    mean_log2 = if_else(
      is_full_window,
      weighted_mean,
      NA_real_
    )
  ) %>%
  fill(mean_log2, .direction = "down") %>%
  mutate(
    chr_mean = mean(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    chr_sd = sd(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    z_score = if_else(
      is.na(chr_sd) | chr_sd == 0,
      0,
      (mean_log2 - chr_mean) / chr_sd
    )
  ) %>%
  ungroup() %>%
  mutate(group = "WT-H3K4me3") %>%
  select(
    chr_id,
    bin_start,
    bin_end,
    mid_Mb,
    window_len,
    is_full_window,
    mean_log2,
    z_score,
    group
  )


# ebs-H3K4me3
h43_ebs_df = read.table(file.path(script_dir, "ebs-H3K4me3_clean_log2_nosm.bedgraph"), header = F)
colnames(h43_ebs_df) <- c(
  "chr_id", "start", "end", "log2_signal"
)

chr_info <- data.frame(
  chr_id = paste0("Chr", 1:5),
  length_bp = c(
    32637894,
    27782540,
    26149336,
    25768865,
    30142610
  )
)

bin_size_bp   <- 1000 * 1000
slide_step_bp <- 50 * 1000

bin_df <- lapply(seq_len(nrow(chr_info)), function(i) {
  
  chr <- chr_info$chr_id[i]
  chr_len <- chr_info$length_bp[i]
  
  starts <- seq(
    from = 0,
    to = chr_len - 1,
    by = slide_step_bp
  )
  
  data.frame(
    chr_id = chr,
    bin_start = starts,
    bin_end = pmin(starts + bin_size_bp, chr_len)
  )
  
}) %>%
  bind_rows() %>%
  mutate(
    window_len = bin_end - bin_start,
    is_full_window = window_len == bin_size_bp,
    mid_Mb = (bin_start + bin_end) / 2 / 1e6
  )

overlap_df <- h43_ebs_df %>%
  inner_join(bin_df, by = "chr_id") %>%
  filter(start < bin_end, end > bin_start) %>%
  mutate(
    overlap_bp = pmin(end, bin_end) -
      pmax(start, bin_start)
  )

bin_signal <- overlap_df %>%
  group_by(chr_id, bin_start, bin_end) %>%
  summarise(
    weighted_mean = sum(
      log2_signal * overlap_bp,
      na.rm = TRUE
    ) / sum(overlap_bp),
    .groups = "drop"
  )

h43_ebs_final <- bin_df %>%
  left_join(
    bin_signal,
    by = c("chr_id", "bin_start", "bin_end")
  ) %>%
  arrange(chr_id, bin_start) %>%
  group_by(chr_id) %>%
  mutate(
    mean_log2 = if_else(
      is_full_window,
      weighted_mean,
      NA_real_
    )
  ) %>%
  fill(mean_log2, .direction = "down") %>%
  mutate(
    chr_mean = mean(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    chr_sd = sd(
      mean_log2[is_full_window],
      na.rm = TRUE
    ),
    z_score = if_else(
      is.na(chr_sd) | chr_sd == 0,
      0,
      (mean_log2 - chr_mean) / chr_sd
    )
  ) %>%
  ungroup() %>%
  mutate(group = "ebs-H3K4me3") %>%
  select(
    chr_id,
    bin_start,
    bin_end,
    mid_Mb,
    window_len,
    is_full_window,
    mean_log2,
    z_score,
    group
  )




# Genome-wide signal profiles
# Choose the tracks and matching colors for each comparison. The last assignment is active.
combined_df <- bind_rows(spo11_final, ebs_final, h43_wt_final)
combined_df <- bind_rows(h43_ebs_final, h43_wt_final)
combined_df <- bind_rows(spo11_final, ebs_spo11_final)

combined_df <- combined_df %>%
  # Start Chr2 and Chr4 profiles at 5 Mb and 4 Mb, respectively.
  filter(!(chr_id == "Chr2" & mid_Mb < 5)) %>%
  filter(!(chr_id == "Chr4" & mid_Mb < 4))

vline_data <- tibble(
  chr_id = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
  pos1 = c(13.05,	7.3,	11.2,	5.6,	10.15),
  pos2 = c(14.377135,	9.388026,	13.595696,	6.774918,	11.947852),
  pos3 = c(17.931424,	12.146757,	15.825034,	11.168629,	15.173109),
  pos4 = c(19.05,	14.6,	18.5,	12.75,	17.2)
)

ggplot(combined_df, aes(x = mid_Mb, y = z_score, color = group)) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c("ebs-SPO11-1-FLAG" = "red", "SPO11-1-FLAG" = "darkslateblue")) + 
  scale_x_continuous(breaks = seq(0, 35, 5)) +
  scale_y_continuous(limits = c(-6, 4)) +
  
  facet_wrap(~chr_id, ncol = 5, scales = "free_x") +
  facetted_pos_scales(
    x = list(
      chr_id == "Chr2" ~ scale_x_continuous(limits = c(5, NA), breaks = seq(5, 35, 5)),
      chr_id == "Chr4" ~ scale_x_continuous(limits = c(4, NA), breaks = seq(4, 35, 5))
    )
  )+
  geom_vline(data = vline_data, aes(xintercept = pos1), 
             color = "dodgerblue", linewidth = 0.75, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos2), 
             color = "magenta", linewidth = 0.75, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos3), 
             color = "magenta", linewidth = 0.75, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos4), 
             color = "dodgerblue", linewidth = 0.75, linetype = 2) +
  coord_cartesian(clip = "on") +
  labs(
    x = "Chromosome position (Mb)",
    y = "Z-score",
    color = "Group"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

