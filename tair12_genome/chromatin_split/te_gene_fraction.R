library(tidyverse)
library(data.table)


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

# Sliding windows for gene and transposable-element coverage.
chr_info <- data.frame(
  chr_id = paste0("Chr", 1:5),
  length_bp = c(32637894, 27782540, 26149336, 25768865, 30142610)
)

window_size_kb <- 1000
slide_step_kb  <- 50
bin_size_bp    <- window_size_kb * 1000
slide_step_bp  <- slide_step_kb * 1000

bin_df <- lapply(1:nrow(chr_info), function(i) {
  chr <- chr_info$chr_id[i]
  chr_total_len <- chr_info$length_bp[i]
  all_starts <- seq(from = 0, to = chr_total_len, by = slide_step_bp)
  all_ends <- all_starts + bin_size_bp
  all_ends[all_ends > chr_total_len] <- chr_total_len
  data.frame(chr_id = chr, bin_start = all_starts, bin_end = all_ends)
}) %>% bind_rows()

bin_df <- bin_df %>%
  mutate(
    mid_bp = (bin_start + bin_end) / 2,
    mid_Mb = mid_bp / 1e6,
    win_len = bin_end - bin_start
  )

setDT(bin_dt <- bin_df)

gff_raw <- read_delim(
  file.path(script_dir, "..", "TAIR12_1Feb26.gff3"),
  delim = "\t",
  comment = "#",
  col_names = c("chr_id","source","type","start","end","score","strand","phase","attr"),
  show_col_types = FALSE
)

feature_dt <- gff_raw %>%
  filter(type %in% c("gene", "mobile_element")) %>%
  select(chr_id, type, feat_start = start, feat_end = end) %>%
  setDT()

overlap_stat <- function(dt_window, dt_feature, select_type) {
  sub_dt <- dt_feature[type == select_type]
  
  setkey(sub_dt, chr_id, feat_start, feat_end)
  setkey(dt_window, chr_id, bin_start, bin_end)
  
  overlap_result <- foverlaps(
    dt_window,
    sub_dt,
    by.x = c("chr_id", "bin_start", "bin_end"),
    by.y = c("chr_id", "feat_start", "feat_end"),
    nomatch = 0
  )
  
  overlap_result[, overlap_len := pmin(bin_end, feat_end) - pmax(bin_start, feat_start)]
  overlap_result <- overlap_result[overlap_len > 0]
  
  sum_dt <- overlap_result[, .(total_overlap = sum(overlap_len)),
                           by = .(chr_id, bin_start, bin_end, mid_bp, mid_Mb, win_len)]
  
  out <- dt_window[sum_dt, on = .(chr_id, bin_start, bin_end, mid_bp, mid_Mb, win_len)]
  out[is.na(total_overlap), total_overlap := 0]
  
  return(out)
}

gene_dt <- overlap_stat(bin_dt, feature_dt, "gene")
mobile_dt   <- overlap_stat(bin_dt, feature_dt, "mobile_element")

stat_final <- merge(
  gene_dt[, .(chr_id, bin_start, bin_end, mid_bp, mid_Mb, win_len, gene_bp = total_overlap)],
  mobile_dt[, .(chr_id, bin_start, bin_end, mid_bp, mid_Mb, win_len, mobile_bp = total_overlap)],
  by = c("chr_id", "bin_start", "bin_end", "mid_bp", "mid_Mb", "win_len")
) %>%
  mutate(
    gene_frac = gene_bp / win_len,
    mobile_frac   = mobile_bp / win_len
  )

# fwrite(stat_final, file.path(script_dir, "window_gene_mobile_fraction.txt"), sep = "\t")


# Find transitions in the dominant feature along each chromosome.
library(data.table)
setDT(stat_final)

stat_final[, dom := fcase(
  gene_frac > mobile_frac, "gene",
  gene_frac < mobile_frac, "mobile",
  default = "equal"
)]

stat_final <- stat_final[order(chr_id, mid_Mb)]
stat_final[, prev_dom := shift(dom), by = chr_id]

switch_points <- stat_final[dom != prev_dom & !is.na(prev_dom)]

res_boundary <- switch_points[, .(
  chr_id,
  boundary_Mb = mid_Mb,
  left_side = prev_dom,
  right_side = dom,
  gene_frac,
  mobile_frac
)]

print(res_boundary)


display_step_Mb <- 5

plot_df <- as_tibble(stat_final) %>%
  pivot_longer(
    cols = c(gene_frac, mobile_frac),
    names_to = "feature",
    values_to = "fraction"
  )

vline_data <- tibble(
  chr_id = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
  pos1 = c(13.05,	7.3,	11.2,	5.6,	10.15),
  pos2 = c(14.189619,	8.404717,	13.538464,	8.10951,	12.349536),
  pos3 = c(17.781818,	12.038746,	16.979273,	11.354181,	15.487525),
  pos4 = c(19.05,	14.6,	18.5,	12.75,	17.2)
)

ggplot(plot_df, aes(x = mid_Mb, y = fraction, color = feature)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = 3, color = "orange", linewidth = 0.7) +
  scale_color_manual(
    values = c("gene_frac" = "salmon", "mobile_frac" = "dodgerblue"),
    labels = c("Gene fraction", "TE fraction")
  ) +
  scale_x_continuous(breaks = seq(0, 35, display_step_Mb)) +
  scale_y_continuous(limits = c(0, 1), position = "left") +
  facet_wrap(~chr_id, ncol = 5, scales = "free_x") +
  labs(
    x = "Chromosome position (Mb)",
    y = "Sequence fraction",
    color = "Feature"
  ) +
  geom_vline(data = vline_data, aes(xintercept = pos1), 
             color = "limegreen", linewidth = 0.5, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos2), 
             color = "red", linewidth = 0.5, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos3), 
             color = "red", linewidth = 0.5, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos4), 
             color = "limegreen", linewidth = 0.5, linetype = 2) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )


