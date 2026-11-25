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

# Sliding-window crossover rates
co_df <- read_excel(file.path(script_dir, "CO-results.xlsx"), sheet = 1)

# Use distinct sample IDs as the denominator for each group.
spore_stat <- co_df %>%
  group_by(group) %>%
  summarise(
    unique_spore = n_distinct(sample_id),
    male_CO_records = n(),
    .groups = "drop"
  )
print(spore_stat)

wt_info <- spore_stat %>% filter(grepl("WT_male", group))
ebs_info <- spore_stat %>% filter(grepl("ebs_male", group))   

wt_spore <- wt_info$unique_spore
ebs_spore <- ebs_info$unique_spore

cat("WT_male distinct sample count: ", wt_spore, "\n")
cat("ebs_male distinct sample count: ", ebs_spore, "\n")

wt_co_df  <- co_df %>% filter(grepl("WT_male", group))   
ebs_co_df  <- co_df %>% filter(grepl("ebs_male", group))   

chr_info <- data.frame(
  chr_id = paste0("Chr", 1:5),
  length_bp = c(32637894, 27782540, 26149336, 25768865, 30142610)
)

window_size_kb <- 1000
slide_step_kb  <- 50
bin_size_bp    <- window_size_kb * 1000
slide_step_bp  <- slide_step_kb * 1000

# Use 0-based windows clipped at chromosome ends.
bin_df <- lapply(1:nrow(chr_info), function(i) {
  chr <- chr_info$chr_id[i]
  chr_male_len <- chr_info$length_bp[i]
  all_starts <- seq(from = 0, to = chr_male_len, by = slide_step_bp)
  all_ends <- all_starts + bin_size_bp
  all_ends[all_ends > chr_male_len] <- chr_male_len
  data.frame(chr_id = chr, bin_start = all_starts, bin_end = all_ends)
}) %>% bind_rows()

bin_df <- bin_df %>%
  mutate(
    mid_bp = (bin_start + bin_end) / 2,
    mid_Mb = mid_bp / 1e6,
    bin_width_Mb = (bin_end - bin_start) / 1e6
  )

count_co <- function(co_data, bin_data){
  co_data %>%
    left_join(bin_data, by = "chr_id") %>%
    filter(start < bin_end & stop > bin_start) %>%
    group_by(chr_id, bin_start, bin_end, mid_bp, mid_Mb, bin_width_Mb) %>%
    summarise(co_count = sum(co_num, na.rm = TRUE), .groups = "drop")
}

wt_binned <- count_co(wt_co_df, bin_df)
ebs_binned <- count_co(ebs_co_df, bin_df)

wt_final <- bin_df %>%
  left_join(wt_binned, by = c("chr_id","bin_start","bin_end","mid_bp","mid_Mb","bin_width_Mb")) %>%
  replace_na(list(co_count = 0)) %>%
  mutate(
    cM = co_count / wt_spore * 100,
    rate = cM,
    type = "WT_male",
    group = "WT_male"
  )

ebs_final <- bin_df %>%
  left_join(ebs_binned, by = c("chr_id","bin_start","bin_end","mid_bp","mid_Mb","bin_width_Mb")) %>%
  replace_na(list(co_count = 0)) %>%
  mutate(
    cM = co_count / ebs_spore * 100,
    rate = cM,
    type = "ebs_male",
    group = "ebs_male"
  )



# Crossover fractions in non-overlapping bins
wt_plot_df <- wt_final
names(wt_plot_df)[names(wt_plot_df) == "cM"] <- "wt_co_sm"
wt_plot_df <- wt_plot_df[, !colnames(wt_plot_df) %in% c("male_co", "cM_per_Mb", "group")]

ebs_plot_df <- ebs_final
names(ebs_plot_df)[names(ebs_plot_df) == "cM"] <- "ebs_co_sm"
ebs_plot_df <- ebs_plot_df[, "ebs_co_sm", drop = FALSE]

stack_df <- cbind(wt_plot_df, ebs_plot_df)

#=========
window_size_Mb <- 1  # Mb

stack_wide <- stack_df %>%
  group_by(chr_id) %>%
  mutate(
    window_id = floor((bin_start / 1e6) / window_size_Mb),
    window_mid_Mb = (window_id * window_size_Mb) + (window_size_Mb / 2)
  ) %>%
  group_by(chr_id, window_id, window_mid_Mb) %>%
  summarise(
    wt_sum = sum(wt_co_sm, na.rm = TRUE),
    ebs_sum = sum(ebs_co_sm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    male = wt_sum + ebs_sum,
    wt_prop = ifelse(male > 0, wt_sum / male * 100, 0),
    ebs_prop = ifelse(male > 0, ebs_sum / male * 100, 0)
  ) %>%
  select(chr_id, window_mid_Mb, wt_prop, ebs_prop)

print(head(stack_wide))

plot_df_long <- stack_wide %>%
  pivot_longer(
    cols = c(wt_prop, ebs_prop),
    names_to = "group",
    values_to = "percent"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("wt_prop", "ebs_prop"),
      labels = c("WT_male", "EBS_male")
    )
  )

ebs_line_data <- stack_wide %>%
  select(chr_id, window_mid_Mb, ebs_prop) %>%
  mutate(
    ebs_prop_percent = ebs_prop / 100
  )

# Centromere positions in Mb
vline_data <- tibble(
  chr_id = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),
  pos1 = c(13.05,	7.3,	11.2,	5.6,	10.15),
  pos2 = c(14.377135,	9.388026,	13.595696,	6.774918,	11.947852),
  pos3 = c(17.931424,	12.146757,	15.825034,	11.168629,	15.173109),
  pos4 = c(19.05,	14.6,	18.5,	12.75,	17.2)
)

ggplot() +
  geom_col(
    data = plot_df_long,
    aes(x = factor(window_mid_Mb), y = percent/100, fill = group),
    position = "fill",
    width = 0.9,
    color = "white", size = 0.1
  ) +
  geom_line(
    data = ebs_line_data,
    aes(x = factor(window_mid_Mb), y = ebs_prop/100, group = 1),
    color = "gray40",
    linewidth = 0.8,
    group = 1
  ) +
  geom_point(
    data = ebs_line_data,
    aes(x = factor(window_mid_Mb), y = ebs_prop/100),
    color = "gray40",
    size = 1.5
  ) +
  geom_hline(
    yintercept = 0.5,  # 50% = 0.5
    color = "gray50",
    linetype = "dashed",
    linewidth = 0.5,
    alpha = 0.6
  ) +
  scale_fill_manual(
    values = c("WT_male" = "lightblue", "EBS_male" = "salmon")
  ) +
  scale_y_continuous(
    labels = scales::percent_format(),
    position = "left",
    expand = c(0, 0),
    limits = c(0, 1)
  ) +
  scale_x_discrete(
    breaks = function(x) {
      x[seq(1, length(x), by = 3)]
    }
  ) +
  facet_wrap(~chr_id, ncol = 5, scales = "free_x") +
  labs(
    x = "Chromosome position (Mb)",
    y = "Proportion (%)",
    fill = "Group"
  ) +
  geom_vline(data = vline_data, aes(xintercept = pos1), 
             color = "dodgerblue", linewidth = 1, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos2), 
             color = "magenta", linewidth = 1, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos3), 
             color = "magenta", linewidth = 1, linetype = 2) +
  geom_vline(data = vline_data, aes(xintercept = pos4), 
             color = "dodgerblue", linewidth = 1, linetype = 2) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 15),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y.left = element_line(),
    axis.text.y.left = element_text()
  )


