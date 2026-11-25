library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(zoo)
library(ggh4x)
library(corrplot)
library(ggsignif)

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

co_df <- read_excel(file.path(script_dir, "CO-results.xlsx"), sheet = 1)

# Use distinct sample IDs as the denominator for each group.
spore_stat <- co_df %>%
  group_by(group) %>%
  summarise(
    unique_spore = n_distinct(sample_id),
    total_CO_records = n(),
    .groups = "drop"
  )
print(spore_stat)

wt_info <- spore_stat %>% filter(group == "WT_female")
ebs_info <- spore_stat %>% filter(group == "ebs_female")

wt_spore <- wt_info$unique_spore
ebs_spore <- ebs_info$unique_spore

cat("WT_female distinct sample count: ", wt_spore, "\n")
cat("ebs_female distinct sample count: ", ebs_spore, "\n")

# LRZ2 intervals define heterochromatin in this analysis.
hetero_regions <- data.frame(
  chr = c("Chr1","Chr2","Chr3","Chr4","Chr5"),
  LRZ2_left = c(13050000,	7300000,	11200000,	5600000,	10150000),
  LRZ2_right = c(19050000,	14600000,	18500000,	12750000,	17200000)
)

co_df_mark <- co_df %>%
  filter(group %in% c("WT_female", "ebs_female")) %>%
  rowwise() %>%
  mutate(
    region_type = if(any(
      hetero_regions$chr == chr_id &
      start <= hetero_regions$LRZ2_right &
      stop >= hetero_regions$LRZ2_left
    )){"hete"}else{"euch"}
  ) %>%
  ungroup()

co_summary <- co_df_mark %>%
  group_by(sample_id, group, region_type) %>%
  summarise(co_sum = sum(co_num), .groups = "drop") %>%
  pivot_wider(
    names_from = region_type,
    values_from = co_sum,
    values_fill = 0
  )

# Express chromatin-specific CO counts as fractions per sample.
co_eu_he_df <- co_summary %>%
  mutate(
    total = euch + hete,
    euch = euch / total,
    hete = hete / total
  ) %>%
  select(sample_id, group, euch, hete)

co_eu_he_df$group <- factor(co_eu_he_df$group, levels = c("WT_female", "ebs_female"))
co_eu_he_df$group <- as.factor(co_eu_he_df$group)

library(ggplot2)
library(ggsignif)

# Euchromatic CO fraction
ggplot(co_eu_he_df, aes(x = group, y = euch, color = group)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  stat_summary(
    fun.data = mean_se,   # mean ± SEM
    geom = "errorbar",
    width = 0.25,
    linewidth = 0.25,
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 1.5,
    color = "black"
  ) +
  geom_signif(
    comparisons = list(c("WT_female", "ebs_female")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    y_position = 1.04,
    color = "black",
    textsize = 4.5
  ) +
  scale_color_manual(values = c("WT_female"="#3498db", "ebs_female"="#e74c3c")) +
  scale_x_discrete(expand = c(0.4, 0.4)) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1.1, 0.2),
    labels = paste0(seq(0, 110, 20), "%")
  ) +
  labs(
    x = "",
    y = "",
    title = " "
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# ========== Pericentromeric ==========
ggplot(co_eu_he_df, aes(x = group, y = hete, color = group)) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.8) +
  stat_summary(
    fun.data = mean_se,   # mean ± SEM
    geom = "errorbar",
    width = 0.25,
    linewidth = 0.25,
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 1.5,
    color = "black"
  ) +
  geom_signif(
    comparisons = list(c("WT_female", "ebs_female")),
    test = "wilcox.test",
    map_signif_level = TRUE,
    y_position = 1.04,
    color = "black",
    textsize = 4.5
  ) +
  scale_color_manual(values = c("WT_female"="#3498db", "ebs_female"="#e74c3c")) +
  scale_x_discrete(expand = c(0.4, 0.4)) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1.1, 0.2),
    labels = paste0(seq(0, 110, 20), "%")
  ) +
  labs(
    x = "",
    y = "",
    title = " "
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )


# Euchromatin
stat_euch <- co_eu_he_df %>%
  group_by(group) %>%
  summarise(
    mean_val = mean(euch),
    sd = sd(euch),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )
cat("==== Euchromatin (euch) mean ± SEM ====\n")
print(stat_euch)

# Heterochromatin
stat_hete <- co_eu_he_df %>%
  group_by(group) %>%
  summarise(
    mean_val = mean(hete),
    sd = sd(hete),
    n = n(),
    sem = sd / sqrt(n),
    .groups = "drop"
  )
cat("\n==== Pericentromeric mean ± SEM ====\n")
print(stat_hete)




