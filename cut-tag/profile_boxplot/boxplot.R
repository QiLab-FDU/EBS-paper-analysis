library(ggplot2)
library(dplyr)

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


# WT H3K4me3 signal by EBS-FLAG and SPO11-1-FLAG gene-set overlap.
# Average signal across the 0.5 kb on either side of each TSS.
overlap_log_df = read.csv(file.path(script_dir, "wt-43_spo11-ebs-overlap_tss.csv"),header = TRUE, stringsAsFactors = FALSE)
overlap_log_df$log2Fold <- rowMeans(overlap_log_df[, names(overlap_log_df) %in% paste0("X", 6:15)] , na.rm = TRUE)
overlap_log_df <- overlap_log_df[, c("gene", "log2Fold")]

spo11_flag_log_df = read.csv(file.path(script_dir, "wt-43_spo11-ebs-nooverlap_tss.csv"),header = TRUE, stringsAsFactors = FALSE)
spo11_flag_log_df$log2Fold <- rowMeans(spo11_flag_log_df[, names(spo11_flag_log_df) %in% paste0("X", 6:15)] , na.rm = TRUE)
spo11_flag_log_df <- spo11_flag_log_df[, c("gene", "log2Fold")]

ebs_flag_log_df = read.csv(file.path(script_dir, "wt-43_ebs-spo11-nooverlap_tss.csv"),header = TRUE, stringsAsFactors = FALSE)
ebs_flag_log_df$log2Fold <- rowMeans(ebs_flag_log_df[, names(ebs_flag_log_df) %in% paste0("X", 6:15)] , na.rm = TRUE)
ebs_flag_log_df <- ebs_flag_log_df[, c("gene", "log2Fold")]

overlap_log_df$group <- factor("overlap", levels = c("SPO11-1-FLAG", "overlap", "EBS-FLAG" ))
spo11_flag_log_df$group <- factor("SPO11-1-FLAG", levels = c("SPO11-1-FLAG", "overlap", "EBS-FLAG"))
ebs_flag_log_df$group <- factor("EBS-FLAG", levels = c("SPO11-1-FLAG", "overlap", "EBS-FLAG"))

combined_df <- bind_rows(overlap_log_df, spo11_flag_log_df, ebs_flag_log_df)

ggplot(combined_df, aes(x = group, y = log2Fold, fill = group)) +
  geom_boxplot(
    width = 0.4, 
    alpha = 0.7,
    coef = 1.5,
    lwd = 0.6,
    outlier.size = 1
  ) +
  scale_y_continuous(limits = c(-0.2, 1)) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("SPO11-1-FLAG"="darkslateblue", "overlap"="darkorchid", "EBS-FLAG"= "#FF66CC" )) +
  labs(title = "", 
       y = "Log2FoldChange (wt_43/wt_43_igg)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, lty = 3, col = "black", lwd = 0.5)


# One-sided Wilcoxon rank-sum tests, with x greater than y.
x <- overlap_log_df$log2Fold   
y <- spo11_flag_log_df$log2Fold  

wilcox.test(x, y, alternative = "greater")


x <- overlap_log_df$log2Fold   
y <- ebs_flag_log_df$log2Fold  

wilcox.test(x, y, alternative = "greater")

x <- spo11_flag_log_df$log2Fold   
y <- ebs_flag_log_df$log2Fold  

wilcox.test(x, y, alternative = "greater")


# SPO11-1-FLAG signal at genes shared by EBS-FLAG and WT H3K4me3.
# Average signal across the 0.5 kb on either side of each TSS.
wt_log_df = read.csv(file.path(script_dir, "spo11_alone-ebs_wt-43-overlap_tss.csv"),header = TRUE, stringsAsFactors = FALSE)
wt_log_df$log2Fold <- rowMeans(wt_log_df[, names(wt_log_df) %in% paste0("X", 6:15)] , na.rm = TRUE)
wt_log_df <- wt_log_df[, c("gene", "log2Fold")]

mut_log_df = read.csv(file.path(script_dir, "spo11-ebs_alone-ebs_wt-43-overlap_tss.csv"),header = TRUE, stringsAsFactors = FALSE)
mut_log_df$log2Fold <- rowMeans(mut_log_df[, names(mut_log_df) %in% paste0("X", 6:15)] , na.rm = TRUE)
mut_log_df <- mut_log_df[, c("gene", "log2Fold")]

wt_log_df$group <- factor("SPO11-1-FLAG", levels = c("SPO11-1-FLAG", "ebs-SPO11-1-FLAG"))
mut_log_df$group <- factor("ebs-SPO11-1-FLAG", levels = c("SPO11-1-FLAG", "ebs-SPO11-1-FLAG"))

combined_df <- bind_rows(wt_log_df, mut_log_df)

ggplot(combined_df, aes(x = group, y = log2Fold, fill = group)) +
  geom_boxplot(
    width = 0.4, 
    alpha = 0.7,
    coef = 1.5,
    lwd = 0.6,
    outlier.size = 1
  ) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  stat_summary(fun = "median", geom = "point", shape = 23, size = 3, fill = "white") +
  scale_fill_manual(values = c("SPO11-1-FLAG" = "darkslateblue", "ebs-SPO11-1-FLAG" = "red")) +
  labs(title = "", 
       y = "Log2FoldChange (sample/igg)") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 0, lty = 3, col = "black", lwd = 0.5)


# One-sided Wilcoxon rank-sum tests, with x greater than y.
x <- wt_log_df$log2Fold   
y <- mut_log_df$log2Fold  

wilcox.test(x, y, alternative = "greater")

