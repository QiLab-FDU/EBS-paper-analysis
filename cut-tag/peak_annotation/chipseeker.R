library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)
library(grid)
library(futile.logger)
library(VennDiagram)
library(dplyr)
library(fuzzyjoin)
library(writexl)

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

# Annotate peaks against TAIR12.
merge = readLines(file.path(script_dir, "peak_list.txt"))
txdb_tair12 <- makeTxDbFromGFF(
  file = file.path(script_dir, "..", "..", "tair12_genome", "TAIR12.gtf"),
  format = "gtf",
  organism = "Arabidopsis thaliana"
)

chr_df = read.csv(file.path(script_dir, "region_type.csv"), header = T)

# narrowPeak V9 stores -log10(q-value).

for (name in merge) {
  path <- file.path(script_dir, "peak_files", paste0(name, "_peaks.narrowPeak"))
  peaks <- readPeakFile(path)
  
  cat(name, " | readPeak peak number:", length(peaks), "\n")
  
  peakAnno <- annotatePeak(peaks, tssRegion = c(-1000, 1000), TxDb = txdb_tair12)
  anno_df <- as.data.frame(peakAnno)
  cat(name, " | after annotatePeak:", nrow(anno_df), "\n")
  
  sorted_df <- anno_df[order(-anno_df$V9), ]
  
  name_vector <- rep(name, nrow(sorted_df))
  new_sorted_df <- cbind(name_vector, sorted_df)
  new_sorted_df$rank <- 1:nrow(new_sorted_df)
  
  cat(name, " | before overlap assignment:", nrow(new_sorted_df), "\n")
  
  peak_gr <- GRanges(
    seqnames = new_sorted_df$seqnames,
    ranges = IRanges(start = new_sorted_df$start, end = new_sorted_df$end)
  )
  chr_gr <- GRanges(
    seqnames = chr_df$chr_ref,
    ranges = IRanges(start = chr_df$start_ref, end = chr_df$end_ref),
    type = chr_df$type
  )
  
  hits <- findOverlaps(peak_gr, chr_gr)
  hit_peak_idx <- queryHits(hits)
  hit_chr_idx <- subjectHits(hits)
  
  overlap_width <- width(pintersect(peak_gr[hit_peak_idx], chr_gr[hit_chr_idx]))
  peak_width <- width(peak_gr[hit_peak_idx])
  overlap_ratio <- overlap_width / peak_width
  
  hit_tb <- tibble(
    peak_row = hit_peak_idx,
    chr_row = hit_chr_idx,
    ratio = overlap_ratio,
    type = mcols(chr_gr)$type[hit_chr_idx]
  )
  
  # Keep the region covering the largest fraction of each peak.
  hit_best <- hit_tb %>%
    group_by(peak_row) %>%
    slice_max(ratio, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  new_annotated_df <- new_sorted_df %>%
    mutate(peak_row = row_number()) %>%
    left_join(hit_best %>% select(peak_row, ratio, type), by = "peak_row") %>%
    # Assign NRZ or LRZ when it covers at least half of the peak.
    mutate(
      type = case_when(
        !is.na(ratio) & ratio >= 0.5 ~ type,
        TRUE ~ "euchromatin"
      )
    ) %>%
    select(-peak_row) %>%
    arrange(type)
  
  cat(name, " | after overlap ratio assignment:", nrow(new_annotated_df), "\n")
  
  df_list <- list(
    "all_annotated" = new_annotated_df
  )
  
  write_xlsx(
    x = df_list,
    path = file.path(script_dir, paste0(name, "_peak_tair12_annotated.xlsx"))
  )
}


# Plot annotation proportions for euchromatic peaks.
name_list <- readLines(file.path(script_dir, "peak_list.txt"))

anno_df <- data.frame()

for (name in name_list) {
  file_path <- file.path(script_dir, paste0(name, "_peak_tair12_annotated.xlsx"))
  df <- readxl::read_excel(file_path, sheet = 1)
  
  df_eu <- df %>%
    filter(type == "euchromatin") %>%
    mutate(sample = name)
  
  anno_df <- bind_rows(anno_df, df_eu)
}

anno_df <- anno_df %>%
  mutate(anno_merge = case_when(
    grepl("Promoter", annotation) ~ "Promoter",
    grepl("5' UTR|3' UTR|Exon|Intron", annotation) ~ "Gene body",
    grepl("Downstream|Distal Intergenic", annotation) ~ "Intergenic",
    TRUE ~ "Other"
  ))

color_map <- c(
  "Promoter" = "#f8c8d0",
  "Gene body" = "#a8d8ea",
  "Intergenic" = "#ffffa0",
  "Other" = "#d1d1e0"
)

draw_donut <- function(sample_name){
  sub_df <- anno_df %>% filter(sample == !!sample_name)
  
  stat_tb <- sub_df %>%
    count(anno_merge, name = "cnt") %>%
    mutate(
      raw_pct = cnt / sum(cnt) * 100
    ) %>%
    arrange(desc(cnt))
  
  # Add the rounding remainder to the category with the largest fractional part.
  raw <- stat_tb$raw_pct
  int_part <- floor(raw * 100) / 100
  remain <- 100 - sum(int_part)
  decimal <- raw - int_part
  add_idx <- order(decimal, decreasing = TRUE)[1]
  int_part[add_idx] <- int_part[add_idx] + remain
  stat_tb$percent <- round(int_part, 2)
  
  stat_tb <- stat_tb %>%
    mutate(
      y_end = cumsum(cnt),
      y_start = lag(y_end, default = 0),
      y_mid = (y_start + y_end) / 2
    )
  
  cat("\n==========", sample_name, " annotation summary ==========\n")
  print(stat_tb %>% select(anno_merge, cnt, raw_pct, percent))
  cat("Rounded percentage total =", sum(stat_tb$percent), "\n")
  
  p <- ggplot(stat_tb, aes(x = 2, y = cnt, fill = anno_merge)) +
    geom_col(width = 1) +
    coord_polar("y", start = 0) +
    annotate("rect", xmin = 1, xmax = 1.6, ymin = -Inf, ymax = Inf, fill = "white") +
    xlim(1, 3.3) +
    geom_text(
      aes(y = y_mid, label = paste0(percent, "%")),
      x = 2.35, size = 0, inherit.aes = TRUE
    ) +
    scale_fill_manual(values = color_map) +
    labs(
      title = paste0(sample_name, " Peak Annotation Proportion"),
      fill = "Annotation"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      legend.position = "bottom"
    )
  return(p)
}

sample_order <- name_list

plot_list <- list()
for(s in sample_order){
  p <- draw_donut(s)
  plot_list[[s]] <- p
  print(p)

}


