library(rtracklayer)
library(dplyr)
library(fuzzyjoin)
library(tidyr)

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

gtf <- import(file.path(script_dir, "..", "TAIR12.gtf"), format = "gtf")
gtf_df <- as.data.frame(gtf)
gtf_df$seqnames <- gsub("Chr", "", gtf_df$seqnames)
gtf_df <- gtf_df[gtf_df$type == "transcript", ]

chr_df = read.csv(file.path(script_dir, "chr_region.csv"), header = T)

# Require the full transcript to fall inside a region.
gtf_df <- gtf_df %>%
  ungroup() %>%
  fuzzy_left_join(
    chr_df %>% ungroup(),
    by = c(
      "seqnames" = "chr",   
      "start" = "start",
      "end" = "end"
    ),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
  select(everything(), region) %>%
  mutate(region = replace_na(region, "euchromatin")) 

euch_df = gtf_df %>% filter(region == "euchromatin") 

write.csv(euch_df, file.path(script_dir, "euch-genes.csv"))
