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

# Compare peak-associated gene sets in euchromatin.
library(readxl)
library(tidyverse)

name_list <- readLines(file.path(script_dir, "peak_list.txt"))
for (name in name_list) {
  file_path <- file.path(script_dir, paste0(name, "_peak_tair12_annotated.xlsx"))
  df <- read_excel(file_path, sheet = 1)  
  assign(paste0(name, "_df"), df)
}

`EBS-FLAG_df` <- `EBS-FLAG_df`[`EBS-FLAG_df`$type == "euchromatin", ]
`SPO11-1-FLAG_df` <- `SPO11-1-FLAG_df`[`SPO11-1-FLAG_df`$type == "euchromatin", ]
`ebs-SPO11-1-FLAG_df` <- `ebs-SPO11-1-FLAG_df`[`ebs-SPO11-1-FLAG_df`$type == "euchromatin", ]
`WT-H3K4me3_df` <- `WT-H3K4me3_df`[`WT-H3K4me3_df`$type == "euchromatin", ]
`ebs-H3K4me3_df` <- `ebs-H3K4me3_df`[`ebs-H3K4me3_df`$type == "euchromatin", ]


spo11 = `SPO11-1-FLAG_df`$geneId
ebs = `EBS-FLAG_df`$geneId
wt_43 = `WT-H3K4me3_df`$geneId
ebs_43 = `ebs-H3K4me3_df`$geneId
ebs_spo11 = `ebs-SPO11-1-FLAG_df`$geneId

length(unique(spo11))
length(unique(ebs))
length(unique(wt_43))
length(unique(ebs_43))
length(unique(ebs_spo11))


# Export gene sets for the TSS profiles.
ebs_only <- setdiff(ebs, spo11)
tibble(geneId = ebs_only) %>%
  write.csv(file.path(script_dir, "gene_set", "EBS-FLAG_setdiff_SPO11-1-FLAG.csv"), row.names = FALSE)

spo11_only <- setdiff(spo11, ebs)
tibble(geneId = spo11_only) %>%
  write.csv(file.path(script_dir, "gene_set", "SPO11-1-FLAG_setdiff_EBS-FLAG.csv"), row.names = FALSE)

ebs_spo11_overlap <- intersect(ebs, spo11)
tibble(geneId = ebs_spo11_overlap) %>%
  write.csv(file.path(script_dir, "gene_set", "EBS-FLAG_intersect_SPO11-1-FLAG.csv"), row.names = FALSE)

ebs_wt_43_overlap <- intersect(ebs, wt_43)
tibble(geneId = ebs_wt_43_overlap) %>%
  write.csv(file.path(script_dir, "gene_set", "EBS-FLAG_intersect_WT-H3K4me3.csv"), row.names = FALSE)

# EBS-FLAG, SPO11-1-FLAG and WT H3K4me3 overlap
inter1 =  intersect(ebs,spo11)
inter2 = intersect(inter1,wt_43)
length(inter2)
write.csv(data.frame(gene = inter2), "3_overlap.csv", row.names = FALSE)

# EBS-FLAG, SPO11-1-FLAG and WT H3K4me3
p = venn.diagram(
  x = list(spo11,ebs,wt_43),
  category.names = c(" ", " "," "),
  filename = NULL,
  fill = c(  "darkslateblue", "#FF66CC", "lightblue"),  
  alpha = 0.5,
  lwd = 2,
  lty = "solid",
  col = c("white", "white", "white"),
  label.col = "black",
  fontfamily = "Arial",
  cex = 2,
  cat.cex = 1.2,
  main = " "
)
grid.draw(p)

total_genes <- 24493  

# RF is the observed overlap divided by the expected overlap under independence.
intersection = intersect(spo11, wt_43)
RF <- length(intersection) / ((length(spo11) * length(wt_43)) / total_genes)
# Test for enrichment using the upper tail of the hypergeometric distribution.
p_value <- phyper(length(intersect(spo11, wt_43))-1, 
                  length(spo11), 
                  total_genes-length(spo11), 
                  length(wt_43), 
                  lower.tail = FALSE)
p_value_formatted <- ifelse(p_value < 0.001, 
                            format(p_value, scientific = TRUE, digits = 3),
                            round(p_value, 3))
p_value_formatted
RF

intersection = intersect(wt_43, ebs)
RF <- length(intersection) / ((length(wt_43) * length(ebs)) / total_genes)
p_value <- phyper(length(intersect(wt_43, ebs))-1, 
                  length(wt_43), 
                  total_genes-length(wt_43), 
                  length(ebs), 
                  lower.tail = FALSE)
p_value_formatted <- ifelse(p_value < 0.001, 
                            format(p_value, scientific = TRUE, digits = 3),
                            round(p_value, 3))
p_value_formatted
RF

intersection = intersect(spo11, ebs)
RF <- length(intersection) / ((length(spo11) * length(ebs)) / total_genes)
p_value <- phyper(length(intersect(spo11, ebs))-1, 
                  length(spo11), 
                  total_genes-length(spo11), 
                  length(ebs), 
                  lower.tail = FALSE)
p_value_formatted <- ifelse(p_value < 0.001, 
                            format(p_value, scientific = TRUE, digits = 3),
                            round(p_value, 3))
p_value_formatted
RF


# WT and ebs H3K4me3
p = venn.diagram(
  x = list(ebs_43,wt_43),
  category.names = c(" ", " "),
  filename = NULL,
  fill = c(  "orchid", "lightblue"), 
  alpha = 0.5,
  lwd = 2,
  lty = "solid",
  col = c("white", "white"),
  label.col = "black",
  fontfamily = "Arial",
  cex = 2,
  cat.cex = 1.2,
  main = " "
)
grid.draw(p)

total_genes <- 24493  

intersection = intersect(wt_43, ebs_43)
RF <- length(intersection) / ((length(wt_43) * length(ebs_43)) / total_genes)
p_value <- phyper(length(intersect(wt_43, ebs_43))-1, 
                  length(wt_43), 
                  total_genes-length(wt_43), 
                  length(ebs_43), 
                  lower.tail = FALSE)
p_value_formatted <- ifelse(p_value < 0.001, 
                            format(p_value, scientific = TRUE, digits = 3),
                            round(p_value, 3))
p_value_formatted
RF

