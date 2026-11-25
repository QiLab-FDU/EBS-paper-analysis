library(Biostrings)
library(rtracklayer)
library(dplyr)
library(readr)
library(openxlsx)
library(stringr)
library(GenomicRanges)

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
reference_dir <- file.path(script_dir, "..", "..", "tair12_genome")

overlap_df <- read.csv(
  file.path(script_dir, "..", "peak_annotation", "gene_set", "3_overlap.csv"),
  header = TRUE
)

gene_vector <- unique(overlap_df$gene)
length(gene_vector)

tair12_gff <- read.table(
  file.path(reference_dir, "TAIR12.gtf"),
  sep = "\t",
  quote = "",
  header = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

tair12_gff$gene <- NA

for(g in gene_vector){
  idx <- grepl(g, tair12_gff$V9, fixed = TRUE)
  tair12_gff$gene[idx] <- g
}

# Keep the first matching GTF record for each gene.
target_df <- tair12_gff[!is.na(tair12_gff$gene), ]
target_df <- target_df[!duplicated(target_df$gene), ]

unmatched_genes <- setdiff(gene_vector, unique(target_df$gene))
print(unmatched_genes)


target_df$V4 <- as.integer(target_df$V4)
target_df$V5 <- as.integer(target_df$V5)

target_df$TSS <- ifelse(target_df$V7 == "+", target_df$V4, target_df$V5)

# Define the 2 kb upstream promoter on each strand.
target_df$prom_start <- ifelse(
  target_df$V7 == "+",
  target_df$TSS - 2000,
  target_df$TSS + 1
)

target_df$prom_end <- ifelse(
  target_df$V7 == "+",
  target_df$TSS - 1,
  target_df$TSS + 2000
)

target_df$genebody_start <- ifelse(
  target_df$V7 == "+",
  target_df$TSS,
  target_df$V4
)

target_df$genebody_end <- ifelse(
  target_df$V7 == "+",
  target_df$V5,
  target_df$TSS
)

target_df$prom_start <- pmax(target_df$prom_start, 1)

# Use one contiguous interval for the promoter and gene body.
target_df$prom_gb_start <- ifelse(
  target_df$V7 == "+",
  target_df$prom_start,
  target_df$genebody_start
)
target_df$prom_gb_end <- ifelse(
  target_df$V7 == "+",
  target_df$genebody_end,
  target_df$prom_end
)

target_df$prom_gb_start <- pmax(target_df$prom_gb_start,1)


table(target_df$prom_start < target_df$prom_end)
table(target_df$genebody_start < target_df$genebody_end)
table(target_df$prom_gb_start < target_df$prom_gb_end)

write.csv(target_df, file.path(script_dir, "3_overlap_tair12_pos.csv"), row.names = FALSE)

# Extract strand-oriented sequences for motif discovery.
genome <- readDNAStringSet(file.path(reference_dir, "TAIR12.fna.gz"))
names(genome) <- sub(" .*", "", names(genome))

gr_prom_gb <- GRanges(
  seqnames = target_df$V1,
  ranges   = IRanges(start = target_df$prom_gb_start,
                     end   = target_df$prom_gb_end),
  strand   = target_df$V7
)

prom_gb_seqs <- DNAStringSet()
for(i in seq_along(gr_prom_gb)){
  chr <- as.character(seqnames(gr_prom_gb)[i])
  s <- start(gr_prom_gb)[i]
  e <- end(gr_prom_gb)[i]
  st <- as.character(strand(gr_prom_gb)[i])
  
  seq_sub <- subseq(genome[[chr]], start = s, end = e)
  
  if(st == "-"){
    seq_sub <- reverseComplement(seq_sub)
  }
  prom_gb_seqs <- c(prom_gb_seqs, DNAStringSet(seq_sub))
}
names(prom_gb_seqs) <- target_df$gene
writeXStringSet(prom_gb_seqs, file.path(script_dir, "3_overlap_promoter_genebody.fasta"))


