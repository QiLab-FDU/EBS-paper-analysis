# EBS bioinformatics analysis

Analysis scripts and input tables for the manuscript *Histone reader EBS acts as a chromatin-to-axis tether directing SPO11-1 localization and crossover patterning in plants*.

## Data and reference

- The CUT&Tag SRA BioProject accession is `PRJNA1567450`. Raw sequencing reads have not been publicly released and are excluded from this repository.
- TAIR12 is the reannotation of the *Arabidopsis thaliana* Col-CC.v2 genome. The reference assembly and the TAIR12 annotation have separate records. [NCBI Col-CC.v2 assembly (GCA_028009825.2)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028009825.2/) · [NCBI TAIR12 assembly (GCA_978657495.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_978657495.1/) · [NCBI TAIR12 BioProject (PRJEB100887)](https://www.ncbi.nlm.nih.gov/bioproject/1358071)
- Use the bundled FASTA, GTF, GFF3 and gene intervals to reproduce the reference inputs used by the scripts. The FASTA uses the chromosome names `Chr1` through `Chr5`. [Bundled reference files](tair12_genome/)
- The crossover workbook contains recorded CO events. Gametes with zero events are absent from the workbook. The boxplot script uses the experimental group totals to include them. See [CO-results.xlsx](Crossover/CO-results.xlsx).

## Repository layout

Paths shown below and in the script tables are relative to the repository root, the directory containing this README.

```text
Crossover/                            Crossover scripts and CO-results.xlsx
cut-tag/
  rawdata_processing.sh              Original read-processing commands
  peak_annotation/                   Peak calls, annotation scripts and gene sets
    peak_files/                       MACS3 narrowPeak inputs
    gene_set/                         Gene lists for the motif and TSS analyses
  motif/                             Sequence extraction script and FASTA
  metaplot/                          Genome-wide signal scripts and bedGraph inputs
  peak_fraction/                     Peak overlap calculation and plotting scripts
    data/                             Included peak overlap CSV results
  profile_boxplot/                   TSS signal tables and boxplot workflow
    matrix/                           Included deepTools matrices
    plot_bed/                         BED inputs and the BED extraction script
tair12_genome/                       TAIR12.fna.gz, TAIR12.gtf and TAIR12_1Feb26.gff3
  chromatin_split/                   Chromatin region files and scripts
```

## Software

| Software | Version used |
| --- | --- |
| R | 4.4.2 |
| Trim Galore | 0.6.10 |
| Bowtie2 | 2.5.4 |
| SAMtools | 1.21 |
| Sambamba | 1.0.1 |
| BEDTools | 2.31.1 |
| MACS3 | 3.0.2 |
| deepTools, including `bigwigCompare` | 3.5.6 |
| UCSC `bedGraphToBigWig` | 2.10 |

Coverage-track conversion also uses UCSC `bigWigToBedGraph`. The BED extraction script uses Python 3 and its standard library.

The R scripts load these packages as needed.

- **CRAN** — `corrplot`, `data.table`, `dplyr`, `futile.logger`, `fuzzyjoin`, `ggh4x`, `ggplot2`, `ggsignif`, `openxlsx`, `readr`, `readxl`, `scales`, `stringr`, `tidyr`, `tidyverse`, `VennDiagram`, `writexl`, `zoo`
- **Bioconductor** — `Biostrings`, `ChIPseeker`, `clusterProfiler`, `GenomicRanges`, `IRanges`, `org.At.tair.db`, `rtracklayer`, `TxDb.Athaliana.BioMart.plantsmart28`

## Crossover analysis

| Figure | Script | Analysis |
| --- | --- | --- |
| 3B | [`Crossover/CO_boxplot.R`](Crossover/CO_boxplot.R) | CO counts per gamete |
| 3D | [`Crossover/CO-fraction_male.R`](Crossover/CO-fraction_male.R) | Male crossover distribution |
| 3F | [`Crossover/CO-scatter_male.R`](Crossover/CO-scatter_male.R) | Male CO fractions by chromatin type |
| S7A | [`Crossover/CO-fraction_female.R`](Crossover/CO-fraction_female.R) | Female crossover distribution |
| S7B | [`Crossover/CO-scatter_female.R`](Crossover/CO-scatter_female.R) | Female CO fractions by chromatin type |

## CUT&Tag analysis

| Step | Scripts and inputs | Figures |
| --- | --- | --- |
| Read preprocessing | Trim and align paired-end reads, filter alignments, remove duplicates, and exclude rDNA regions using [`cut-tag/rawdata_processing.sh`](cut-tag/rawdata_processing.sh). | |
| Peak calling and annotation | Call peaks with [`cut-tag/peak_annotation/masc3.sh`](cut-tag/peak_annotation/masc3.sh), annotate them against TAIR12 with [`cut-tag/peak_annotation/chipseeker.R`](cut-tag/peak_annotation/chipseeker.R), and compare euchromatic peak-associated gene sets with [`cut-tag/peak_annotation/venn.R`](cut-tag/peak_annotation/venn.R). | 4A, 4B, S11A, S11B |
| Motif discovery | [`cut-tag/motif/target_sequence.R`](cut-tag/motif/target_sequence.R) extracts strand-oriented sequences spanning the 2 kb upstream promoter and gene body for the 1,273 genes in the supplied [three-way overlap](cut-tag/peak_annotation/gene_set/3_overlap.csv). The resulting [FASTA](cut-tag/motif/3_overlap_promoter_genebody.fasta) was submitted to MEME with default settings. | 4C |
| Genome-wide profiles | Generate bedGraph tracks with [`cut-tag/metaplot/bedgraph_generate.sh`](cut-tag/metaplot/bedgraph_generate.sh) and plot windowed signals with [`cut-tag/metaplot/metaplot.R`](cut-tag/metaplot/metaplot.R). The final active comparison is SPO11-1-FLAG versus ebs-SPO11-1-FLAG. For the other comparisons, select the relevant tracks and matching colors in the final plotting section. | 4D, 5E, S11C |
| Peak overlap fractions | [`cut-tag/peak_fraction/calculate_fraction.R`](cut-tag/peak_fraction/calculate_fraction.R) computes direct interval overlaps and a chromosome-wise circular-shift control. [`cut-tag/peak_fraction/plot_fraction.R`](cut-tag/peak_fraction/plot_fraction.R) reads the included CSV results. | 4E, S12C |
| TSS signal profiles and boxplots | Generate coverage tracks and deepTools matrices, then use the included TSS signal tables for the boxplots. See the [TSS profile workflow](cut-tag/profile_boxplot/README.md). | S9, S12A, S12B |

The reference utilities identify euchromatic transcripts and calculate gene and transposable-element coverage in sliding windows. [`tair12_genome/chromatin_split/euchromatin_genes.R`](tair12_genome/chromatin_split/euchromatin_genes.R) · [`tair12_genome/chromatin_split/te_gene_fraction.R`](tair12_genome/chromatin_split/te_gene_fraction.R)

## Running the scripts

Keep the directory structure intact. R and Python resolve the bundled input files relative to their script directories. The included workbook, annotated peaks, bedGraph tracks and CSV tables are sufficient for the plotting examples below. From the repository root, run these commands in an interactive R session.

```r
source("Crossover/CO_boxplot.R")
source("cut-tag/peak_fraction/plot_fraction.R")
source("cut-tag/profile_boxplot/boxplot.R", print.eval = TRUE)
```

The main file handoffs are listed below. The paths refer to files already included in the repository unless an output directory is marked as generated.

| Step | Input path | Output or next step |
| --- | --- | --- |
| Crossover plots | `Crossover/CO-results.xlsx` | Plots in the interactive R session |
| Peak annotation | `cut-tag/peak_annotation/peak_list.txt`, `cut-tag/peak_annotation/peak_files/`, `cut-tag/peak_annotation/region_type.csv` and `tair12_genome/TAIR12.gtf` | `cut-tag/peak_annotation/*_peak_tair12_annotated.xlsx` |
| Gene-set overlap | `cut-tag/peak_annotation/*_peak_tair12_annotated.xlsx` | TSS gene lists in `cut-tag/peak_annotation/gene_set/` |
| Motif input | `cut-tag/peak_annotation/gene_set/3_overlap.csv`, `tair12_genome/TAIR12.gtf` and `tair12_genome/TAIR12.fna.gz` | `cut-tag/motif/3_overlap_promoter_genebody.fasta` |
| Genome-wide profiles | `cut-tag/metaplot/*_clean_log2_nosm.bedgraph` | Plots in the interactive R session |
| Peak overlap fractions | `cut-tag/peak_annotation/*_peak_tair12_annotated.xlsx` | CSVs in `cut-tag/peak_fraction/data/` and a generated text report in `cut-tag/peak_fraction/reports/` |
| TSS profiles and boxplots | `cut-tag/profile_boxplot/plot_bed/`, `cut-tag/profile_boxplot/matrix/` and `cut-tag/profile_boxplot/*_tss.csv` | See the TSS [input paths and workflow](cut-tag/profile_boxplot/README.md). |

When regenerating the three-way gene list, `venn.R` writes `3_overlap.csv` to the current R working directory. The motif script reads the copy in `cut-tag/peak_annotation/gene_set/`.

The shell scripts retain the original preprocessing and track-generation commands for reference. Running them requires separately supplied paired FASTQs or BAMs, working directories matching the commands, a genome index, chromosome sizes and a Chr2/Chr4 rDNA exclusion BED. Replace example sample names and server-specific reference paths for a local run. Use LF line endings when running the commands with Bash.

For tools that require an uncompressed FASTA, unpack the reference from the repository root.

```sh
gzip -dk tair12_genome/TAIR12.fna.gz
```

The motif extraction script can read the compressed FASTA directly.

## License

The analysis code and accompanying documentation are available under the MIT License. Reference data retain their original terms. [MIT License](LICENSE)
