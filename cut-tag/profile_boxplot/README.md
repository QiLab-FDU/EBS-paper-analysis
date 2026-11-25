# TSS signal profiles and boxplots

This workflow compares signal around transcription start sites (TSSs) for the gene sets used in Figures S9 and S12. See the [software versions and data sources](../../README.md).

## Workflow

The paths shown here are relative to the repository root. The BED extraction script reads gene lists from `cut-tag/peak_annotation/gene_set/` and reference intervals from `tair12_genome/tair12_final.bed`. It writes the profile BED files to `cut-tag/profile_boxplot/plot_bed/`. The included matrices are in `cut-tag/profile_boxplot/matrix/`, and the CSV files used by the R script are in `cut-tag/profile_boxplot/`.

1. Call peaks for each sample with MACS3 3.0.2 using `-g 142481245 -q 0.01 --keep-dup all --SPMR`.
2. Annotate peaks with ChIPseeker and retain euchromatic peaks. Use the gene sets exported by [`cut-tag/peak_annotation/venn.R`](../peak_annotation/venn.R) and convert them to TAIR12 BED intervals with [`cut-tag/profile_boxplot/plot_bed/tair12-extract-genes.py`](plot_bed/tair12-extract-genes.py).
3. Generate 10 bp coverage tracks with [`cut-tag/profile_boxplot/bw_generate.sh`](bw_generate.sh). Use [`cut-tag/profile_boxplot/deeptools.sh`](deeptools.sh) to build TSS-centered matrices with 100 bp bins over 1 kb on either side of the TSS.
4. Use the included CSV tables with [`cut-tag/profile_boxplot/boxplot.R`](boxplot.R) to average signal within 0.5 kb of each TSS and run one-sided Wilcoxon rank-sum tests.

When regenerating matrices from the original commands, use `plot_bed/EBS-FLAG_setdiff_SPO11-1-FLAG.bed` for the EBS-FLAG-only S9 TSS matrix. For the S12 combined profile, use the `computeMatrix` output `spo11-ebs-separate_wt-43_log2_tss_matrix.gz` as the `plotProfile -m` input.

## Included signal tables

| Figure | Signal table | Gene set |
| --- | --- | --- |
| S9 | [WT H3K4me3](wt-43_ebs-spo11-nooverlap_tss.csv) | EBS-FLAG-only genes |
| S9 | [WT H3K4me3](wt-43_spo11-ebs-nooverlap_tss.csv) | SPO11-1-FLAG-only genes |
| S9 | [WT H3K4me3](wt-43_spo11-ebs-overlap_tss.csv) | EBS-FLAG and SPO11-1-FLAG shared genes |
| S12 | [SPO11-1-FLAG in WT](spo11_alone-ebs_wt-43-overlap_tss.csv) | EBS-FLAG and WT H3K4me3 shared genes |
| S12 | [SPO11-1-FLAG in ebs](spo11-ebs_alone-ebs_wt-43-overlap_tss.csv) | EBS-FLAG and WT H3K4me3 shared genes |
