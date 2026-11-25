# Generate CPM-normalized coverage in 10 bp bins.
prefixes=("WT-H3K4me3-1" "WT-H3K4me3-2" "WT-IgG" "SPO11-1-FLAG-1" "SPO11-1-FLAG-2" "SPO11-1-FLAG-IgG" "ebs-SPO11-1-FLAG-1" "ebs-SPO11-1-FLAG-2" "ebs-SPO11-1-FLAG-IgG")
for prefix in "${prefixes[@]}"
do
bamCoverage -b ../bam/"${prefix}_final_clean.bam" -o ../bedgraph/"${prefix}_10bp.bedgraph" -p 30 --normalizeUsing CPM --centerReads --extendReads --binSize 10 --smoothLength 30 --outFileFormat bedgraph
done


# Merge replicates from the directory containing the bedGraph files.
bedtools unionbedg -i WT-H3K4me3-1_10bp.bedgraph WT-H3K4me3-2_10bp.bedgraph > merged.tmp
awk 'BEGIN {OFS="\t"} {
  sum = $4 + $5 + $6;
  avg = sum / 3;
  print $1, $2, $3, avg
}' merged.tmp > WT-H3K4me3_10bp.bedgraph
rm merged.tmp

bedtools unionbedg -i SPO11-1-FLAG-1_10bp.bedgraph SPO11-1-FLAG-2_10bp.bedgraph > merged.tmp
awk 'BEGIN {OFS="\t"} {
  sum = $4 + $5 + $6;
  avg = sum / 3;
  print $1, $2, $3, avg
}' merged.tmp > SPO11-1-FLAG_10bp.bedgraph
rm merged.tmp

bedtools unionbedg -i ebs-SPO11-1-FLAG-1_10bp.bedgraph ebs-SPO11-1-FLAG-2_10bp.bedgraph > merged.tmp
awk 'BEGIN {OFS="\t"} {
  sum = $4 + $5 + $6;
  avg = sum / 3;
  print $1, $2, $3, avg
}' merged.tmp > ebs-SPO11-1-FLAG_10bp.bedgraph
rm merged.tmp


prefixes=("WT-H3K4me3" "SPO11-1-FLAG" "ebs-SPO11-1-FLAG" "WT-IgG" "SPO11-1-FLAG-IgG" "ebs-SPO11-1-FLAG-IgG")
for prefix in "${prefixes[@]}"
do
bedGraphToBigWig "${prefix}_10bp.bedgraph" tair12/TAIR12.sizes ../bw/"${prefix}_10bp.bw"
done


# Calculate log2 signal relative to each matched IgG control.
bigwigCompare -b1 ../bw/WT-H3K4me3_10bp.bw -b2 ../bw/WT-IgG_10bp.bw --operation log2 -o ../bw/WT-H3K4me3_log2_10bp.bw -p 24 -v --binSize 10
bigwigCompare -b1 ../bw/SPO11-1-FLAG_10bp.bw -b2 ../bw/SPO11-1-FLAG-IgG_10bp.bw --operation log2 -o ../bw/SPO11-1-FLAG_log2_10bp.bw -p 24 -v --binSize 10
bigwigCompare -b1 ../bw/ebs-SPO11-1-FLAG_10bp.bw -b2 ../bw/ebs-SPO11-1-FLAG-IgG_10bp.bw --operation log2 -o ../bw/ebs-SPO11-1-FLAG_log2_10bp.bw -p 24 -v --binSize 10
