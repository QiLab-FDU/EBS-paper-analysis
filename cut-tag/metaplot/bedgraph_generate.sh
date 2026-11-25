# Generate CPM-normalized coverage in 100 kb bins.
prefixes=("sample-1" "sample-2" "IgG" )
for prefix in "${prefixes[@]}"
do
bamCoverage -b ../bam/"${prefix}_final_clean.bam" -o "${prefix}_cpm_nosm.bedgraph" -p 24 -v --normalizeUsing CPM --binSize 100000 --outFileFormat bedgraph
done

bedtools unionbedg -i sample-1_cpm_nosm.bedgraph sample-2_cpm_nosm.bedgraph > merged.tmp
awk 'BEGIN {OFS="\t"} {
  sum = $4 + $5 + $6;
  avg = sum / 3;
  print $1, $2, $3, avg
}' merged.tmp > sample_cpm_nosm.bedgraph
rm merged.tmp

prefixes=("sample" "IgG" )
for prefix in "${prefixes[@]}"
do
bedGraphToBigWig "${prefix}_cpm_nosm.bedgraph" /public/home/xwzheng/ldn/tair12/TAIR12.sizes ../bw/"${prefix}_cpm_nosm.bw"
done

# Calculate log2 signal relative to the IgG control.
bigwigCompare -b1 ../bw/sample_cpm_nosm.bw -b2 ../bw/IgG_cpm_nosm.bw --operation log2 -o ../bw/sample_log2_nosm.bw -p 24 -v 

prefixes=("sample" )
for prefix in "${prefixes[@]}"
do
bigWigToBedGraph ../bw/"${prefix}_log2_nosm.bw" "${prefix}_log2_nosm.bedgraph"
done