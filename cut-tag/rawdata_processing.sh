
# Decompress tair12_genome/TAIR12.fna.gz before indexing the reference.
# The alignment command uses the Bowtie2 index path from the original server.
samtools faidx TAIR12.fna
cut -f1,2 TAIR12.fna.fai > TAIR12.sizes

bowtie2-build TAIR12.fna TAIR12  

prefixes=("sample")
for prefix in "${prefixes[@]}"
do
trim_galore --paired --quality 20 --length 20 ../raw_data/"${prefix}_1.fq.gz"       ../raw_data/"${prefix}_2.fq.gz" -j 8
cd ../bam
bowtie2 --very-sensitive --no-discordant --no-mixed --threads 30 -k 10 --xeq -x /public/home/xwzheng/ldn/tair12/TAIR12 -1 ../trim/"${prefix}_1_val_1.fq.gz" -2 ../trim/"${prefix}_2_val_2.fq.gz" -S ../sam/"$prefix.sam" &> ../txt/"${prefix}_bowtie2.txt"
samtools sort -@ 16 -o "${prefix}_sorted.bam" ../sam/"${prefix}.sam"
samtools flagstat -@ 16 "${prefix}_sorted.bam" > ../txt/"${prefix}_sorted.txt"
rm ../sam/"$prefix.sam"
# Keep alignments with at most two mismatches (XM tag).
samtools view -h "${prefix}_sorted.bam" | awk '/^@/ || (match($0,/XM:i:[0-9]+/) && substr($0,RSTART+5,RLENGTH-5)<=2)' | samtools view -bS > "${prefix}_mismatch2.bam"
samtools flagstat -@ 16 "${prefix}_mismatch2.bam" > ../txt/"${prefix}_mismatch2.txt"

# Retain unique alignments without an XS:i tag and with MAPQ >= 42.
samtools view -h "${prefix}_mismatch2.bam" | awk '$1~/^@/ || ($0!~/XS:i:/ && $5>=42)' | samtools view -b > "${prefix}_unique.bam"
samtools flagstat -@ 16 "${prefix}_unique.bam" > ../txt/"${prefix}_unique.txt"
samtools view -h "${prefix}_mismatch2.bam" | awk '$1~/^@/ || $0~/XS:i:/' | samtools view -b > "${prefix}_multimap.bam"
samtools flagstat -@ 16 "${prefix}_multimap.bam" > ../txt/"${prefix}_multimap.txt"

samtools sort -@ 16 "${prefix}_unique.bam" -o "${prefix}_unique_sorted.bam"
sambamba markdup -r "${prefix}_unique_sorted.bam" "${prefix}_unique_rmdup.bam"
samtools flagstat -@ 16 "${prefix}_unique_rmdup.bam" > ../txt/"${prefix}_unique_rmdup.txt"

samtools view -h "${prefix}_multimap.bam" | awk '$5 >= 10 || /^@/' | samtools view -bS -o "${prefix}_multimap_q10.bam"
samtools flagstat -@ 16 "${prefix}_multimap_q10.bam" > ../txt/"${prefix}_multimap_q10.txt"
samtools view -b -F 0x100 "${prefix}_multimap_q10.bam" > "${prefix}_multimap_bwbest.bam"
samtools flagstat -@ 16 "${prefix}_multimap_bwbest.bam" > ../txt/"${prefix}_multimap_bwbest.txt"
samtools sort -@ 16 "${prefix}_multimap_bwbest.bam" -o "${prefix}_multimap_bwbest_sorted.bam"
sambamba markdup -r "${prefix}_multimap_bwbest_sorted.bam" "${prefix}_multimap_rmdup.bam"
samtools flagstat -@ 16 "${prefix}_multimap_rmdup.bam" > ../txt/"${prefix}_multimap_rmdup.txt"
rm "${prefix}_multimap_q10.bam"
rm "${prefix}_multimap_bwbest.bam"

samtools merge -@ 24 "${prefix}_final_rmdup.bam" "${prefix}_unique_rmdup.bam" "${prefix}_multimap_rmdup.bam" 
samtools sort -@ 24 "${prefix}_final_rmdup.bam" -o "${prefix}_final_rmdup_sorted.bam"
samtools index "${prefix}_final_rmdup_sorted.bam"
samtools flagstat -@ 16 "${prefix}_final_rmdup_sorted.bam" > ../txt/"${prefix}_final_rmdup_sorted.txt"
# Exclude alignments overlapping the Chr2 and Chr4 rDNA regions.
samtools view -b "${prefix}_final_rmdup_sorted.bam" | bedtools intersect -v -a - -b ../../tair12/nor_chr2_4_rDNA.bed > "${prefix}_final_clean.bam"
samtools index "${prefix}_final_clean.bam"
rm "${prefix}_mismatch2.bam"
rm "${prefix}_unique.bam"
rm "${prefix}_unique_sorted.bam"
rm "${prefix}_unique_rmdup.bam"
rm "${prefix}_multimap.bam"
rm "${prefix}_multimap_rmdup.bam"
rm "${prefix}_multimap_bwbest_sorted.bam"
rm "${prefix}_final_rmdup.bam"
cd ../trim
done