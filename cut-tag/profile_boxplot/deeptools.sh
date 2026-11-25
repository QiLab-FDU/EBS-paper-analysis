# Fig. S9. WT H3K4me3 around TSSs of the EBS-FLAG and SPO11-1-FLAG gene sets.
computeMatrix scale-regions -b 2000 -a 2000 --regionBodyLength 2000 -R plot_bed/EBS-FLAG_intersect_SPO11-1-FLAG.bed plot_bed/SPO11-1-FLAG_setdiff_EBS-FLAG.bed plot_bed/EBS-FLAG_setdiff_SPO11-1-FLAG.bed -S WT-H3K4me3_log2_10bp.bw --samplesLabel "WT-H3K4me3" -out wt-43_spo11-ebs-region_matrix.gz -p 24  --skipZeros --binSize 100  --missingDataAsZero 

plotProfile -m wt-43_spo11-ebs-region_matrix.gz -o ../plot/wt-43_spo11-ebs_34-region_profile.png --dpi 300 --colors darkorchid darkslateblue deeppink  --plotHeight 8 --plotWidth 8 --legendLocation center-left --regionsLabel "1" "2" "3"

computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R plot_bed/EBS-FLAG_intersect_SPO11-1-FLAG.bed -S WT-H3K4me3_log2_10bp.bw --samplesLabel "WT-H3K4me3" -out wt-43_spo11-ebs-overlap_tss_matrix.gz -p 24  --skipZeros --binSize 100  --missingDataAsZero 
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R plot_bed/SPO11-1-FLAG_setdiff_EBS-FLAG.bed -S WT-H3K4me3_log2_10bp.bw --samplesLabel "WT-H3K4me3" -out wt-43_spo11-ebs-nooverlap_tss_matrix.gz -p 24  --skipZeros --binSize 100  --missingDataAsZero 
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R plot_bed/EBS-FLAG_intersect_SPO11-1-FLAG.bed -S WT-H3K4me3_log2_10bp.bw --samplesLabel "WT-H3K4me3" -out wt-43_ebs-spo11-nooverlap_tss_matrix.gz -p 24  --skipZeros --binSize 100  --missingDataAsZero 

# Fig. S12. SPO11 signal at genes shared by EBS-FLAG and WT H3K4me3.
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R plot_bed/EBS-FLAG_intersect_WT-H3K4me3.bed -S SPO11-1-FLAG_log2_10bp.bw ebs-SPO11-1-FLAG_log2_10bp.bw --samplesLabel "SPO11-1-FLAG" "ebs-SPO11-1-FLAG" -out spo11-ebs-separate_wt-43_log2_tss_matrix.gz -p 24  --skipZeros --binSize 100  --missingDataAsZero 

plotProfile -m spo11-ebs-34-separate_wt-43_log2_tss_matrix.gz -o ../plot/spo11-ebs-34-separate_wt-43_log2_tss_profile.png --dpi 300 --colors salmon darkslateblue --perGroup --plotHeight 8 --plotWidth 8 --legendLocation center-left

computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R plot_bed/EBS-FLAG_intersect_WT-H3K4me3.bed -S SPO11-1-FLAG_log2_10bp.bw --samplesLabel "SPO11-1-FLAG" -out spo11_alone-ebs_wt-43-overlap_tss_matrix.gz -p 24  --skipZeros --binSize 100  --missingDataAsZero 
computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R plot_bed/EBS-FLAG_intersect_WT-H3K4me3.bed -S ebs-SPO11-1-FLAG_log2_10bp.bw --samplesLabel "ebs-SPO11-1-FLAG" -out spo11-ebs_alone-ebs_wt-43-overlap_tss_matrix.gz -p 24  --skipZeros --binSize 100  --missingDataAsZero 