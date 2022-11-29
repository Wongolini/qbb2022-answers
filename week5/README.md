```
samtools view -q 10 D2_Sox2_R2.bam -o D2_Sox2_R2_filtered.bam
macs2 callpeak -t D2_Sox2_R1_filtered.bam -c D2_Sox2_R1_input.bam -f BAM -n Sox2_R1 -g 8.3e7 -B --outdir macs2_peaks_out/
macs2 callpeak -t D2_Sox2_R2_filtered.bam -c D2_Sox2_R2_input.bam -f BAM -n Sox2_R2 -g 8.3e7 -B --outdir macs2_peaks_out/
bedtools intersect -a D2_Klf4_peaks.bed -b macs2_peaks_out/Sox2_R1_R2_intersect.tsv | wc -l
42 lines in this intersection
582 Sox2_R1_R2_intersect.tsv
42/582 = .07216
mkdir bdg_files
cp all bdg files in bdg_files
for x in *; do python ../scale_bdg.py $x scale_$x; done
for x in scale*; do awk '{ if ($2 < 35507055 && $3 > 35502055) print $0 }' $x > cropped_$x; done
sort -k5nr Sox2_combined.narrowPeak > sorted_Sox2_combined.narrowPeak
samtools faidx mm10.fa -r macs2_peaks_out/faidx_sorted300_Sox2_combined.narrowPeak -o faidx_sorted300_Sox3_combined.fasta
meme-chip faidx_sorted300_Sox3_combined.fasta -maxw 7
tomtom combined.meme /Users/cmdb/Downloads/motif_databases/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme 
grep KLF4|SOX2 tomtom.tsv > MEME_extracted_SOX2_KLF4.tsv

```
