#use flanking TSS upstream of TSS
awk -F '\t' '{if ($4==2) print}' ~/data/bed_files/chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21.bed > chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21_promoters.bed
bedtools intersect -b chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21_promoters.bed -a ~/data/vcf_files/random_snippet.vcf > putative_promoters.vcf
echo "Most likely SNP in place of Cytosine"
grep -v "#" putative_promoters.vcf | awk '{if ($4=="C") {print $5}}' | sort | uniq -c

