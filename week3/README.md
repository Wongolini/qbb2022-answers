To run bwa mem
bwa mem -t 4 -R "@RG\tID:A0_11\tSM:A0_11" reference_genome/sacCer3.fa fastq/A01_11.fastq
run_bwa_mem.sh
run_samtools_sort.sh
for x in ./*.bam; do samtools index -b $x $x.bai; done
ls *.bam > samples_bam.txt
freebayes -L samples_bam.txt -p 1 -f reference_genome/sacCer3.fa > sacCer_variants.vcf 
vcffilter -f "QUAL > 20" sacCer_variants.vcf > filtered_sacCer_variants.vcf
snpeff ann R64-1-1.99 filtered_sacCer_primitives.vcf > annotated_sacCer_variants.vcf

