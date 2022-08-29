
genefile=/Users/cmdb/data/bed_files/genes.bed
vcf_file=~/data/vcf_files/random_snippet.vcf
bedtools_out=./intersect.out
bedtools intersect -a $genefile -b $vcf_file > $bedtools_out
wc -l $bedtools_out

