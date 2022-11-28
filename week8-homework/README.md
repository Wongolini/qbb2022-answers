hg38.fa is the human ref genome
nanopore data in nanopore

Part 1: Call and phase variant for each region

To call variants, you will be using medaka_variant. Unfortunately, like many newer tools medaka_variant does not have good documentation. In fact, I have found the help statement (argument -h) to be the most helpful.

To call variants, Medaka needs to know which model was used to do the basecalling, both for variant calling and phasing. You should use the default for both -s and -m. You will also need the phased vcf file as output (look at -p). Only one region can be specified for Medaka at a time so you will need to generate a phased vcf file for each region in regions.bed using the format chr:start-end to specify the region.

'''
run_medaka_variant.sh
whatshap_haplotag.sh
for x in ./*.bam; do
	samtools index $x;
done

whatshap split chr11_medaka_variant_out_whatshap_haplotag.bam chr11_medaka_variant_out_whatshap_haplotag.list --output-h1 chr11_h1.bam --output-h2 chr11_h2.bam 

whatshap split chr14_medaka_variant_out_whatshap_haplotag.bam chr14_medaka_variant_out_whatshap_haplotag.list --output-h1 chr14_h1.bam --output-h2 chr14_h2.bam 

whatshap split chr15_medaka_variant_out_whatshap_haplotag.bam chr15_medaka_variant_out_whatshap_haplotag.list --output-h1 chr15_h1.bam --output-h2 chr15_h2.bam 

whatshap split chr20_medaka_variant_out_whatshap_haplotag.bam chr20_medaka_variant_out_whatshap_haplotag.list --output-h1 chr20_h1.bam --output-h2 chr20_h2.bam 

'''