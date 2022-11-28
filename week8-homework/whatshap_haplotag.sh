regs=( "chr11:1900000:2800000" "chr14:100700000:100990000" "chr15:23600000:25900000" "chr20:58800000:58912000" )
#regs=( "one" "two" "three" "four" )
chrdirs=( "chr11_medaka_variant_out" "chr14_medaka_variant_out" "chr15_medaka_variant_out" "chr20_medaka_variant_out" )
for i in ${!regs[@]}; do
	chr=${chrdirs[i]}
	#echo ${regs[i]}
	
	 whatshap haplotag -o ${chr}_whatshap_haplotag.bam --output-haplotag-list ${chr}_whatshap_haplotag.list --reference=hg38.fa --regions ${regs[i]} ${chr}/round_0_hap_mixed_phased.vcf.gz nanopore/methylation.bam
done
