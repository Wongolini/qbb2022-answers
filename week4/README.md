plink --vcf gwas_data/genotypes.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out cichlid
plink --vcf gwas_data/genotypes.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract cichlids.prune.in --make-bed --pca --out cichlids
python plot_pca.py  cichlids.eigenvec
plink --freqx -vcf ../gwas_data/genotypes.vcf
python ../manhatten_plotter.py GS451_IC50_gwas_results.qassoc ../gwas_data/genotypes.vcf ../gwas_data/GS451_IC50.txt 
echo "Causal gene on chr19:20372113 could be zinc-finger protein RefSeq:NR_036455.1, start:20370466 end:20399611."
# To figure out the reference if it wasn't known (which would be horrible) you could:
# 1. Email the author of the vcf file you are working with
# 2. Be more organized
# 3. Check if the VCF file tells you in the metadata
# 4. Painstakingly check the meta data for contigs to see if they match up with any of the publsihed UCSC genome browser reference genomes. And verify the REF positions and corresponding bp [A,T,G,C].
