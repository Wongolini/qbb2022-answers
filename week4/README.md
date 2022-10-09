plink --vcf gwas_data/genotypes.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out cichlid
plink --vcf gwas_data/genotypes.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract cichlids.prune.in --make-bed --pca --out cichlids
python plot_pca.py  cichlids.eigenvec
plink --freqx -vcf ../gwas_data/genotypes.vcf
python ../manhatten_plotter.py GS451_IC50_gwas_results.qassoc ../gwas_data/genotypes.vcf ../gwas_data/GS451_IC50.txt 
echo "Causal gene on chr19:20372113 could be zinc-finger protein RefSeq:NR_036455.1, start:20370466 end:20399611."