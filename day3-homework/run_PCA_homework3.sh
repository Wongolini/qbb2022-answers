plink --vcf ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out cichlids
plink --vcf ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract cichlids.prune.in --make-bed --pca 3 --out cichlids

./homework3.py
./exercise3.py cichlids.eigenvec

