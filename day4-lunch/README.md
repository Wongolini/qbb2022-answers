# 1. Reproduce plots using bxlab/cmdb-plot-vcfs
## Portion of do_all.sh output that reports how many bp each feature covers
```
*** Subsetting .vcf for each feature
--- Subsetting exons.chr21.bed.vcf
    + Covering 1107407 bp
--- Subsetting processed_pseudogene.chr21.bed.vcf
    + Covering 956640 bp
--- Subsetting protein_coding.chr21.bed.vcf
    + Covering 13780687 bp
```
## Three other gene_type present in the GENCODE.gtf

gene_type "transcribed_unprocessed_pseudogene" - this is a psueodogene, but it is still transcribed. Does its RNA serve a function? The RNA is not processed, so it isn't an miRNA?

gene_type "lncRNA" - long non coding RNAs regulate gene expression at epigenetic/chromatin remodeling and trancsriptional levels. 

gene_type "miRNA" - micro RNAs are small fragments from hairpin loops that have homology to the 3'UTR of coding RNAs. They help modulate gene expression at the post-transcriptional level.

# 2. Modify Workflow
## Work inside cmdb-plot-vcfs direcotry to:
	For Plots:
		- log transform data <-- add all data with .001 to make all data nonzero
		- set same y-axis <-- use ax.set_ylabels()/ax.set_ylim()/ax.set_xlim()
		- add titles  <-- ax.set_title()
	Create lncRNA.chr21.bed.vcf.png plot

# Add documentation for bxlab/cmdb-plot-vcfs:

# Synopsis
Parse VCF files and plot the distribution of allele counts for a given gene_type; gene_type= protein_coding, lncRNA, miRNA,...

# Dependencies:
	1. bedtools
	2. python3
	3. matplotlib 

# Usage
1. Make a copy of code and input data
```
cd
git clone https://github.com/bxlab/cmdb-plot-vcfs
cd cmdb-plot-vcfs
cp ~/data/vcf_files/random_snippet.vcf .
cp ~/data/gtf_files/gencode.v41.annotation.gtf.gz .
gunzip *.gz
```
2. Create a new conda environment with required software
```
conda create --name day4-lunch
conda activate day4-lunch
conda install bedtools matplotlib
```

3. Run Workflow 
[do_all.sh] is a wrapper script that will run subset_regions.sh, bedtools and plot_vcf_ac.py to produce a series of .png files of the distribution of alleles.
```
bash do_all.sh <random_snippet.vcf> <gencode.v41.annotation.gtf>
```

# Example output
For chrom 21 and gene_type lncRNA, the following image file will be produced from do_all.sh
lncRNA.exons.chr21.bed.vcf.png



