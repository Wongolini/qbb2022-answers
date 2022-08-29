 # QBB2022 - Day 1 - Lunch Exercises Submission

 1. I’m excited to learn <your_short_answer>.

2.B 
Average number of exons per gene: 62.342465753424655
I found this with "wc -l" for exons and genes, then divided the number genes by exons

2.C. Describe how you would find the median
I would write a python script that maps every exon beginning and end indices to the gene indices
and count for each gene how many exons exist. Then take the median of those counts.

3.A.
cut -f 4 chromHMM.E116_15_coreMarks_hg38lift_stateno.chr21.bed | sort -n | uniq -c

 305 1
 678 2
  79 3
 377 4
 808 5
 148 6
1050 7
 156 8
 654 9
  17 10
  17 11
  30 12
  62 13
 228 14
 992 15

3.B.
I would write a script. 
1. Define the unique states and use these unique states to define a dictionary in python
2. For each line, find the size of the genome fragment 
3. Use the state as a key for the dictionary and += the size


4.A.
grep AFR integrated_call_samples.panel | cut -f 2 | sort | uniq -c
 123 ACB
 112 ASW
 173 ESN
 180 GWD
 122 LWK
 128 MSL
 206 YRI

4.B.
for x in $(cut -f 3 integrated_call_samples.panel | grep -v "pop" | sort | uniq); do echo $x; grep $x integrated_call_samples.panel | cut -f 2 | sort | uniq -c; done
AFR
 123 ACB
 112 ASW
 173 ESN
 180 GWD
 122 LWK
 128 MSL
 206 YRI
AMR
 148 CLM
 107 MXL
 130 PEL
 150 PUR
EAS
 109 CDX
 142 CHB
 171 CHS
 127 JPT
 124 KHV
EUR
 184 CEU
 105 FIN
 107 GBR
 162 IBS
 112 TSI
SAS
 144 BEB
 113 GIH
 118 ITU
 158 PJL
 128 STU

5.B.
cut -f 1-9,13 random_snippet.vcf > HG00100.vcf

5.C.
cut -f 10 HG00100.vcf | grep -v '##' | grep -v "HG00100" | sort | uniq -c
9514 0|0
 127 0|1
 178 1|0
 181 1|1

5.D.
How many rows does AF=1 occur
cut -f 8 HG00100.vcf | grep ;AF=1; | wc -l
      15

5.E. 
grep -v ## HG00100.vcf | grep -v INFO | cut -f 8 | head -1 | grep -o AF= | uniq -c
   6 AF=


