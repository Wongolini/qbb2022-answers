## Week 5 -- 10 points possible

0.66 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 1  = 9.66 points out of 10

1. Filter reads, Call peaks, and intersect peaks across Sox2 replicates (0.33 points each)

* should filter the reads for the input bam files as well; only see one `samtools view` and there should be more than 1 if you filtered the R1 and the R2.
* the genome size for mac2 callpeak should be different, [as seen here for chr17](https://github.com/igvteam/igv/blob/master/genomes/sizes/mm10.chrom.sizes)
* how did you make `macs2_peaks_out/Sox2_R1_R2_intersect.tsv` ? --> -0.33

2. Find the number of total peaks and overlapping peaks for Klf4 and Sox2 (0.5 for commands, 0.5 for result)



3. scale bedgraph files (4 different datasets, 0.25 each)

4. crop bedgraph files (4 different datasets, 0.25 each)

5. python script for plotting


6. 4 panel plot of read pile ups

* good plot overall, especially good job of using a single maximum for the y-axis limit,  but I'm not sure that the plot title makes sense

7. motif finding sort intersected sox2 replicate narrow peak by 5th columm, keep first 300 lines, awk command for reformatting (0.33 each)

* didn't keep the first 300 lines; was this a purposeful choice?

8. use samtools faidx to extract peak sequences and meme-chip to perform motif finding (0.5 each)

9. download and unpack motif database

* code for having done this? I see that you have based off of the tomtom command; just please record such things

10. match profiles from tomtom for klf4 and sox2 (0.5 for commands, 0.5 for result)
