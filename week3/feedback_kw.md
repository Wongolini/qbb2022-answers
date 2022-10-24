# Week 3 Variant Calling -- Feedback

0 + 1 + 1 + 0.75 + 1 + 0 + 0.75 + 1 + 0.75 + 1 = 7.25 points out of 10 possible points

1. Index genome

  * what code did you use for this? --> +0

2. align reads

  * --> +1

3. sort bam files and index sorted bam files (0.5 points each)

  * --> +1

4. variant call with freebayes

  * should include the `--genotype-qualities` flag since we're asking you to plot a histogram of the GQ (variant genotype qualities)
  * --> +0.75

5. filter variants

  * --> +1

6. decompose complex haplotypes

  * --> +0
  * not listed in your README and the uploaded 1000 line vcf does not have decomposed complex haplotypes (e.g., there are variants with spans of nucleotides)

7. variant effect prediction

  * --> +0.75
  * your README specifies that you used a file that I imagine was from the vcfallelicprimitives result, but given that I have no record of you running that, the recorded command in the README is not consistent.

8. python plotting script

  * you want the GQ not the QA for the format field
  * good script overall --> +1

9. 4 panel plot (0.25 points each panel)

  * instead of 4 separate plots, you want a multipanel plot. You can make one of these in python by setting the number of rows and number of columns when you call `plt.subplots()` --> +0.75
  * interesting choice to separate out the samples into multiple panels for read depth and genotype quality. You were just supposed to plot a single pooled distribution for each. If you still want to separate them out by sample, you could utilize color and semi-transparency to separate out the sample distributions by ID within the same panel, but not necessary to do so.

10. 1000 line vcf

  * --> +1
