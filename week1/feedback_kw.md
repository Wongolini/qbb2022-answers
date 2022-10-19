# Week 1 Genome Assembly -- Feedback

1 + 1 + 1 + 1 + 1 + 1 + 1 + 1 + 0.25 + 0 = 8.25 points out of 10 possible points

1. Question 1.1, 1.4 how many reads (0.5 pts each)

  * good

2. Question 1.2, 1.4 simulation script(s)

  * good work on setting a random seed to make it reproducible
  * good use of functions
  * `random_sequence` doesn't need to be an array when you're just adding one and not noise. `genome_vec[begin:end] += 1` works just as well and should require less memory. Curious how specifically you would want to add noise to make it more realistic?


3. Question 1.2, 1.4 plotting script(s)

  * I like the double plotting of the poisson, both with the theoretical mu and the observed average. Very nice touch

4. Question 1.2, 1.4 histograms with overlaid Poisson distributions (0.5 pts each)

  * fantastic plots!

5. Question 1.3, 1.4 how much of genome not sequenced/comparison to Poisson expectations (0.5 pts each, 0.25 for answer and 0.25 for code)

  * you've gone above and beyond here, but specifically, we were asking how many genomic locations would you expect to have 0 coverage given the poisson distribution. You could find this from `y2[0]`.

6. Question 2 De novo assembly (0.5 pts each, 0.25 for answer and 0.25 for code)

  * number of contigs --> +0.5
  * total length of contigs --> +0.5

7. Question 2 De novo assembly cont (0.5 pts each, 0.25 for answer and 0.25 for code)

  * size of largest contig --> +0.5
  * contig n50 size --> +0.5, your answer is correct, but should be reported as the length of the contig, not the number of the contig.

8. whole genome alignment (0.33 pts each, 0.33/2 for answer and 0.33/2 for code)

  * average identity --> +0.33
  * length of longest alignment --> +0.33
  * how many insertions and deletions in assembly --> +0.33

9. decoding the insertion (0.5 pts each, 0.25 for answer and 0.25 for code)

  * position of insertion in assembly --> you have a typo here in your answer that you seem to correct when you run samtools, but even then, the insertion doesn't include 27498: "The insertion begins at 26788 and ends at 27998." --> +0.25
  * length of novel insertion --> didn't see this; +0

10. decoding the insertion cont (0.5 pts each, 0.25 for answer and 0.25 for code)

  * DNA sequence of encoded message
  * secret message --> I wouldn't expect you to get this secret message given your incorrect coordinates for samtools faidx; please show me that you did indeed get this message and I will give the points back.
