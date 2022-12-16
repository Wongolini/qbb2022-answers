## Week 6 -- 10 points possible

1 + 1 + 1 + 1 + 1 + 1 + 1 + 0.75 + 0.75 + 1 = 9.5 of 10 points possible

1. Given data question: What percentage of reads are valid interactions?

2. Given data question: What constitutes the majority of invalid 3C pairs?/What does it mean?

3. Script set up to (0.5 pts each)

  * read in data  
  * Filter data based on location  

4. Script set up to log transform the scores

5. Script set up to shift the data by subtracting minimum value

* I didn't see this. Please correct me if I'm wrong and missed it --> **fixed in regrade**

6. Script set up to convert sparse data into square matrix

7. Script set up to (0.33 pts each)

  * remove distance dependent signal
  * smooth
  * subtract

8. Turned in the plot of the 3 heatmaps (ddCTCF, dCTCF, and difference) for subset dataset (0.33 pts each ddCTCF/dCTCF/difference)

* please make a 3 panel plot rather than turning in separate plots

9. Turned in the plot of the 3 heatmaps (ddCTCF, dCTCF, and difference) for full dataset (0.33 pts each ddCTCF/dCTCF/difference)

* please make a 3 panel plot rather than turning in separate plots

10. Heatmap questions (0.33 pts each)

  * See the highlighted difference from the original figure?
  * impact of sequencing depth?
  * highlighted signal indicates?

Didn't see an answer to this question about heatmap interpretation; but I remember discussing this in our meeting, so I'm going to give you points for it

Possible Bonus points:

1. Insulation script set up to (0.25 pts each)

  * read in data
  * filter data based on location
  * log transform the data
  * shift the data by subtracting minimum value

2. Insulation script set up to (0.5 pts each)

  * convert sparse data into square matrix
  * find the insulation score by taking mean of 5x5 squares of interactions around target

3. Turned in the plot of the heatmap + insulation scores below (0.5 pts each panel)

**regrade 12/15:<br />
points for subtracting minimum score given back<br/>3D Genome plot you referenced has 6 panels, but repeats ddCTCF - dCTCF full_6400 rather than including ddCTCF_full.6400. It also repeats ddCTCF - dCTCF target_6400 rather than including dCTCF_ontarget_6400. I therefore stand by the grade where you did not turn in the multi-panel plots we asked for and you will not receive full points.**
