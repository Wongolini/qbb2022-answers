np.around(np.arange(.55,1.05,.05), decimals=2)[::-1] 
line creates a numpy array starting at .55 and ending at 1.05 incrementing by .05. Each entry is rounded to 2 decimal points. The [::-1] reverses the order of the array.

The experiment tests the hypothesis that the results of a 100 iterations of a series of coin tosses of n tosses 
follows a binomial distribution of p=.5. The assumption to use binom(p=.5,n=100) is beause if the coin is truly random it should follow a normal bell curve with E(x)=.5*n. However, we are testing it for a range of tosses (small to large sampling sizes) and a range of probabilities; hence the usage of the bonferroni power to assess 
the probability of a false positive/false negative resulting from the binomial test. We see that for small 
numbers of tosses power has a greater sensitivity to changing p, likewise, for a fixed p=.55, power has sensitivity
to changing n. 

If power is 1, we correctly rejected 100% the null hypothesis. We see that for a small n, p=.55 could be mistaken 
as a good fit for a binomial distribution. However, when increasing n, we see the probability of properly rejecting the 
null hypothesis increases. Similarily, when the p increases, the probability of properly rejecting the null hypothesis
increases despite n=10 because it becomes more unlikely for the assumed E(x)=.5*n as the skew grows sharper with increasing p. 
We see this pattern remains consistent inside the body of the heatmap.

Using multiple hypothesis testing streamlines the results due to the removal of "False positives" using the tailed p-val test.
This can be seen with the No_Multiple_Hypothesis_Testing image where the data shows a gradient along the axis of tosses and probs. 
Whereas with multiple hypothesis testing done by bonferroni, the results are more binary due to the conservative adjustments made 
to the pvalues. Thus we see that for any small sample size n=10, we almost always accept the null hypothesis even if p=1. We 
also see that even for large n, p=.55 fails to reject the null hypothesis.

HW4 notes

Genetic variation typically arrives from random events in meiosis, but certain recombination mechanisms (transmission distortion (TD)) have a selfish bias to propogate into the next generation despite providing negative fitness for the offspring.

"A Hidden Markov Model (HMM), with transition probabilities determined by rates of meiotic crossover and emission probabilities determined by rates of genotyping error, is then used to infer the most likely path along the phased haplotypes for each gamete (Fig. 1c), thereby imputing missing genotype data (Fig. 1d)."

"In order to benchmark rhapsodi’s performance, we developed a generative model to construct input data with varying sample sizes, rates of meiotic recombination, sequencing coverages, and genotyping error rates (Fig. S1)"

Uses hidden markov model with transition probabilities determined by rates of genotyping error. The model showed rigorous statistical significance in detecting genetic variation from various input data.

However, even with real world sperm input data, the model failed to detect evidence of TD in male meiotic cells. However the authors could not exclude TD as a possibility, just not within the 25 sperm donors used in the model.

A markov model is inherently the same as a coin flip model, just with multiple states with probabilities for each transition from state n to state m. Essentially it is similar to rolling a die, but instead of a 1/6 probability for each side, or even a single weight for each side, the die has conditional probabilities of rolling to a side given the die's prior state. Such a model can be visualized by fig.1.c. graph. A Hidden Markov Model was used to condense the model because an entire genome was considered, and thus many spatial coordinates required a repeating concurrent model. The number of times the sample is run is the number of sperm considered in the experiment. The number of iterations to detect TD or SNPs is the size of the genome (bp). The transition probabilities are determined by real data for genetic variation. A binomial test is used because truly random variation should follow a bell curve, but the data is discrete. If TD were to be true, there would be a skew in the observed experimental output that would not agree with the binomial test. The same is true for our experiment, explained above.


