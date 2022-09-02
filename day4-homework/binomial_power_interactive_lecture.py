#!/usr/bin/env python
#%%
import numpy as np
from scipy.stats import binomtest
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

def simulate_coin_toss(n_tosses, prob_heads = 0.5, seed=None):
    '''
    Simulates a coin toss (fair or unfair depending on prob_heads)
    Input: n_tosses, an integer, number of coin tosses to simulate
           prob_heads, a float, the probability the coin toss will return heads; default is 0.5, meaning a fair coin
           seed, an integer, the random seed that will be used in the simulation; default is None
    Output: results_arr, an array of 1's and 0's with 1's representing heads and 0's representing tails
    Resources: https://np.org/doc/stable/reference/random/generated/np.random.choice.html
    '''
    if seed is not None:
        np.random.seed(seed)
    results_arr = np.random.choice([0,1], size=n_tosses, p = [1-prob_heads, prob_heads])
    return (results_arr)

def perform_hypothesis_test(n_heads, n_tosses):
    '''
    Performs a two-sided binomial test
    Input: n_heads, an integer, number of coin tosses that resulted in heads/success
           n_tosses, an integer, total number of coin tosses/trials
    Output: pval, a float reporting the final pvalue from the binomial test
    Resources: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binomtest.html#scipy.stats.binomtest
    '''
    binom_result = binomtest(n_heads, n_tosses) 
    # The binomial test [1] is a test of the null hypothesis 
    # that the probability of success in a Bernoulli experiment is p.
    # by default binomtest uses p=.5
    pval = binom_result.pvalue
    return(pval)

def correct_pvalues(pvals):
    '''
    Will apply the bonferroni multiple hypothesis testing correction method to input pvalues
    Input: pvals, an array-like object containing uncorrected pvalues
    Output: corrected_pvalues[1], an array containing corrected pvalues
    Resources: https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    '''
    corrected_pvalues = multipletests(pvals, method='bonferroni')
    return(corrected_pvalues[1])

def interpret_pvalues(pvals):
    '''
    Will interpret or convert pvalues from floats to booleans (Trues or Falses) reporting whether or not you reject the null hypothesis
    True -- reject null Hypothesis
    False -- fail to reject null hypothesis
    Input: pvals, an array-like object containing pvalues (floats)
    Output: interpreted, an array containing booleans
    '''
    interpreted = np.array(pvals) < 0.05
    return (interpreted)

def compute_power(n_rejected_correctly, n_tests):
    '''
    Will compute the power, defined as the number of correctly rejected null hypothesis divided by 
    the total number of tests (AKA the True Positive Rate or the probability of detecting something 
    if it's there)
    Input: n_rejected_correctly, an integer, the total number of tests in which the null hypothesis was correctly rejected
           n_tests, an integer, the total number of hypothesis tests which were performed
    output: power, a float, the power
    '''
    power = n_rejected_correctly / n_tests
    return(power)

def run_experiment(prob_heads, n_toss, n_iters = 100, seed = 389, correct_the_pvalues = False):
    '''
    Input: prob_heads, a float, the probability of a simulated coin toss returning heads
           n_toss, an integer, the number of coin tosses to simulate
           n_iters, an integer, the number of iterations for each simulation. default is 100
           seed, an integer, the random seed that you will set for the whole simulation experiment. default is 389
           correct_the_pvalues, a boolean, whether or not to perform multiple hypothesis testing correction
    Output: power, a float, the power of the experiment
    '''
    np.random.seed(seed)
    pvals = []
    for k in range(n_iters):
        results_arr = simulate_coin_toss(n_toss, prob_heads = prob_heads)
        n_success = np.sum(results_arr)
        pvals.append(perform_hypothesis_test(n_success, n_toss))
    if correct_the_pvalues:
        pvals = correct_pvalues(pvals)
    pvals_translated_to_bools = interpret_pvalues(pvals) # array of T/F
    power = compute_power(np.sum(pvals_translated_to_bools), n_iters) 
    return(power)
#%%
def heatmap(results, xlab, ylab,title):
    fig,ax = plt.subplots(figsize=(8,7))
    ax = sns.heatmap(results,annot=True,
                     cmap=plt.get_cmap('cividis'),
                     cbar_kws={'label': 'power'})
                     #cividis is a colorblind friendly palette
    ax.set_xlabel('probs')
    ax.set_ylabel('tosses')
    ax.set_xticklabels(xlab)
    ax.set_yticklabels(ylab)
    ax.set_title('Simulation Power Heat Map ' + title,
                 fontsize=13)
    #plt.show()
    plt.savefig(title+'.png')
    plt.close()
#%%
def simulate(title, multiple_hyp_test):
    results = np.zeros((num_tosses.shape[0], probs.shape[0])) #tosess-> rows, probs-->columns
    for i,n_toss in enumerate(num_tosses):
        for j,p in enumerate(probs):
            if multiple_hyp_test==True:
                results[i,j]=run_experiment(p, n_toss, correct_the_pvalues=True)
            else:
                results[i,j]=run_experiment(p, n_toss, correct_the_pvalues=False)
    heatmap(results, probs, num_tosses, title)

if __name__ == "__main__":
    '''
    Make a plot that displays power, computed without 
    Multiple Hypothesis Testing Correction and a plot 
    that displays power with Multiple Hypothesis Testing Correction. 
    It is up to you whether these are two separate plots or whether 
    you have one plot with two columns/panels

    The power of a hypothesis test is between 0 and 1; if the power is 
    close to 1, the hypothesis test is very good at detecting a false null hypothesis

    As the power increases, there is a decreasing probability of a type II error, 
    also called the false negative rate (β) since the power is equal to 1 - β.
    '''
    #Power is the probability of rejecting the null hypothesis when, in fact, it is false.
    num_tosses = np.array([10, 50, 100, 250, 500, 1000])
    probs = np.around(np.arange(0.55, 1.05, 0.05), decimals=2)[::-1]
    simulate('Multiple_Hypothesis_Testing', True)
    simulate('No_Multiple_Hypothesis_Testing', False)


# %%
