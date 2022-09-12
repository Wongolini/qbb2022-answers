#%%
from scipy.stats import poisson
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats

'''
Question 1.2. Write a program (in Python) to 
simulate sequencing 5x coverage of a 1Mbp genome 
with 100bp reads. The output of this simulation should 
be an array of length 1 million, where each element in 
the array is the coverage at that base (i.e. a count of 
the number of reads that overlapped that base’s position). 
You do not actually need to consider the sequence of the genome 
or the strand of the reads. Using this array, plot a histogram of 
the coverage. Then, overlay the histogram with a Poisson distribution
with lambda=5.

'''
import random

random.seed(25)
#%%
def simulate_sequencing(coverage):
    genome_length = int(1e6)
    read_length = 100
    n_reads = int((coverage/100)*genome_length)
    genome_vec = np.zeros((genome_length),dtype=np.int32)
    fragment_size = 100
    for i in range(n_reads): # iterate thru number of reads s.t. coverage=5
        begin=np.random.randint(0,genome_length-100) #randomly choose a place in the genome
        end=begin+100 
        random_sequence = np.ones(fragment_size,dtype=np.int32) # random array of hits with length of fragment size
        # one could add noise to random_sequence to make the simulation more realistic
        genome_vec[begin:end]+=random_sequence # add this array to the genome_vec
    return genome_vec

# Q 1.2
def plot_simulation(genome_vec,avg):
    x = np.arange(0,max(genome_vec),.001)
    y1 = poisson.pmf(x, mu=np.average(genome_vec))
    y2 = poisson.pmf(x, mu=avg)
    plt.title('Plot genome coverage simulation')
    plt.plot(x,y1,label='mu=empirical lambda {}'.format(np.average(genome_vec)))
    plt.plot(x,y2,'r--',label='mu={}'.format(avg),alpha=.5)
    plt.hist(genome_vec,density=True,bins=15,alpha=.3,label='aligned_hits_freq')
    plt.xlabel('Number of hits')
    plt.ylabel('Frequency')
    plt.legend()
    plt.show()
    plt.close()

def compute_0_cov(genome_vec):
    num_coverage_0 = len(np.where(genome_vec==0)[0])
    print('# bps with 0x coverage: {}'.format(num_coverage_0))


def QQ(genome_vec,avg):
    summary=stats.probplot(genome_vec, dist='poisson', sparams=(avg,), plot=plt, rvalue=True)
    slope, intercept, r_value, p_value, std_err = stats.linregress(summary[0][0], summary[0][1])
    plt.title('QQ plot of simulated sequencing data vs Poisson(mu={})'.format(avg))
    plt.show()
    plt.close()
    print('y={}x+{}'.format(slope,intercept))
    return r_value
#%%
genome_vec = simulate_sequencing(5)
plot_simulation(genome_vec, 5)

# %%
# Q 1.3
compute_0_cov(genome_vec)
# %%
#Q 1.4
'''
Question 1.4. Now repeat the analysis with 15x coverage: 
simulate the appropriate number of reads and compute coverage,
make a histogram, 
overlay a Poisson distribution with lambda=15,
compute the number of bases with 0x coverage, and
evaluate how well it matches the Poisson expectation.
'''
genome_vec = simulate_sequencing(15)
plot_simulation(genome_vec, 15)
compute_0_cov(genome_vec)
r_value = QQ(genome_vec,15)
print('The simulation data fits a poisson distribution with lambda={} with r^2={}'.format(15,r_value))
# %%
