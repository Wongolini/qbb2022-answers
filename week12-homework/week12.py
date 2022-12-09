#%%
import numpy as np
import matplotlib.pyplot as plt 


# %%

def WrightFisher(allele_freq, pop_size):
    p = allele_freq
    generations = []
    while p not in [0,1]:
        p = np.random.binomial(2*pop_size,p)/(2*pop_size)
        generations.append(p)
    return generations

def plotWrightFisher(allele_freq, pop_size):
    generations = WrightFisher(allele_freq, pop_size)
    plt.plot(generations)
    plt.xlabel('Generation')
    plt.ylabel('allele freq')
    plt.title(r'Wright Fisher $N={}, p_0={}$'.format(pop_size, allele_freq))
    plt.savefig('WrightFisher_N{}_p{}.png'.format(allele_freq, pop_size))

plotWrightFisher(.5,1000)
#%%
fixation_list = np.zeros((1000))
for iter in range(1000):
    fixation_list[iter] = len(WrightFisher(.5, 100))
plt.hist(fixation_list)
plt.title(r'Histogram: Time to Fixation $N=100,p_0=.5$')
plt.ylabel('Frequency')
plt.xlabel('Generations')
plt.savefig('Histogram_time_to_fixation.png')
# Implement a Wright-Fisher simulation of allele 
# frequencies for an arbitrary starting allele 
# frequency and population size.
#%%
pop_sizes = np.linspace(100,1e6,6)
fixation_time = []
for pop_size in pop_sizes:
    fixation_time.append(len(WrightFisher(.5, pop_size)))
plt.plot(np.log10(pop_sizes), np.log10(fixation_time), '-o')
plt.title('WrightFisher popsizes vs fixation time; p=.5, pop_sizes[100,1e6]')
plt.xlabel(r'$log_{10}$(popsizes)')
plt.ylabel(r'$log_{10}$(fixation_time)')
plt.tight_layout()
plt.savefig('WrightFisher_fixation_time_varied_popsizes.png')
# %%
pop_size = 100
allele_freqs = np.linspace(.01,.99,20)
num_iters = 100
output_mat = np.zeros((num_iters,(len(allele_freqs))))
# columns are allele_freq
# rows are nth iteration
for i,a in enumerate(allele_freqs):
    fixation = np.zeros((num_iters))
    for iter in range(num_iters):
        fixation[iter] = len(WrightFisher(a, pop_size))
    output_mat[:,i] = fixation
#%%
plt.violinplot(output_mat)
plt.xlabel('allele_freqs')
plt.ylabel('Time to Fixation (generations)')
plt.title('Variability of Allele Freq vs Time to Fixation')
plt.savefig('ViolinPlot_Allele_freq_vs_fixation.png')
# %%
