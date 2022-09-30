#%%

import matplotlib.pyplot as plt
import pandas as pd
from yaml import parse
from vcfParser import *
import numpy as np
#%%
def plot_stats(vcf_reader):
    # histogram for DP=depth --> rec.INFO['DP']
    # histogram for QUAL --> rec.QUAL
    # histogram for allele frequency --> rec.INFO['AF']
    # barplot for predicted effects -->
    pass

vcf_file = 'annotated_sacCer_variants.vcf'
reader = parse_vcf(vcf_file)

'''
  'GT': 'Genotype',
  'GQ': 'Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype',
  'GL': 'Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy',
  'DP': 'Read Depth',
  'AD': 'Number of observation for each allele',
  'RO': 'Reference allele observation count',
  'QR': 'Sum of quality of the reference observations',
  'AO': 'Alternate allele observation count',
  'QA': 'Sum of quality of the alternate observations',
  'MIN_DP': 'Minimum depth in gVCF output block.'
'''
samples = reader[0][9:19]
fields = ['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL']
sample_data = np.empty((len(samples),len(fields),len(reader)-1),dtype=object)
# 3d matrix structure
# rows = samples
# columns = fields
# depth = length of vcf file (variants)
#%%
#qualities = np.zeros((len(reader)-1))
for i,record in enumerate(reader):
    if i == 0:
        continue 
    #print(samples)
    #print(reader[i][8])
    data = np.array(record[9:19])
    sample_data[:,:,i-1] = data

class Vis:
    def __init__(self, data):
        self.data = data 
    
    def make_hist(self,field,title):
        # create 10 subplots for all samples
        # create 2 x 5 plots
        fig,axs = plt.subplots(nrows=2,ncols=5,figsize=[20,10])
        for s,ax in zip(samples,axs.ravel()):
            x = self.data[samples.index(s),fields.index(field),:]
            for i,v in enumerate(x):
                if v=='.':
                    x[i] = np.nan
                else:
                    x[i] = np.float32(v)
            ax.hist(x,density=True,bins=125)
            ax.tick_params(labelrotation=45)
            #for label in ax.xaxis.get_ticklabels()[::2]:
            #    label.set_visible(False)
            #ax.set_xlim(0,100)
            ax.set_title(s)
            ax.set_xlabel(field)
            ax.set_ylabel('frequency')
        fig.suptitle(title)
        plt.tight_layout()
        plt.show()
        plt.close()
'''

The read depth distribution of variant genotypes (histogram)
    This information can be found in the sample specific FORMAT field for each variant/line. Check the file header to decide which ID is appropriate.

The quality distribution of variant genotypes (histogram)
    This information can be found in the sample specific FORMAT field for each variant/line. Check the file header to decide which ID is appropriate.

The allele frequency spectrum of your identified variants (histogram)
    This information is pre-calculated for you and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate.

A summary of the predicted effect(s) of each variant as determined by snpEff (barplot)
    This information was added to the VCF by snpEff and can be found in the variant specific INFO field. Check the file header to decide which ID is appropriate and how to parse the information.
We encourage you to consider every possible effect for each variant, but feel free to just grab the first one.

'''

V=Vis(sample_data)
V.make_hist('DP','Depth') # sample specific FORMAT field
V.make_hist('QA','Quality ALT') # sample specific FORMAT field
af = [reader[i][7]['AF'] for i in range(1,len(reader))]
plt.figure(figsize=[10,10])
plt.title('Allele Frequency Spectrum')
plt.xlabel('Allele Frequency')
plt.ylabel("Frequency")
plt.hist(af,density=True,bins=55)
#%%
effects = {'LOF':0, 'NMD':0}
for i,r in enumerate(reader):
    if i == 0:
        continue
    k = list(r[7].keys()) 
    if 'LOF' in k:
        effects['LOF'] += 1 
    if 'NMD' in k:
        effects['NMD'] += 1
plt.title('Predicted Variants Effects')
plt.ylabel('Counts')
plt.bar(effects.keys(), effects.values(), color ='maroon',
        width = 0.4)
'''
'LOF': "Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'",
'NMD': "Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'"
'''
#V.make_hist('AO','Allele Frequency Spectrum') # variant specific INFO field
#%%


# %%
