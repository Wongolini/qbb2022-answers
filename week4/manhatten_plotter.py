import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import sys 
from vcfParser import *
# P = p-val
def manhattan(df):
    num_plots = len(set(df['CHR']))
    fig,axs = plt.subplots(ncols=num_plots,sharey=True)
    axs[0].invert_yaxis()
    axs[0].set_ylabel('log(p)')
  
    fig.suptitle('MANHATTAN PLOT {}'.format(qaassoc_file))
    for i,ax in enumerate(axs.ravel()):
        
        if i == 0:
            ax.spines['right'].set_visible(False)
     
        if i>0 and i<num_plots:
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().set_ticks([])

        if i == num_plots-1:
            ax.spines['left'].set_visible(False)
            ax.spines['right'].set_visible(True)
            ax.get_xaxis().set_ticks([])

        p_df = df[df['CHR']==i+1]
        ax.scatter(p_df[p_df['P']>=np.log10(1e-5)]['BP'],p_df[p_df['P']>=np.log10(1e-5)]['P'],alpha=.8,s=.5)
        ax.scatter(p_df[p_df['P']<np.log10(1e-5)]['BP'],p_df[p_df['P']<np.log10(1e-5)]['P'],color='red',alpha=.8,s=.5)
        ax.axhline(y=np.log10(1e-5),xmin=0,xmax=max(p_df['BP']),c='orange',linewidth=.5)
        ax.set_xlabel('Chr{}'.format(i+1))
        ax.get_xaxis().set_ticks([])
     
    plt.subplots_adjust(wspace=0.1,hspace=0.1)
    plt.show()

def boxplot_phenotype(df,vcf_file,ic50_file):
    #find top associated SNP
    ic50_df = pd.read_csv(ic50_file,sep='\t')
    ic50_df.index = ['{}_{}'.format(x,x) for x in ic50_df.FID]
    ic50_df = ic50_df.drop(['FID','IID'],axis=1)

    best = np.argmin(df['P'])
    slice = df.iloc[best]
    snp = slice['SNP']
    chr = slice['CHR']
    bp = slice['BP']
    p = slice['P']
    vcf = parse_vcf(vcf_file)
    for i,record in enumerate(vcf):
        if i == 0:
            header = record
        if record[2]==snp:
            store = record
            break 
    samples = header[9:]
    gen = store[9:]
    genotype_counts={'het':[],'hom0':[],'hom1':[]}
    for i,g in enumerate(gen):
        if '.' not in g:
            g_=g.split('/')
            g_=list(map(int, g_))
            if sum(g_) == 1:
                genotype_counts['het'].append(samples[i])
            elif sum(g_) == 0:
                genotype_counts['hom0'].append(samples[i])
            elif sum(g_) == 2:
                genotype_counts['hom1'].append(samples[i])
            else:
                pass
                

    #print(ic50_df.loc['1001_1001'])
    het = ic50_df.loc[genotype_counts['het']].to_numpy().flatten()
    het = het[~np.isnan(het)]
    hom0 = ic50_df.loc[genotype_counts['hom0']].to_numpy().flatten()
    hom0 = hom0[~np.isnan(hom0)]
    hom1 = ic50_df.loc[genotype_counts['hom1']].to_numpy().flatten()
    hom1 = hom1[~np.isnan(hom1)]
    all_genotypes = np.array([het,hom0,hom1],dtype=object)
    plt.boxplot(all_genotypes)
    plt.title('{} EFFECT SIZE BOXPLOT {}'.format(snp, ic50_file))
    plt.xticks([1,2,3], ['het','hom 0/0', 'hom 1/1'])
    plt.savefig('GWAS_effect_size.png')



if __name__ == "__main__":
    qaassoc_file = sys.argv[1]
    vcf_file = sys.argv[2]
    ic50_file = sys.argv[3]
    #df = pd.read_csv(file,sep='\t')
    #names = np.loadtxt(file)[0]
    data = pd.read_csv(qaassoc_file,delim_whitespace=True)
    data.sort_values(by=['CHR','BP'])
    data['P'] = np.log10(data['P'])
    manhattan(data)
    boxplot_phenotype(data, vcf_file, ic50_file)

