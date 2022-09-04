#%%
import numpy as np
import csv
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import pandas as pd 

def get_header(file):
    with open(file) as f:
        header = f.readlines()[0].rstrip('\n').split(',')
    return header

def count_mutations_per_probad(data):
    # Proband_id
    # find unique proband_ids
    '''
    The Phase_combined column records the inferred parent of origin 
    of each de novo mutation. 
    Break the counts of de novo mutations down into maternally inherited 
    versus paternally inherited de novo mutations (ignore mutations of unknown parental origin). 

    Store these counts in a new file with 3 fields containing: 
        proband ID, number of paternally inherited de novo mutations for that proband, and number of
     maternally inherited de novo mutations for that proband.
    proband_id, # de novo muts, # maternal
    '''
    probands = np.unique(data['Proband_id'])
    mutations_per_proband = {p:np.where(p==data['Proband_id'])[0]
                            for p in probands}
    proband_data = []
    for p_list in mutations_per_proband.keys():

        proband_id_data = [p_list,len(np.where(p_list==data['Proband_id'])[0])]
        proband_slice = data[mutations_per_proband[p_list]]
        paternal_slice = proband_slice['Phase_combined']
        for parent in ['father','mother']:
           
            proband_id_data.append(len(np.where(parent==paternal_slice)[0]))
        proband_data.append(proband_id_data)
    return proband_data,('Proband_id','n_muts','n_father','n_mother')

def write_out(to_out, header):
    with open('day5_lunch.txt','w') as f:
        writer = csv.writer(f,delimiter=',')
        writer.writerow(header)
        for line in to_out:
            writer.writerow(line)

def merge_data(d1,d2):
    d1n = d1[np.argsort(d1['Proband_id'])]
    d2n = d2[np.argsort(d2['Proband_id'])]
    #assert np.all(d1n['Proband_id'] == d2n['Proband_id'])
    all_dat = np.zeros((len(d1n), 6),dtype=np.int32)
    
    all_dat[:,0] = d1n['Proband_id']
    all_dat[:,1] = d1n['n_muts']
    all_dat[:,2] = d1n['n_father']
    all_dat[:,3] = d1n['n_mother']
    all_dat[:,4] = d2n['Father_age']
    all_dat[:,5] = d2n['Mother_age']

    #return all_dat, ['Proband_id','n_muts','n_father','n_mother','Father_age','Mother_age']
    # I gave up and used dataframes from pandas :(
    # I couln't get structure array to work from scratch
    return pd.DataFrame(all_dat,columns=['Proband_id','n_muts','n_father','n_mother','Father_age','Mother_age'])


class Stats_Visualise:
    def __init__(self, df):
        self.df = df 
        self.record_stats={}
        

    def compare_hist(self,f1,f2):
        #f1 field 1
        #f2 field 2
        plt.title('{} vs {} hist'.format(f1,f2))
        plt.hist(self.df[f1],density=True,alpha=.25,label=f1)
        plt.hist(self.df[f2],density=True,alpha=.25,label=f2)
        plt.legend()
        plt.ylabel('frequency')
        
        plt.savefig('{}_{}_hist.png'.format(f1,f2))
        plt.show()
        plt.close()

    def cartesian(self,f1,f2):
        plt.title('{} vs {}'.format(f1,f2))
        plt.scatter(self.df[f1],self.df[f2])
        plt.xlabel(f1)
        plt.ylabel(f2)
        
        plt.savefig('{}_{}_scatter.png'.format(f1,f2))
        plt.show()
        plt.close()

    def OLS(self,f1,f2):
        '''
        Use ordinary least squares smf.ols() to test for an 
        association between maternal age and maternally inherited 
        de novo mutations.

        Is this relationship significant?
        What is the size of this relationship?
        '''
        model = smf.ols(formula = "{} ~ 1 + {}".format(f1,f2), data = self.df)
        results = model.fit()
        self.record_stats['{}_{}_ols_pval'.format(f1,f2)] = results.pvalues
        print('OLS on fields: {}, {}'.format(f1,f2))
        print(results.summary())


if __name__ == "__main__":
    #file1 = sys.argv[1] 
    #file2 = sys.argv[2]
    file1 = 'aau1043_dnm.csv'
    file2 = 'aau1043_parental_age.csv'
    data1 = np.genfromtxt(file1,delimiter=',',
                        encoding='utf-8',dtype=None, skip_header=1,
                        names=get_header(file1))
    proband_data,header=count_mutations_per_probad(data1)
    write_out(proband_data,header)

    data1n = np.genfromtxt('day5_lunch.txt',dtype=None,skip_header=1,
                            encoding="UTF-8",delimiter=',',
                            names=get_header('day5_lunch.txt'))

    data2 = np.genfromtxt(file2,dtype=None,skip_header=1,
                          encoding='utf-8',delimiter=',',
                          names=get_header(file2))

    all_data= merge_data(data1n, data2)
    V= Stats_Visualise(all_data)
    
    V.cartesian('Mother_age','n_muts')
    V.cartesian('Father_age','n_muts')
    V.OLS('Mother_age','n_muts',)
    V.OLS('Father_age','n_muts')
    V.compare_hist('n_father','n_mother')
   

