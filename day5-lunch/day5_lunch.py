#%%
import numpy as np
import csv
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import pandas as pd 
from scipy import stats


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


class Visualise:
    def __init__(self, df):
        self.df = df 

    def compare_hist(self,f1,f2,xlab):
        #f1 field 1
        #f2 field 2
        fig,ax = plt.subplots(figsize=[5,5])
        ax.set_title('{} vs {} hist'.format(f1,f2))
        ax.hist(self.df[f1],density=True,alpha=.25,label=f1)
        ax.hist(self.df[f2],density=True,alpha=.25,label=f2)
        ax.legend()
        ax.set_ylabel('frequency')
        ax.set_xlabel(xlab)
        #plt.savefig('{}_{}_hist.png'.format(f1,f2))
        #plt.show()
        self.hist_plot = ax

    def cartesian(self,f1,f2):
        fig,ax = plt.subplots(figsize=[5,10])
        ax.set_title('{} vs {}'.format(f1,f2))
        ax.scatter(self.df[f1],self.df[f2],alpha=.5)
        ax.set_xlabel(f1)
        ax.set_ylabel(f2)
        self.scatter_plot = ax 

class Stats:
    def __init__(self,all_data):
        self.df = all_data 
        self.record_stats={}

    def OLS(self,f1,f2):
        '''
        predict f2 from f1
        Use ordinary least squares smf.ols() to test for an 
        association between maternal age and maternally inherited 
        de novo mutations.

        Is this relationship significant?
        What is the size of this relationship?
        '''
        self.model = smf.ols(formula = "{} ~ 1 + {}".format(f2,f1), data = self.df).fit()
     
        self.record_stats['{}_{}_ols_pval'.format(f1,f2)] = self.model.pvalues
        V = Visualise(self.df)
        V.cartesian(f1,f2)        
        ax = V.scatter_plot
        x_theory = list(range(max(self.df[f1])))
        y_predict = [self.model.params[1]*x+self.model.params[0] for x in x_theory]
        ax.plot(x_theory, y_predict,'r--')
        print('OLS on fields: {}, {}'.format(f1,f2))
        #plt.show()
        plt.savefig('{}_{}_scatter.png'.format(f1,f2))
        print(self.model.summary())
    
    def OLS_predict(self,f1,theoretical_x):
        # takes model from OLS and produces regression line 
        # predicts only for a single variabless
        theoretical_data = pd.DataFrame(self.df.iloc[0]).transpose()
        for col in theoretical_data:
            theoretical_data[col] = 0
        theoretical_data[f1] = theoretical_x
        self.prediction=self.model.predict(theoretical_data)


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
 
    print('T-TEST')
    print(stats.ttest_ind(all_data["n_mother"],
                      all_data["n_father"]))

    S=Stats(all_data)
    S.OLS('Mother_age','n_muts',)
    print(S.model.params,'Mother_age, n_muts')
    print('Pval: ',S.model.pvalues)
    S.OLS('Father_age','n_muts')
    print(S.model.params,'Father_age, n_muts')

    print('Pval: ',S.model.pvalues)
    S.OLS_predict('Father_age',50.5)
    print('OLS MODEL PREDICTS \n Father_age={} --> n_muts={}'.format(50.5,S.prediction[0]))




# %%
