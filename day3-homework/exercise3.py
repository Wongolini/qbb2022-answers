import numpy as np
import matplotlib.pyplot as plt 
import sys 

# population, superpopulation, and sex (3 separate plots)
def visualize_color(data, qualities):
    # how I would do it without numpy
    pc1=data['PC1']
    pc2=data['PC2']
    qual=list(set(data[qualities]))
    qual_dict = {q:[[],[]] for q in qual}
    for i in range(len(data)):
        #print(data[i])
        qual_dict[data[i][qualities]][0].append(pc1[i]) # ith data entry of quality for PC1
        qual_dict[data[i][qualities]][1].append(pc2[i]) # ith data entry of quality for PC2
        
    
    for q in qual_dict.keys():
        plt.scatter(qual_dict[q][0],qual_dict[q][1],alpha=.2,s=3,label=q)
    plt.legend(ncol=3,bbox_to_anchor=(.5,.4))
    plt.savefig('{}_PCA.png'.format(qualities))
    plt.show()

def visualize_color_2(data, qualities):
    # how I would do it with numpy
    meta = np.unique(data[qualities])
    for m in meta:
        indices = np.where(data[qualities]==m)[0]
        plt.scatter(data[indices]['PC1'], data[indices]['PC2'],label=m,alpha=.2,s=3)
    plt.legend(ncol=3, bbox_to_anchor=(.5,.4))
    plt.savefig('{}_PCA_numpy.png'.format(qualities))
    plt.show()

if __name__ == "__main__":
    data = np.genfromtxt(sys.argv[1], 
                         dtype=None, encoding='UTF-8',
                         names=['sample_n1','pop','super_pop',
                                'gender','sample_n2','PC1','PC2','PC3'])
    for qual in ['pop','super_pop','gender']:
        visualize_color(data,qual)
        visualize_color_2(data,qual)


