#!/usr/bin/env python
#%%
import sys
import seaborn as sns
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class ThreeSeq:
    def __init__(self, file, bin_fname, title):
        chrom = b'chr15'
        start = 11170245
        end = 12070245
        self.title = title
        self.data = numpy.loadtxt(file, dtype=numpy.dtype([
            ('F1', int), ('F2', int), ('score', float)]))
        
        self.frags = numpy.loadtxt(bin_fname, dtype=numpy.dtype([
            ('chr', 'S5'), ('start', int), ('end', int), ('bin', int)]))
        self.start_bin = self.frags['bin'][numpy.where((self.frags['chr'] == chrom) &
                                            (self.frags['start'] <= start) &
                                            (self.frags['end'] > start))[0][0]]
        self.end_bin = self.frags['bin'][numpy.where((self.frags['chr'] == chrom) &
                                        (self.frags['start'] <= end) &
                                        (self.frags['end'] > end))[0][0]] + 1
        keep_ind_data = [i for i,x in enumerate(self.data) if x[0]>=self.start_bin and x[1]<=self.end_bin]
        self.trimmed_data = self.data[keep_ind_data]
        
        self.trimmed_data['score'] = numpy.log10(.001+self.trimmed_data['score'])
        self.square_matrix()
        self.remove_dd_bg()
        self.smooth_matrix()

    def square_matrix(self):
        shape0 = min([min(self.trimmed_data['F1']), min(self.trimmed_data['F2'])])
        shape1 = max([max(self.trimmed_data['F1']), max(self.trimmed_data['F2'])])
        matrix = numpy.zeros((shape1-shape0+1, shape1-shape0+1),dtype=numpy.float16)
        for x in self.trimmed_data:
            s = x[0]-shape0
            e = x[1]-shape0
            score = x[2]
            matrix[s,e] = score
            matrix[e,s] = score
        self.matrix = matrix 
        
    
    def remove_dd_bg(self):
        N = self.matrix.shape[0]
        mat2 = numpy.copy(self.matrix)
        for i in range(N):
            bg = numpy.mean(self.matrix[numpy.arange(i, N), numpy.arange(N - i)])
            mat2[numpy.arange(i, N), numpy.arange(N - i)] -= bg
            if i > 0:
                mat2[numpy.arange(N - i), numpy.arange(i, N)] -= bg
        self.matrix_removed_dd_bg = mat2

    def smooth_matrix(self):
        N = self.matrix_removed_dd_bg.shape[0]
        invalid = numpy.where(self.matrix_removed_dd_bg[1:-1, 1:-1] == 0)
        nmat = numpy.zeros((N - 2, N - 2), float)
        for i in range(3):
            for j in range(3):
                nmat += self.matrix_removed_dd_bg[i:(N - 2 + i), j:(N - 2 + j)]
        nmat /= 9
        nmat[invalid] = 0
        self.final_matrix = nmat
    
    def plot(self):
        fig, ax = plt.subplots(figsize=[13,13])
        ax.set_title(self.title)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        sns.heatmap(self.final_matrix,xticklabels=range(self.start_bin, self.end_bin),
                        yticklabels=range(self.start_bin, self.end_bin),
                        cmap='seismic',ax=ax)
        for i,label in enumerate(ax.xaxis.get_ticklabels()):
            if i%10==0:
                label.set_visible(True)
            else:
                label.set_visible(False)
        for i,label in enumerate(ax.yaxis.get_ticklabels()):
            if i%10==0:
                label.set_visible(True)
            else:
                label.set_visible(False)
        plt.savefig('{}.png'.format(self.title))
        self.ax = ax



    #plt.savefig('{}.png'.format(self.title))

'''
Produce a horizontal 3-panel plot with a heatmap for ddCTCF, dCTCF, and dCTCF - ddCTCF
chr15:11170245-12070245

Filter out data with one or both interaction ends falling 
outside the desired bin range

Log-transform the scores (the dynamic range of data makes 
it hard to visualize the non-transformed data).
Also, shift the data by subtracting the minimum value so 
the new minimum value is zero (this will prevent issues where 
there is missing information)

Convert the sparse data into a square matrix (note that the 
sparse data only contains one entry per interaction with the 
lower-numbered bin in the first column). By converting the 
sparse matrix it into a complete matrix for plotting, you 
have two entries per interaction. For one line of the sparse 
data format, the data relates to the full matrix as follows:
mat[sparse['F1'][i], sparse['F2'][i]] = sparse['score'][i]
Plot the two matrices using the same maximum value 
(set vmax in imshow). I suggest using the magma color map, 
although you need to flip your scores to mimic the paper figure

CTCF = TF that helps with dna looping

'''
in1_fname = './analysis/hic_results/matrix/ddCTCF/iced/6400/ddCTCF_ontarget_6400_iced.matrix'
in2_fname = './analysis/hic_results/matrix/dCTCF/iced/6400/dCTCF_ontarget_6400_iced.matrix'
in3_fname = '/Users/cmdb/qbb2022-answers/week6-homework/matrix/dCTCF_full.6400.matrix'
in4_fname = '/Users/cmdb/qbb2022-answers/week6-homework/matrix/ddCTCF_full.6400.matrix'
bin_fname = './matrix/6400_bins.bed'

T1 = ThreeSeq(in1_fname, bin_fname, 'ddCTCF_ontarget_6400')
T1.plot()
T2 = ThreeSeq(in2_fname, bin_fname, 'dCTCF_ontarget_6400')
T2.plot()
T3 = ThreeSeq(in3_fname, bin_fname, 'dCTCF_full.6400')
T3.plot()
T4 = ThreeSeq(in4_fname, bin_fname, 'ddCTCF_full.6400')
T4.plot()



#%%
obs = [T1,T2,T3,T4]
fig, axs = plt.subplots(2,2,figsize=[13,13])
fig.suptitle('3D Genome ddCTCF-dCTCF 6400 full/on-target',fontsize=15)
for i,ax in enumerate(numpy.ravel(axs)):
    ob = obs[i]
    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    sns.heatmap(ob.final_matrix,xticklabels=range(ob.start_bin, ob.end_bin),
                    yticklabels=range(ob.start_bin, ob.end_bin),
                    cmap='seismic',ax=ax)
    for j,label in enumerate(ax.xaxis.get_ticklabels()):
        if j%10==0:
            label.set_visible(True)
        else:
            label.set_visible(False)
    for j,label in enumerate(ax.yaxis.get_ticklabels()):
        if j%10==0:
            label.set_visible(True)
        else:
            label.set_visible(False)
    ax.set_title(ob.title)
    plt.tight_layout()
plt.savefig('3Dgenome.png')

# %%
# dd-dc
TC = T2
TC.final_matrix = (T1.final_matrix-TC.final_matrix)
TC.title = 'ddCTCF - dCTCF target_6400'
TC.plot()

TD = T4
TD.final_matrix = (T3.final_matrix-TD.final_matrix)
TD.title = 'ddCTCF - dCTCF full_6400'
TD.plot()
# %%
