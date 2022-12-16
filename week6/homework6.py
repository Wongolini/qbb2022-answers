#%%
import numpy
import sys

def remove_dd_bg(mat):
    N = mat.shape[0]
    mat2 = numpy.copy(mat)
    for i in range(N):
        bg = numpy.mean(mat[numpy.arange(i, N), numpy.arange(N - i)])
        mat2[numpy.arange(i, N), numpy.arange(N - i)] -= bg
        if i > 0:
            mat2[numpy.arange(N - i), numpy.arange(i, N)] -= bg
    return mat2

def smooth_matrix(mat):
    N = mat.shape[0]
    invalid = numpy.where(mat[1:-1, 1:-1] == 0)
    nmat = numpy.zeros((N - 2, N - 2), float)
    for i in range(3):
        for j in range(3):
            nmat += mat[i:(N - 2 + i), j:(N - 2 + j)]
    nmat /= 9
    nmat[invalid] = 0
    return nmat


file = '/Users/cmdb/qbb2022-answers/week6-homework/analysis/hic_results/matrix/dCTCF/raw/6400/dCTCF_ontarget_6400.matrix'

#mat columns: bin# first part of intr; bin# second part of intr; norm score of intr
mat = numpy.genfromtxt(file)
print(mat)
#remove_dd_bg(mat)
'''
Produce a horizontal 3-panel plot with a heatmap for ddCTCF, dCTCF, and dCTCF - ddCTCF
chr15:11170245-12070245

You will want to input the sparse format that is provided as 
results from HiCPro. The sparse format data you want should 
come from the iced or normalized data folder. For the input bin 
file, you can use either 6400bp bin file from the raw folder in 
hic_results/matrix/XXXX.

To create the plot, you will need to do the following:

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

For the difference plot, I suggest using the seismic color map 
and norm=colors.CenteredNorm. It helps to remove the distance 
dependent signal and smooth the data first as there is noise. 
You can use the following function to remove the distance dependent signal:
'''
# %%
