#%%
import scanpy as sc
import matplotlib.pyplot as plt
# Read 10x dataset
adata = sc.read_10x_h5("neuron_10k_v3_filtered_feature_bc_matrix.h5")
# Make variable names (in this case the genes) unique
adata.var_names_make_unique()
#%%
sc.pp.pca(adata)
sc.pl.pca(adata,save=True)
# %%
# step 1 filtering
sc.pp.filter_genes(adata, min_counts=1)         # only consider genes with more than 1 count
sc.pp.normalize_per_cell(                       # normalize with total UMI count per cell
     adata, key_n_counts='n_counts_all'
)
filter_result = sc.pp.filter_genes_dispersion(  # select highly-variable genes
    adata.X, flavor='cell_ranger', n_top_genes=1000, log=False
)
adata = adata[:, filter_result.gene_subset]     # subset the genes
sc.pp.normalize_per_cell(adata)                 # renormalize after filtering

sc.pp.log1p(adata)         # log transform: adata.X = log(adata.X + 1)
sc.pp.scale(adata)         # scale to unit variance and shift to zero mean
#%%
# dim reduction pre-processing PCA
sc.pp.pca(adata)
sc.pl.pca(adata,save=True)
# %%
# Clustering
sc.pp.neighbors(adata)
# Use leiden clustering to identify clusters in the data
# Produce t-SNE and UMAP plots showing the clusters
# Note that to create the UMAP transformation, you need to specify maxiter
# I suggest 1000

sc.tl.leiden(adata,key_added='leiden; res=1.0')
#%%
# UMAP
sc.tl.umap(adata,maxiter=1000)
sc.pl.umap(adata,color=['leiden; res=1.0'],save='leiden_res1.0.png',title='UMAP')

# tSNE
sc.tl.tsne(adata)
sc.pl.tsne(adata,color=['leiden; res=1.0'],save='leiden_res1.0.png',title='tSNE')


# %%
# rank genes
# use t-test

sc.tl.rank_genes_groups(adata,method='t-test',groupby='leiden; res=1.0',key_added='t-test')
sc.pl.rank_genes_groups(adata,key='t-test', save='ranked_genes_t-test.png')

# %%

sc.tl.rank_genes_groups(adata,method='t-test',groupby='leiden; res=1.0',key_added='logreg')
sc.pl.rank_genes_groups(adata,key='logreg',save='ranked_genes_logreg.png')


# %%
'''
Step 4: Cell Types?
Now the fun part.

Using your knowledge (or the knowledge of neuroscience aficionados in your cohort, or Google),
identify some marker genes that should distinguish different brain cell types. 
You must identify at least 6 cell types. There are many resources online for 
identifying relationships between marker genes and cell types. Take a look at this 
database for an example of ways to identify cell types.

Support plots that provide evidence for your cell type assignments/what you used to 
diagnose/decide on cell types.
You can color UMAP and t-SNE plots by any gene of your choice, which is helpful for 
visualizing which clusters are enriched for which genes, and which clusters might 
correspond to a specific brain cell type.
Alternatively, you can also produce dotplots and clustermaps that allow you to see 
how a specific set of genes are associated with your clusters. Also, stacked violin plots, 
etc???

Besides these support plots, make an overall t-SNE or UMAP plot that labels your clusters 
with the cell types you think they mostly represent. Make sure to provide the support plots 
you made in order to establish your labeling. See this tutorial for an example of how to
apply labels.

marker genes to use:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6329831/

'''
# Cell types and marker genes 
# Cajal-Retzius cell markers Reln, Trp73, Lhx1, and Lhx5
# excitatory neuron clusters (5-E, 13-E, 3-E, 7-E, and 2-E) broadly expressed Bcl11b
# Lhx6, a transcription factor associated with PV and SST interneurons
# canonical markers of vasoactive intestinal peptide (Vvscode-webview://171u9smekr5kmiek5a5st0s6msqu236qchvjqdauce838r6fi1t2/5e1ce34d-b1a1-41a4-8a50-aad9aec56df5IP) cells, including Htr3a, Npas1, and Adarb2
# RG markers Hes1, Hes5, Pax6, and Ednrb
# oligodendrocytes (Cluster 16-P), marked by Olig2 expression

# https://scanpy-tutorials.readthedocs.io/en/latest/plotting/core.html
def plot_celltypes_umap(key,gene_list,adata):
    # check which genes are in adata.var
    x = [g for g in gene_list if g in list(adata.var.index)]
    sc.pl.umap(adata, color=x,save='{}.png'.format(key))

def dotplot_celltypes(cell_dict, adata):
    for key in cell_dict.keys():
        x = [g for g in cell_dict[key] if g in list(adata.var.index)]
        cell_dict[key] = x
    sc.pl.dotplot(adata, cell_dict, 'leiden; res=1.0', dendrogram=True,save='leiden_res1.0.png')

def stacked_violin_plot(cell_dict, adata):
    for key in cell_dict.keys():
        x = [g for g in cell_dict[key] if g in list(adata.var.index)]
        cell_dict[key] = x
    sc.pl.stacked_violin(adata, cell_dict, groupby='leiden; res=1.0', swap_axes=False, dendrogram=True, save='leiden_res1.0.png')

cell_dict = {
    'Cajal-Retzius':['Reln','Trp73','Lhx1','Lhx5'],
    'Excitatory_neuron':['Bcl11b'],
    'Interneurons':['Lhx6'],
    'VIP-cells':['Htr3a','Npas1','Adarb2'],
    'RG':['Hes1','Pax6','Ednrb'],
    'Oligodendrocytes':['Olig2']
}
cluster2annotation = {
     '0': 'Cajal-Retzius',
     '1': 'Excitatory_neuron',
     '2': 'Interneurons',
     '3': 'VIP-cells',
     '4': 'RG',
     '5': 'Oligodendrocytes'
}

adata.obs['cell type'] = adata.obs['leiden; res=1.0'].map(cluster2annotation).astype('category')
#%%
sc.pl.umap(adata, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10, legend_fontoutline=2,
           save='cellTypeUMAP.png')
#%%
sc.pl.tsne(adata, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10, legend_fontoutline=2,
           save='cellTypeTSNE.png')
#%%
#for key in cell_dict.keys():
#    plot_celltypes_umap(key,cell_dict[key], adata)

#%%
dotplot_celltypes(cell_dict, adata)
stacked_violin_plot(cell_dict, adata)
# %%
'''
Now we just want you to put all that together and reproduce 
the tSNE (or UMAP) with the clusters labeled based on the 
cell-types you assigned them. You have almost all the code 
set up already to do this, you just need to make the plot (-1.5 points)
'''
