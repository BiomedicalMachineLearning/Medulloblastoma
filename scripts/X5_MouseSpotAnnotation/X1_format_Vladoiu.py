"""
    Downloads & compiles the dev. mouse cerebellum from
        Vladiou, et al. to use a reference to annotate mouse
        spots later. Very memory intensive; ran on HPC.
        Ref. Paper: https://www.nature.com/articles/s41586-019-1158-7
        Data link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118068

        INPUT: * NONE, downloads the GEO tar from above link.

        OUTPUT: * data/third_party_data/Vladiou2019_Nature_scRNA/
                                                     Vladoiu2019_lateRef.h5ad ->
                                                             Just P7 cell types.
"""

import os, sys
import numpy, pandas
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from anndata import AnnData
import scanpy as sc

import scripts.utils.helpers as hs
hs.setUp()

import scripts.utils.preprocessing.load_data.load_10X as load_10X
import scripts.utils.visualisation.helpers as vhs

data_dir = hs.data_dir+'third_party_data/Vladiou2019_Nature_scRNA/'
raw_dir = data_dir+'GSE118068_RAW/'
data_out = data_dir
out_plots = 'figure_components/MouseAnnot_figures/'

################################################################################
                        # Reading in the data #
################################################################################
# Getting file prefixes and grouping into file-type (matrix, gene, barcode) #
prefixes = numpy.unique(['_'.join(file_name.split('_')[0:2])
                        for file_name in os.listdir(raw_dir)
                         if not file_name.startswith('.')
			and file_name.endswith('.gz')])
suffixes = ['_matrix.mtx.gz', '_genes.tsv.gz', '_barcodes.tsv.gz']

# Loading in the count matrices as dataframes #
all_counts = None
time_labels = []
for i, prefix in enumerate(prefixes):
    # Loading in the 10X data #
    paths = [raw_dir+prefix+suffix for suffix in suffixes]
    count_df = load_10X.read_10X(paths[0], paths[1], paths[2])
    n_cells = count_df.shape[1]
    time = prefix.split('_')[1] if 'E18' not in prefix else 'E18.5'
    count_df.columns = [f'{time}_{name.replace("-1", "")}' for name in count_df.columns]
    if type(all_counts) == type(None):
        all_counts = count_df
    else:
        all_counts = pandas.concat([all_counts, count_df], axis=1)
        del count_df

    print(all_counts.shape)
    # Saving the time point labels #
    time_labels += [time] * n_cells

################################################################################
        # Mapping the annotations to that provided by the authors #
################################################################################
# Formatting the cell meta data #
cell_meta = pandas.read_csv(data_dir+'cluster_annotations.csv',
                                      sep=',', index_col=0)

#subsetted to just annotated
all_counts = all_counts.loc[:, cell_meta.index]

cell_meta['time'] = [name.split('_')[0] for name in cell_meta.index]

# Converting the cluster labels to cell types, based on the paper extended
# figure provided by the authors.
cluster_to_label = {0: 'Excitatory cerebellar nuclei neurons',
		   1: 'Embryonic and postnatal GCPs-1',
		   2: 'Neural stem cells', 3: 'Unipolar brush cell and GCP progenitor',
		   4: 'Unipolar brush cells', 5: 'GABA interneurons',
		   6: 'Brainstem progenitors', 7: 'Granule cells',
		   8: 'VZ progenitors', 9: 'Unipolar brush cell precursors',
		   9: 'Unipolar brush cell precursors', 
		   10: 'Differentiating Purkinje cells',
		   11: 'Gliogenic progenitors-1', 
		   12: 'Upper rhombic lip progenitors',
		   13: 'Mesenchymal stem cells-1',
		   14: 'Purkinje cells', 15: 'Postnatal GCPs-2',
		   16: 'Post mitotic NTZ neurons', 
		   17: 'Roof plate-like stem cells',
		   18: 'Proliferating VZ progenitors',
		   19: 'Oligodendrocyte precursor cells', 
		   20: 'Gliogenic progenitors-2',
		   21: 'Astrocyte/Bergmann glia precursors',
		   22: 'Endothelial cells', 
		   23: 'Postnatal excitatory cerebellar nuclei neurons',
		   24: 'GABA interneuron precursors',
		   25: 'Pericytes', 26: 'Early proliferating VZ progenitors',
		   27: 'Mesenchymal stem cells-2', 
		   28: 'Microglia', 29: 'Meninges', 30: 'Red blood cells'}

labels = [cluster_to_label[cluster] for cluster in cell_meta['Cluster identity']]

cell_meta['cell_labels'] = labels

###### Removing any genes with litle expression ######
new_counts = pandas.DataFrame(all_counts.values,
                              index=all_counts.index, columns=all_counts.columns)

def getFreq(x):
    return len(numpy.where(x>0)[0])

gene_freqs = numpy.apply_along_axis(getFreq, 1, new_counts.values)
gene_bool = gene_freqs>50
print(sum(gene_bool))

filt_counts = new_counts.loc[gene_bool, :]

# For some reason, there are cases where a gene has a duplicate name.
# Going to filter these rare cases out, and double check downstream that
# the gene expressions match up to those expected from the paper.
index_freqs = [len(numpy.where(filt_counts.index==gene)[0])
               for gene in filt_counts.index]
genes = filt_counts.index.values[numpy.array(index_freqs)==1]

filt_counts = filt_counts.loc[genes, :]

################################################################################
                                # Normalisation #
################################################################################
cell_counts = filt_counts
cell_meta = pd.read_csv(data_out+'Vladoiu2019_cell_meta.txt',
                                                          sep='\t', index_col=0)
cell_ad = AnnData(cell_counts.transpose(), obs=cell_meta)
cell_ad.raw = cell_ad

# Normalising #
sc.pp.filter_genes(cell_ad, min_cells=10)
sc.pp.normalize_total(cell_ad)
sc.pp.log1p(cell_ad)

################################################################################
            # # Subsetting to just the mature cells # #
################################################################################
vlad_data = cell_ad

# Subsetting to just the mature cells #
# NOTE:
vlad_data_sub = vlad_data[vlad_data.obs['time'].values.astype(str) == 'P7', :]
cell_labels_sub = vlad_data_sub.obs['cell_labels'].values.astype(str)
print(np.unique(cell_labels_sub))
remove_cells = ['Mesenchymal stem cells-1',
                'Unipolar brush cell and GCP progenitor',
                'VZ progenitors']
cell_bool = [label not in remove_cells for label in cell_labels_sub]
mature_cells = np.unique(cell_labels_sub[cell_bool])

cell_labels = vlad_data.obs['cell_labels'].values.astype(str)
mature_cell_bool = [label in mature_cells for label in cell_labels]
vlad_data = vlad_data[mature_cell_bool, :]

# Just visualising the data #
sc.pp.highly_variable_genes(vlad_data,
                            min_disp=.4, min_mean=0.02, max_mean=3)
sc.pl.highly_variable_genes(vlad_data)
print(sum(vlad_data.var['highly_variable']))

sc.pp.pca(vlad_data, n_comps=100)
sc.pp.neighbors(vlad_data, n_neighbors=15)
sc.tl.umap(vlad_data)

sc.pl.umap(vlad_data, color=['cell_labels_merged', 'time', 'cell_labels', ])

######## Creating the UMAP visualisation of the data #########
fig, ax = plt.subplots(figsize=(12,6))
sc.pl.umap(vlad_data, color='cell_labels', #legend_loc='on data',
           ax=ax, show=False, size=10)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'vlad_celllabels_umap.pdf', 300)


vlad_data.write_h5ad(data_out + 'Vladoiu2019_lateRef.h5ad', compression='gzip')



