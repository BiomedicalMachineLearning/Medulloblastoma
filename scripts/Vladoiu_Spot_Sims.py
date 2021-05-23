""" The purpose of this script:
Measures the cosine similarity between the human/mix spots &
the pseudobulked cell types !

    OUTPUT: data/spot_meta/vladoiu_spot_sims.txt
"""

################################################################################
                    # Environment setup #
################################################################################

import os
import stlearn as st
import scanpy as sc
import pandas as pd
import pandas
import numpy
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.spatial.distance import cosine
from scipy.stats import zscore

from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

import scripts.helpers as hs
hs.setUp()
import scripts.stlearn_helpers as st_hs
import beautifulcells.clustering.marker_labelling.marker_labelling \
                                                             as marker_labelling
import beautifulcells.normalisation.Discretize.Discretize as discretize

import beautifulcells.visualisation.heatmap.heatmap_helpers as heat_hs
import beautifulcells.visualisation.helpers as vhs

import beautifulcells.preprocessing.load_data.simple_pickle as simple_pickle

################################################################################
                    # Loading data & normalising #
################################################################################
# read in visium dataset downloaded from: support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Breast_Cancer_Block_A_Section_2
data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir+sample for sample in samples_]
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']

filt_files = numpy.array( os.listdir(data_dir+'ryan_ids/') )
filt_files = filt_files[numpy.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir+'ryan_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

datas = [st_hs.load_and_norm(sample_dir, spot_filters[i], species_metas[i],
                             sme=False, scale=True, subset=True)
         for i, sample_dir in enumerate(data_dirs)]

##### Concatenating into one dataframe ########
# Getting overlapping genes #
genes = set(datas[0].var_names)
data_dfs = [datas[0].to_df().transpose()]
for i in range(1, len(datas)):
    genes = genes.intersection(set(datas[i].var_names))
    data_dfs.append( datas[i].to_df().transpose() )
genes = np.array(list(genes))

human_genes = [gene for gene in genes if 'hg38-' in gene]
human_dfs = [data_df.loc[human_genes, :] for data_df in data_dfs]

human_df = pd.concat(human_dfs, axis=1)
human_genes_formatted = [gene.replace('hg38-', '') for gene in human_df.index]
human_df.index = human_genes_formatted

###### Loading in the Vladoiu_PseudoBulk counts to compare against #############
ref = pandas.read_csv('data/third_party_data/Vladoiu2019_Nature_scRNA/'
                      'Vladoiu2019_pseudo_counts.txt', sep='\t', index_col=0)
# Trying to subset to more reliable genes, just top 50 from each cell type #
de_genes = pandas.read_csv('data/third_party_data/Vladoiu2019_Nature_scRNA/'
                                  'vladoiu_de_genes.txt', sep='\t', index_col=0)
top_de = []
n_top = 50
de_results = {}
for cell_type in de_genes.columns:
    present_genes = [gene for gene in de_genes.loc[:,cell_type]
                     if gene.upper() in human_df.index]
    de_results[cell_type] = [gene.upper() for gene in present_genes]
    if len(present_genes) < n_top:
        print(cell_type, len(present_genes))
    top_de.extend( present_genes[0:n_top] )
top_de = np.unique(top_de)

ref = ref.loc[top_de, :]

# Converting gene names to human format #
ref_genes = [gene.upper() for gene in ref.index]
ref.index = ref_genes

overlap = [gene for gene in human_df.index
           if gene in ref_genes and gene!='PISD']

ref = ref.loc[overlap, :]
human_df = human_df.loc[overlap, :]

# normalising #
ref_norm = st_hs.norm_df(ref, scale=True)

# Getting cosine similarity between each tumor cluster & mouse cell types
sims = numpy.zeros((human_df.shape[1], ref_norm.shape[1]))
for i in range(sims.shape[0]):
    tumor_expr = human_df.values[:, i]
    for j in range(sims.shape[1]):
        cell_expr = ref_norm.values[:, j]
        sims[i, j] = 1-cosine(tumor_expr, cell_expr)
sims = pandas.DataFrame(sims, index=human_df.columns,
                        columns=ref.columns)
sims.to_csv(data_dir+'spot_meta/vladoiu_spot_sims.txt', sep='\t')











