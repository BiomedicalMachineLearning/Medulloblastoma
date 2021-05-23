""" Created: 8th April 2021

Creating a UMAP for visualising the clusters and the merged clusters for each
of the samples.

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

import scripts.helpers as hs
hs.setUp()
import scripts.stlearn_helpers as st_hs

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

datas = [sc.read_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}.h5ad') \
         for i in range(len(samples))]

################################################################################
                    # Adding in the cluster labels #
################################################################################
meta = pandas.read_csv(data_dir+'spot_meta/clusters_human_mix.txt',
                                                          sep='\t', index_col=0)

colors = simple_pickle.loadPickle(data_dir+'spot_meta/colors.pkl')

for i, data in enumerate(datas):
    sample = samples[i]
    sample_indices = [i for i in range(meta.shape[0])
                      if sample in meta.index.values[i]]
    ordered_obs = meta.index.values[sample_indices]
    obs_names = numpy.array([f'{sample}-{name}' for name in data.obs_names])
    data.obs['cluster'] = meta.loc[obs_names,'cluster'].values
    tissue_type = meta.loc[obs_names, 'tissue_type'].astype('category')
    tissue_type.index = data.obs.index
    data.obs['tissue_type'] = tissue_type
    tissue_colors = numpy.array([colors[tissue] for tissue in
                      data.obs['tissue_type'].dtypes.categories]).astype(object)
    data.uns['tissue_type_colors'] = tissue_colors

    data.obs['sum_counts'] = data.X.sum(axis=1)
    data.obs['all_genes'] = (data.X>0).sum(axis=1)
    human_genes = [gene for gene in data.var_names if 'hg38-' in gene]
    mouse_genes = [gene for gene in data.var_names if 'mm10-' in gene]
    data.obs['human_counts'] = data[:,human_genes].X.sum(axis=1)
    data.obs['human_genes'] = (data[:, human_genes].X>0).sum(axis=1)
    data.obs['mouse_counts'] = data[:,mouse_genes].X.sum(axis=1)
    data.obs['mouse_genes'] = (data[:, mouse_genes].X > 0).sum(axis=1)

################################################################################
                    # Creating the UMAP #
################################################################################
# for i, data in enumerate(datas):
#     st.em.run_umap(data, spread=5)

for i, data in enumerate(datas):
    #st.em.run_umap(data, spread=5)
    # sc.pl.umap(data, color=['tissue_type'], show=False,
    #            legend_loc='on data')
    # vhs.dealWithPlot(True, True, True, 'figure_components/stlearn_clustering/',
    #                  f'{samples[i]}_umap_tissue_type.png', dpi=300)
    # sc.pl.umap(data, color=['cluster'], show=False)
    # vhs.dealWithPlot(True, True, True, 'figure_components/stlearn_clustering/',
    #                  f'{samples[i]}_umap_clusters.png', dpi=300)

    sc.pl.umap(data, color=['species', 'sum_counts',
                            'human_counts', 'mouse_counts'],
               show=False, size=80)
    vhs.dealWithPlot(True, True, True, 'figure_components/stlearn_clustering/',
                     f'{samples[i]}_umap_species-counts.png', dpi=300)

    sc.pl.umap(data, color=['species', 'all_genes',
                            'human_genes', 'mouse_genes'],
               show=False, size=80)
    vhs.dealWithPlot(True, True, True, 'figure_components/stlearn_clustering/',
                     f'{samples[i]}_umap_species-gene-counts.png', dpi=300)

# Saving the results #
# [data.write_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}.h5ad',compression='gzip')
#  for i, data in enumerate(datas)]


# Some gene plotting!! #

for i, data in enumerate(datas):
    sc.pl.spatial(data, color=['hg38-HRK', 'hg38-UNC5B'], show=True, size=2)


##### Spatial plots ######
for i, data in enumerate(datas):
    sc.pl.spatial(data, color=['species', 'all_genes',
                            'human_genes', 'mouse_genes'],
               show=False, size=2)
    vhs.dealWithPlot(True, True, True, 'figure_components/stlearn_clustering/',
                     f'{samples[i]}_spatial_species-gene-counts.png', dpi=300)



