"""
(08/07/2021) -> Based on the analysis conducted in mix_analysis.py,
                loads in the data & generates the figure panels to explain
                the observations.

                INPUT: data/scanpy_h5ads/integrated_border.h5ad
                       data/thirdy_party_data/Hovestadt2019_Nature_scRNA/
                                                          hov_scrna.h5ad
                        data/thirdy_party_data/Vladoiu2019_Nature_scRNA/
                                                Vladoiu2019_lateRef.h5ad

                OUTPUT: figure_components/species_classify_v2/
                                                       mix_only/panels/*
"""

################################################################################
                                # Environment setup #
################################################################################

import os
import stlearn as st
import scanpy as sc
import pandas as pd
import numpy
import numpy as np
import matplotlib.pyplot as plt

from anndata import AnnData
import scripts.helpers as hs

hs.setUp()
import scripts.stlearn_helpers as st_hs
import scripts.species_classify_v2.label_transfer_helpers as lt_hs
import scripts.species_classify_v2.mix_only.mix_helpers as mix_hs
import beautifulcells.visualisation.heatmap.heatmap_helpers as heat_hs
import beautifulcells.visualisation.helpers as vhs

data_dir = 'data/scanpy_h5ads/'
out_plots = 'figure_components/species_classify_v2/mix_only/panels/'
out_plots2 = 'figure_components/species_classify_v2/mix_only/supp_panels/'
hov_dir = 'data/third_party_data/Hovestadt2019_Nature_scRNA/'
vlad_dir = 'data/third_party_data/Vladoiu2019_Nature_scRNA/'

################################################################################
                        # Loading the Data #
################################################################################
samples = ['A1', 'B1', 'C1', 'D1']
vis_data = sc.read_h5ad(data_dir+'integrated_border.h5ad')
hov_data = sc.read_h5ad(hov_dir+'hov_scrna.h5ad')
hov_data = hov_data[hov_data.obs['leiden']!='6',:]
vlad_data = sc.read_h5ad(vlad_dir+'Vladoiu2019_lateRef.h5ad')

# Adding in treatment labels to the vis_data #
treat_samps = {'treat': ['A1', 'B1'], 'control': ['C1', 'D1']}
treat_labels = np.array([' '*max([len(l_) for l_ in treat_samps])]*len(vis_data))
for sample in samples:
    sample_bool = vis_data.obs['sample'].values == sample
    for treat in treat_samps:
        if sample in treat_samps[treat]:
            treat_labels[sample_bool] = treat
vis_data.obs['treat'] = treat_labels

################################################################################
   # First generating plots to highlight the border cells in visium data #
################################################################################
vis_data.obs['border'] = ['border']*len(vis_data)
vis_data.uns['border_colors'] = np.array(['red'])
# spatial_coords = vis_data.obsm['spatial']
# min_x, max_x, min_y, max_y = min(spatial_coords[:,0]), max(spatial_coords[:,0]), \
#                              min(spatial_coords[:,1]), max(spatial_coords[:,0])
for i in range(len(samples)):
    sc.pl.spatial(vis_data[vis_data.obs['sample']==samples[i],:], color='border',
                  library_id=f'Visium8_{samples[i]}_Hybrid',
                  crop_coord=(100, 1800,
                              100, 1900),
                  show=False
                  )
    vhs.dealWithPlot(True, True, True,
                     out_plots, f'{samples[i]}_border_spatial.pdf', 300)

################################################################################
                # The integrated UMAP with clustering/samples #
################################################################################
vis_data.uns['leiden_colors'] = np.array(['dodgerblue','aquamarine','greenyellow'])
vis_data.uns['sample_colors'] = np.array(['orange','aquamarine',
                                          'hotpink', 'orchid'])
sc.pl.umap(vis_data, color='leiden', show=False)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'border_leiden_umap.pdf', 300)
sc.pl.umap(vis_data, color='sample', show=False)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'border_sample_umap.pdf', 300)

# Showing the leiden clusters spatially #
for i in range(len(samples)):
    sc.pl.spatial(vis_data[vis_data.obs['sample']==samples[i],:], color='leiden',
                  library_id=f'Visium8_{samples[i]}_Hybrid',
                  crop_coord=(100, 1800,
                              100, 1900),
                  show=False
                  )
    vhs.dealWithPlot(True, True, True,
                     out_plots, f'{samples[i]}_leiden_spatial.pdf', 300)

# Showing the UMAPs & spatials of the cell type localisations #
cell_types = ['0', '1', '2', '3', '4', '5', 'GABA interneurons',
       'Glutamatergic neurons', 'Granule cells', 'Microglia', 'Purkinje cells',
       'Unipolar brush cells', 'glial', 'vasculature']
cmap='winter'
for i in range(len(samples)):
    for j in range(len(cell_types)):
        sc.pl.spatial(vis_data[vis_data.obs['sample']==samples[i],:],
                        library_id=f'Visium8_{samples[i]}_Hybrid',
                      color=cell_types[j], cmap=cmap, show=False,
                      crop_coord=(100, 1800,
                                  100, 1900),
                      )
        vhs.dealWithPlot(True, True, True,
                         out_plots, f'{cell_types[j]}_border_spatial.pdf', 300)

# The umaps #
for j in range(len(cell_types)):
    sc.pl.umap(vis_data,
                  color=cell_types[j], cmap=cmap, show=False,
                  )
    vhs.dealWithPlot(True, True, True,
                     out_plots, f'{cell_types[j]}_border_umap.pdf', 300)

################################################################################
                        # The Hovestadt UMAPs #
################################################################################
sc.pl.umap(hov_data, color='leiden', show=False)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'hov_celltype_umap.pdf', 300)

sc.pl.umap(hov_data, color='exper_id', show=False)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'hov_experid_umap.pdf', 300)

sc.pl.umap(hov_data, color='exper_type', show=False)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'hov_expertype_umap.pdf', 300)

######### UMAPs of gene expression to show labelling ##########
markers = {'G2M': ['CDC20', 'CENPA'], 'S': ['H2AFX', 'RRM2'],
           'G1': ['MCM6', 'UNG'], # markers from Zhang (2019) Cancer Cell
            # from third_party_data/Hovestadt2019_Nature_scRNA/
            # marker_genes_logodds_ranked.xlsx
           '0': ['GFAP', 'FABP7', 'SOX9'], '1': ['EBF3', 'STMN4'],
           '2': ['TPM2', 'CCS', 'PRDX4'],
           '3': ['NHLH1'], '4': ['TOP2A', 'TPX2'], '5': ['EBF2', 'POU4F1'],
           '6': ['CDKN2C'] # Suggests G0, not actively dividing.
                           #https://r.search.yahoo.com/_ylt=AwrgDd3tnudgZhYAF2ML5gt.;_ylu=Y29sbwNncTEEcG9zAzEEdnRpZANEMTAyNl8xBHNlYwNzcg--/RV=2/RE=1625821038/RO=10/RU=https%3a%2f%2fpubmed.ncbi.nlm.nih.gov%2f31851939%2f/RK=2/RS=ohqcQgrjU.vAf7JxrKwwlvKx9M4-
           }
genes = []
[genes.extend(markers[cell_type]) for cell_type in markers]

for gene in genes:
    sc.pl.umap(hov_data, color=gene, show=False)
    vhs.dealWithPlot(True, True, True,
                     out_plots2, f'hov_{gene}_umap.pdf', 300)

################################################################################
                        # The Vladoiu UMAPs #
################################################################################
sc.pl.umap(vlad_data, color='cell_labels_merged', show=False)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'vlad_celllabelsmerged_umap.pdf', 300)

fig, ax = plt.subplots(figsize=(12,6))
sc.pl.umap(vlad_data, color='cell_labels', #legend_loc='on data',
           ax=ax, show=False, size=10)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'vlad_celllabels_umap.pdf', 300)

fig, ax = plt.subplots(figsize=(7,6))
sc.pl.umap(vlad_data, color='time', #legend_loc='on data',
           ax=ax, show=False, size=10)
vhs.dealWithPlot(True, True, True,
                 out_plots, f'vlad_time_umap.pdf', 300)

################################################################################
                # Creating cell colocalisation heatmaps #
################################################################################
# Creating cell colocalisation heatmaps #
from sklearn.metrics.pairwise import pairwise_distances
from scipy.stats import spearmanr
import seaborn as sb

def spearman(vals1, vals2):
    rho, p = spearmanr(vals1, vals2)
    return rho

hov_set = np.unique(hov_data.obs['leiden'].values.astype(str))
vlad_set = np.unique(vlad_data.obs['cell_labels_merged'].values.astype(str))

metric = spearman

# got the row/col orders from the clustermap; out.data2d #
row_order = ['0', '5', '1', '3', '2', '4']
col_order = ['Glutamatergic neurons', 'Purkinje cells', 'GABA interneurons',
             'Unipolar brush cells', 'glial', 'Granule cells', 'Microglia',
             'vasculature']

treats = list(treat_samps.keys())
for i in range(len(treats)):
    vis_sub = vis_data[vis_data.obs['treat']==treats[i],:]
    human_scores = vis_sub.obs.loc[:, hov_set].values.transpose()
    mouse_scores = vis_sub.obs.loc[:, vlad_set].values.transpose()

    distances = pairwise_distances(human_scores, mouse_scores, metric=metric)
    dists_comp = pd.DataFrame(distances, index=hov_set, columns=vlad_set)
    dists_comp = dists_comp.loc[row_order, col_order]
    out = sb.clustermap(dists_comp, row_cluster=False, col_cluster=False,
                  cmap='PiYG', vmin=-0.4, vmax=0.4,
                  )
    vhs.dealWithPlot(True, True, True,
                     out_plots, f'{treats[i]}_species_comp_heatmap.pdf', 300)

    # distances = pairwise_distances(human_scores, human_scores, metric=metric)
    # dists_human = pd.DataFrame(distances, index=hov_set, columns=hov_set)
    # sb.clustermap(dists_human, row_cluster=True, col_cluster=True,
    #               cmap='PiYG',
    #               )
    # vhs.dealWithPlot(True, True, True,
    #                  out_plots, 'human_comp_heatmap.pdf', 300)


# Compile this information into a figure.



############ Cross-referencing data from other figures... ########

sc.pl.violin(vis_data, ['hg38-SLC1A2'], groupby='treat', show=False)
vhs.dealWithPlot(True, True, True, out_plots,
                 'border_SLC1A2_treat_violins.pdf', 300)

for samp in samples:
    sc.pl.spatial(vis_data[vis_data.obs['sample']==samp,:],
                  color=['hg38-SLC1A2', 'hg38-TYMS'],
                  library_id=f'Visium8_{samp}_Hybrid', show=False)
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{samp}_SLC1A2_TYMS_spatial.pdf',300)


























