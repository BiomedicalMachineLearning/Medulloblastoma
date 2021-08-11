"""
Generating gene plot panels to show GFAP
                     expression as prognostic for astrocytes.

           INPUT: data/scanpy_h5ad/*_all_species_SME.h5ad
                  data/spot_meta/species_classify_v2/
                                   *Vladoiu_singleR_scores_mouse.txt
           OUTPUT: figure_components/MouseAnnot_figures/
                                                   *Gfap_spatial.pdf
"""

################################################################################
                        # Environment setup #
################################################################################

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import scripts.utils.helpers as hs
hs.setUp()
import scripts.utils.stlearn_helpers as st_hs
import scripts.utils.visualisation.helpers as vhs

data_dir = hs.data_dir
annot_dir = 'data/spot_meta/species_classify_v2/'
out_plots = 'figure_components/MouseAnnot_figures/'

################################################################################
                    # Loading data & normalising #
################################################################################
samples = ['A1', 'B1', 'C1', 'D1']
datas = [sc.read_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}_all_species_SME.h5ad')
         for i in range(len(samples))]

################################################################################
                    # Adding the mouse annotations #
################################################################################
datas_mouse = []
for i, data in enumerate(datas):
    datas_mouse.append( st_hs.species_split(data, species='mouse') )

# Adding in the SingleR annotations #
annot_files = [f'{samp}_Vladoiu_singleR_scores_mouse.txt' for samp in samples]
annot_dfs = [pd.read_csv(annot_dir+annot_file, index_col=0, sep='\t')
                                                  for annot_file in annot_files]
annot_all = pd.concat(annot_dfs)

# Making a copy to include the clustering with the cell type annotations #
#### Finding good colour schemes... #######
labels = annot_all.loc[:,'labels'].values.astype(str)
label_set = np.unique(labels)
colors = vhs.getColors(labels, label_set)
colors['Astrocyte/Bergmann glia precursors'] = 'springgreen'
colors['Endothelial cells'] = 'magenta'
colors['Pericytes'] = 'gold'
endo = colors['Meninges']
colors['Meninges'] = colors['Granule cells']
colors['Granule cells'] = endo

for i, data in enumerate(datas_mouse):
    annot_df = annot_dfs[i]
    annot_df.index = [bc.split('_')[-1] for bc in annot_df.index]

    data.obs['Vladiou'] = annot_df.loc[data.obs_names, 'labels']
    data.obs['Vladiou'] = data.obs['Vladiou'].astype('category')
    label_set_i = data.obs['Vladiou'].cat.categories
    colors_i = [colors[label] for label in label_set_i]
    data.uns['Vladiou_colors'] = colors_i

    astro_bool = data.obs['Vladiou']=='Astrocyte/Bergmann glia precursors'
    top, bottom = max(data.obsm['spatial'][:,0]), max(data.obsm['spatial'][:,1])
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.spatial(data[astro_bool,:], color='Vladiou', ax=ax,
                  )

################################################################################
                # Making the plots of the gene expression #
################################################################################
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("",
                                                 ["mediumpurple","springgreen"])
for i, data in enumerate(datas):
    fig, ax = plt.subplots(figsize=(12,8))
    sc.pl.spatial(data, color='mm10-Gfap', cmap=cmap, vmax=1.5, show=False,
                  ax=ax, frameon=False)
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{samples[i]}_Gfap_spatial.pdf', 300)














