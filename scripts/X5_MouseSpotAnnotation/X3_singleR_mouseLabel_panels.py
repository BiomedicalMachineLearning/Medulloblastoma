"""
Label the data according to the mouse labels from the Vladiou & output figure.

                          INPUT: data/scanpy_h5ad/*_all_species_SME.h5ad
                                 data/spot_meta/species_classify_v2/
                                               *Vladoiu_singleR_scores_mouse.txt
                          OUTPUT: figure_components/MouseAnnot_figures/
                                                         *VladLabels_spatial.pdf
"""

################################################################################
                        # Environment setup #
################################################################################

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import stlearn as st

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
# read in visium dataset downloaded from: support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Breast_Cancer_Block_A_Section_2
samples = ['A1', 'B1', 'C1', 'D1']
datas = [sc.read_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}_all_species_SME.h5ad')
         for i in range(len(samples))]

################################################################################
                    # Visualising for the mouse #
################################################################################
datas_mouse = []
for i, data in enumerate(datas):
    datas_mouse.append( st_hs.species_split(data, species='mouse') )

# Adding in the SingleR annotations #
#annot_files = [file_ for file_ in os.listdir(annot_dir) if 'STABH7' in file_]
#annot_files = np.array(annot_files)[np.argsort(annot_files)]
annot_files = [f'{samp}_Vladoiu_singleR_scores_mouse.txt' for samp in samples]
annot_dfs = [pd.read_csv(annot_dir+annot_file, index_col=0, sep='\t')
                                                  for annot_file in annot_files]

# Making a copy to include the clustering with the cell type annotations #
all_labels = []
for i, data in enumerate(datas_mouse):
    annot_df = annot_dfs[i]
    annot_df.index = [bc.split('_')[-1] for bc in annot_df.index]

    data.obs['Vladiou'] = annot_df.loc[data.obs_names, 'labels']
    all_labels.extend(np.unique(data.obs['Vladiou'].values.astype(str)))

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.spatial(data, color='Vladiou', ax=ax)

#### Finding good colour schemes... #######
labels = all_labels
label_set = np.unique(labels)
colors = vhs.getColors(labels, label_set)
colors['Astrocyte/Bergmann glia precursors'] = 'springgreen'
#colors['Microglia'] = 'brown'
colors['Endothelial cells'] = 'magenta'
colors['Pericytes'] = 'gold'
endo = colors['Meninges']
colors['Meninges'] = colors['Granule cells']
colors['Granule cells'] = endo

for i, data in enumerate(datas_mouse):
    label_set_i = data.obs['Vladiou'].cat.categories
    colors_i = [colors[label] for label in label_set_i]
    data.uns['Vladiou_colors'] = colors_i

    fig, ax = plt.subplots(figsize=(12,8))
    #sc.pl.spatial(data, color='Vladiou', ax=ax, show=False)
    st.pl.cluster_plot(data, use_label='Vladiou', ax=ax, size=20, crop=False)
    vhs.dealWithPlot(True, True, True,
                     out_plots, f'{samples[i]}_VladLabels_spatial.pdf', 300)






































