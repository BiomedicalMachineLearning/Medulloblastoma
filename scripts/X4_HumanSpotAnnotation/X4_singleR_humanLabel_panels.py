"""
    Label the data according to the human labels from the
                                        Fetal Human 3 reference & output figure.

      INPUT: * data/scanpy_h5ads/*_all_species_SME.h5ad
             * data/spot_meta/species_classify_v2/
                                           *FetalBrain3_singleR_scores_human.txt

      OUTPUT: * figure_components/HumanAnnot_figures/
                                  *FetalBrain3Labels_spatial.pdf
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
out_plots = 'figure_components/HumanAnnot_figures/'

################################################################################
                    # Loading data & normalising #
################################################################################
samples = ['A1', 'B1', 'C1', 'D1']
datas = [sc.read_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}_all_species_SME.h5ad')
         for i in range(len(samples))]

################################################################################
                    # Visualising for the mouse #
################################################################################
datas_human = []
for i, data in enumerate(datas):
    datas_human.append( st_hs.species_split(data, species='human') )

# Adding in the SingleR annotations #
annot_files = [f'{samp}_FetalBrain3_singleR_scores_human.txt'
                                                            for samp in samples]
annot_dfs = [pd.read_csv(annot_dir+annot_file, index_col=0, sep='\t')
                                                  for annot_file in annot_files]

# Making a copy to include the clustering with the cell type annotations #
all_labels = []
for i, data in enumerate(datas_human):
    annot_df = annot_dfs[i]
    annot_df.index = [bc.split('_')[-1] for bc in annot_df.index]

    data.obs['FetalBrain3'] = annot_df.loc[data.obs_names, 'labels']
    all_labels.extend(np.unique(data.obs['FetalBrain3'].values.astype(str)))

    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.spatial(data, color='FetalBrain3', ax=ax)

#### Finding good colour schemes... #######
labels = all_labels
label_set = np.unique(labels)
colors = vhs.getColors(labels, label_set)
colors['Proliferating cell_KIAA0101_high'] = 'aqua'
colors['Proliferating cell_UBE2C_high'] = 'navy'
colors['Astrocyte'] = 'yellow'
colors['Unknown'] = 'grey'
colors['Radial glia_HES1_high'] = 'salmon'

for i, data in enumerate(datas_human):
    label_set_i = data.obs['FetalBrain3'].cat.categories
    colors_i = [colors[label] for label in label_set_i]
    data.uns['FetalBrain3_colors'] = colors_i

    fig, ax = plt.subplots(figsize=(12,8))
    st.pl.cluster_plot(data, use_label='FetalBrain3', ax=ax, size=20,
                       crop=False, show_color_bar=False)
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{samples[i]}_FetalBrain3Labels_spatial.pdf', 300)

""" NOTE that in the paper:
        * Proliferating cell_KIAA0101_high = Proliferating cell_1
        * Proliferating cell_UBE2C_high = Proliferating cell_2
        * Neuron_NEUROD6 high = Neuron
        * Radial glia_HES1_high = Radial glia
        
 These changes were to prevent the reader from over-interpreting the results;
    i.e. that we were saying spots labelled Proliferating cell_KIAA0101_high were
        high in KIAA0101. When this analysis can only tell us if more similiar to 
        proliferating cell versus the other cell types. 
"""



































