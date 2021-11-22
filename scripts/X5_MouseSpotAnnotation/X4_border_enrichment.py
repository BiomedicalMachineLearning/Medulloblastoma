"""
Perform FET to test for enrichment of
    dominant spot cell types for the border spots to provide
    statistical evidence for astrocytes at the border.

    INPUT: * data/Visium8*/
           * data/scanpy_h5ad/*_all_species_SME.h5ad
           * data/spot_meta/species_classify_v2/
                           *Vladoiu_singleR_scores_mouse.txt

    OUTPUT: * data/spot_meta/astro_enrich_stats.xlsx
            * figure_components/MouseAnnot_figures/
                                         *border_spatial.pdf
"""


################################################################################
                        # Environment setup #
################################################################################

import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

import seaborn as sb

import scripts.utils.helpers as hs
hs.setUp()
import scripts.utils.stlearn_helpers as st_hs
import scripts.utils.visualisation.helpers as vhs
import scripts.utils.postprocessing.format.format as form

data_dir = hs.data_dir
out_dir = 'data/spot_meta/'
out_dir2 = 'data/scanpy_h5ads/'
annot_dir = 'data/spot_meta/species_classify_v2/'
out_plots = 'figure_components/MouseAnnot_figures/'

################################################################################
                    # Loading data & normalising #
################################################################################
samples = ['A1', 'B1', 'C1', 'D1']
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir + sample for sample in samples_]
species_metas = ['species_classify_v2/A1_species.txt',
                 'species_classify_v2/B1_species.txt',
                 'species_classify_v2/C1_species.txt',
                 'species_classify_v2/D1_species.txt']

filt_files = np.array(os.listdir(data_dir + 'filter_ids/'))
filt_files = filt_files[np.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir + 'filter_ids/' + filt_file,
                          header=None).values[:, 0] for filt_file in filt_files]

datas = [st_hs.load_and_norm(sample_dir, spot_filters[i], species_metas[i],
                             subset=True, sme=False,
                             subset_labels='mix')
         for i, sample_dir in enumerate(data_dirs)]

################################################################################
                    # Subsetting to just the border spots #
################################################################################
###### Subsetting to just border spots ######
row_coords = [[[20, 60]], [[0, 100], [62, 100], [30, 50], [0, 15]], #treat
              [[0, 45], [40,77], [34,63], [20,30], [55,100]], #C1
              [[0,45],[27,45], [15,45], [38, 58], [58,75], [65,75]]] #D1
col_coords = [[[0, 50]], [[0, 35], [0, 200], [0, 50], [70, 200]], #treat
              [[30, 112], [20,66], [2,40], [110,200], [95,150]], #C1
              [[80,200], [30,200], [66,200], [26,76], [26,62], [10,26]]] #D1
for i in range(len(datas)):
    sc.pl.spatial(datas[i], color=['array_row', 'array_col'])
    keep_bool = np.array([True] * len(datas[i]))
    for j in range(len(row_coords[i])):
        row_bool = np.logical_and(
            datas[i].obs['array_row'].values > row_coords[i][j][0],
            datas[i].obs['array_row'].values < row_coords[i][j][1])
        col_bool = np.logical_and(
            datas[i].obs['array_col'].values > col_coords[i][j][0],
            datas[i].obs['array_col'].values < col_coords[i][j][1])
        outlier_bool = np.logical_and(row_bool, col_bool)
        keep_bool = np.logical_and(keep_bool, outlier_bool == False)

        sc.pl.spatial(datas[i][keep_bool, :], color=['array_row', 'array_col'])

    datas[i] = datas[i][keep_bool,:]

# Concatenating the border spots #
border_data = datas[0].concatenate(datas[1:],
                                   batch_key="sample", uns_merge="first",
                                    batch_categories=samples
                                   )

# Saving for use downstream #
border_data.write_h5ad(out_dir2+'integrated_border.h5ad', compression='gzip')

for samp in samples:
    sc.pl.spatial(border_data[border_data.obs['sample']==samp,:], color='sample',
                  library_id=f'Visium8_{samp}_Hybrid')

# Getting list of barcodes in each sample corresponding to the border spots #
border_bcs = {}
for samp in samples:
    border_bcs[samp] = [bc.split(f'-{samp}')[0] for bc in border_data.obs_names
                        if samp in bc]
border_totals = [len(border_bcs[samp]) for samp in border_bcs]

################################################################################
   # First generating plots to highlight the border cells in visium data #
################################################################################
border_data.obs['border'] = ['border']*len(border_data)
border_data.uns['border_colors'] = np.array(['red'])
for i in range(len(samples)):
    sc.pl.spatial(border_data[border_data.obs['sample']==samples[i],:],
                  color='border', library_id=f'Visium8_{samples[i]}_Hybrid',
                  crop_coord=(100, 1800, 100, 1900), show=False
                  )
    vhs.dealWithPlot(True, True, True,
                     out_plots, f'{samples[i]}_border_spatial.pdf', 300)

################################################################################
                    # Adding the mouse annotations #
################################################################################
datas = [sc.read_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}_all_species_SME.h5ad')
                                                    for i in range(len(samples))]
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

    top, bottom = max(data.obsm['spatial'][:,0]), max(data.obsm['spatial'][:,1])
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.spatial(data, color='Vladiou', ax=ax,
                  )

################################################################################
          # Performing FET cell type enrichment for border spots #
################################################################################
samp_totals = [annot_df.shape[0] for annot_df in annot_dfs]
astro_bcs = {}
for i, samp in enumerate(samples):
    bc_labels = annot_dfs[i].loc[:,'labels']
    astro_bcs[samp] = bc_labels[bc_labels=='Astrocyte/Bergmann glia precursors']\
                                                       .index.values.astype(str)

cont_tables = {}
cont_dfs = []
for i, samp in enumerate(samples):
    cont_table = np.zeros((2,2))
    astros_border = sum([bc in astro_bcs[samp] for bc in border_bcs[samp]])
    astros_notBorder = sum([bc not in border_bcs[samp]
                            for bc in astro_bcs[samp]])
    notAstro_border = sum([bc not in astro_bcs[samp]
                           for bc in border_bcs[samp]])
    notAstro_notBorder = samp_totals[i] - \
                                (astros_border+astros_notBorder+notAstro_border)
    cont_table[0,:] = [astros_border, astros_notBorder]
    cont_table[1,:] = [notAstro_border, notAstro_notBorder]

    cont_tables[samp] = cont_table
    cont_dfs.append( pd.DataFrame(cont_table,
                                 index=['Astrocyte spot', 'Not Astrocyte spot'],
                                 columns=['Border spot', 'Not Border spot']) )

"""
{'A1': array([[  78.,  205.],
              [ 182., 1473.]]),
 'B1': array([[  53.,  119.],
              [ 101., 1337.]]),
 'C1': array([[  90.,  208.],
              [ 133., 1007.]]),
 'D1': array([[ 141.,  245.],
              [ 160., 1495.]])}
"""

# Performing FET #
stats = {}
for i, samp in enumerate(samples):
    stats[samp] = fisher_exact(cont_tables[samp], alternative='greater')

"""
{'A1': (3.0794425087108013, 3.4787115483051444e-12), #oddsratio, p-value
 'B1': (5.895748398369249, 2.242040950590036e-17),
 'C1': (3.276098901098901, 1.44149353016164e-13),
 'D1': (5.377423469387755, 2.832076127957724e-34)}
"""

samp_names = ['Palbo A', 'Palbo B', 'Control C', 'Control D']
stats_df = pd.DataFrame(index=['oddsratio', 'p-value'],
                        columns=samp_names)
for i, samp in enumerate(stats):
    stats_df.iloc[:,i] = stats[samp]

################################################################################
                     # Saving as excel sheet #
################################################################################
sheet_names = ['Border Astro FET-Stats']+\
              [samp+' continency table' for samp in samp_names]
dfs = [stats_df]+cont_dfs

form.writeDFsToExcelSheets(out_dir+'astro_enrich_stats.xlsx',
                             dfs, sheet_names)



