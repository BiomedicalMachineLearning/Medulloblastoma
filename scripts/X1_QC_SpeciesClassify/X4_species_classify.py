"""
 Uses the species scores to classify spots by into human/mix/mouse.

         INPUT: * data/Visium8*/
                * spot_meta/*_species.txt
         OUTPUT: * data/spot_meta/species_classify_v2/
                 * figure_components/species_figure/
"""

################################################################################
                    # Environment setup #
################################################################################

import os
import stlearn as st
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib

import scripts.utils.helpers as hs
hs.setUp()
import scripts.utils.stlearn_helpers as st_hs
import scripts.X1_QC_SpeciesClassify.species_helpers as sp_hs
import scripts.utils.visualisation.helpers as vhs

out_dir = 'data/spot_meta/species_classify_v2/'
out_plots = 'figure_components/species_figure/'

################################################################################
                    # Loading data & normalising #
################################################################################
data_dir = hs.data_dir
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir+sample for sample in samples_]

filt_files = np.array( os.listdir(data_dir+'filter_ids/') )
filt_files = filt_files[np.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir+'filter_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

# Loading all the data that pass QC without subsetting #
datas = [st_hs.load_and_norm(sample_dir, spot_filters[i], None,
                             subset=False, sme=False)
         for i, sample_dir in enumerate(data_dirs)]

################################################################################
            # Getting the number of genes of each type per spot #
################################################################################
for i, data in enumerate(datas):
    # Loading in the species scores #
    treat = 'treated' if samples[i] in ['A1', 'B1'] else 'untreated'
    species_df = pd.read_csv(f'{data_dir}spot_meta/'
                             f'{samples[i]}_{treat}_species.txt', sep=' ',
                             index_col=0)
    data.obs['human_scores'] = species_df.loc[data.obs_names,'human_scores'].values
    data.obs['mouse_scores'] = species_df.loc[data.obs_names,'mouse_scores'].values

    sc.pl.spatial(data, color=['human_scores', 'mouse_scores'])

################################################################################
            ### Now creating the species classification plots ###
################################################################################
matplotlib.rc('xtick', labelsize=12)
matplotlib.rc('ytick', labelsize=12)
colors = {'human': 'dodgerblue', 'mouse': 'mediumseagreen', 'mix': 'gold'}
human_cutoffs = [[250, 400], [250, 400], [250, 300], [250, 250]]
mouse_cutoffs = [[250, 250], [400, 400], [325, 250], [560, 250]]
spot_sizes = [15, 20, 15, 15]
for i, data in enumerate(datas):
    species = np.array([' '*max([len('human'), len('mix'), len('mouse')])
                        ]*data.shape[0])
    human_counts = data.obs['human_scores'].values
    mouse_counts = data.obs['mouse_scores'].values
    human_bool = np.logical_and(human_counts>human_cutoffs[i][0],
                                mouse_counts<human_cutoffs[i][1])
    mouse_bool = np.logical_and(human_counts<mouse_cutoffs[i][0],
                                mouse_counts>mouse_cutoffs[i][1])
    mix_bool = np.logical_and(human_bool==False, mouse_bool==False)

    species[human_bool] = 'human'
    species[mouse_bool] = 'mouse'
    species[mix_bool] = 'mix'

    data.obs['species'] = pd.DataFrame(species, index=data.obs_names).iloc[:,0].astype('category')
    data.uns['species_colors'] = [colors[label] for label in
                                  data.obs['species'].cat.categories.values]

    sp_hs.plot_species_scatters(data.obs['species'].values,
                                human_counts, mouse_counts, colors
                                )
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{samples[i]}_species_gene_counts_scatter.pdf', 300)

    st.pl.cluster_plot(data, use_label='species', size=spot_sizes[i],
                       cell_alpha=.75, show_color_bar=False)
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{samples[i]}_species_spatial.pdf', 300)

################################################################################
 # For sample C1, going to spatial sub-cluster in order to isolate the distal #
                    # region & call it mouse. #
################################################################################
data = datas[2]
data_mix = data[data.obs['species'].values=='mix',:]
sc.pl.spatial(data_mix, color='species')

coords = data_mix.obsm['spatial']
data_mix.obs['coords0'] = coords[:,0]
data_mix.obs['coords1'] = coords[:,1]
sc.pl.spatial(data_mix, color=['coords0', 'coords1'])

sub_clusts = np.array(['outlier' if (coords[i,0]<450 and coords[i,1]<600) or \
                                 (coords[i,0]<450 and coords[i,1]>1630) else ''
              for i in range(coords.shape[0])])

data_mix.obs['sub_clusts'] = sub_clusts
sc.pl.spatial(data_mix, color='sub_clusts')

# Updating the species information #
mix_indices = np.where(data.obs['species'].values=='mix')[0]
outlier_indices = mix_indices[sub_clusts=='outlier']

old_species = data.obs['species'].values
old_species[outlier_indices] = 'mouse'
data.obs['species'] = old_species
sc.pl.spatial(data, color='species')

# Making plots of the new classifications #
i = 2
data = datas[i]
sp_hs.plot_species_scatters(data.obs['species'].values,
                            data.obs['human_scores'].values,
                            data.obs['mouse_scores'].values, colors
                            )
vhs.dealWithPlot(True, True, True, out_plots,
                 f'{samples[i]}_species_gene_counts_scatter.pdf', 300)

st.pl.cluster_plot(data, use_label='species', size=spot_sizes[i],
                   cell_alpha=.75, show_color_bar=False)
vhs.dealWithPlot(True, True, True, out_plots,
                 f'{samples[i]}_species_spatial.pdf', 300)

################################################################################
                # Saving the new species classifications #
################################################################################
for i, data in enumerate(datas):
    new_species_df = data.obs.loc[:, ['species', 'human_scores', 'mouse_scores']]
    new_species_df.to_csv(out_dir+f'{samples[i]}_species.txt', sep='\t')

################################################################################
            # Plotting the human/mouse genes per spot #
################################################################################
prefixes = ['A1-', 'B1-', 'C1-', 'D1-']

## Adding the human/mouse counts to the datas ##
species_gene_counts = {species_prefix:{prefix: {} for prefix in prefixes}
                       for species_prefix in ['hg38-', 'mm10-']}
species_counts = {species_prefix:{prefix: {} for prefix in prefixes}
                       for species_prefix in ['hg38-', 'mm10-']}
for i, data in enumerate( datas ):
    expr = data.X.toarray()
    bool_expr = expr > 0
    for species_prefix in ['hg38-', 'mm10-']:
        gene_bool = np.array([gene.startswith(species_prefix)
                                      for gene in data.var_names])
        gene_counts = np.apply_along_axis(hs.count_true, 1,
                                                        bool_expr[:, gene_bool])
        species_gene_counts[species_prefix][prefixes[i]] = gene_counts
        data.obsm[f'{species_prefix}gene-counts'] = pd.Series(gene_counts,
                                             index=data.obs_names).astype(float)

## Plotting the human/mouse counts per spot ##
spot_sizes = [15, 20, 15, 15]
for i, (data, prefix) in enumerate( zip(datas, prefixes) ):
    for species_prefix in ['hg38-', 'mm10-']:
        st.pl.het_plot(data, use_het=f'{species_prefix}gene-counts',
                       size=spot_sizes[i], cell_alpha=.75,
                       show_color_bar=True, vmin=0, vmax=4000)
        vhs.dealWithPlot(True, True, True, out_plots,
                        f'{prefix}{species_prefix}gene_counts_spatial.pdf', 300)












