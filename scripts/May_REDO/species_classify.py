""" Create the figures for the species classification.
"""

################################################################################
                    # Environment setup #
################################################################################
import os, sys
import stlearn as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scripts.May_REDO.helpers as hs
import matplotlib
import beautifulcells.visualisation.helpers as vhs
import beautifulcells.preprocessing.QC.cell_QC as cell_QC
import scanpy as sc

out_plots = 'figure_components/species_out/'

hs.setUp()

samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
samples_2 = ['A/', 'B/', 'C/', 'D/']

################################################################################
                    # Loading in the data #
################################################################################
prefixes = ['A1-', 'B1-', 'C1-', 'D1-']

# Filtering the spots #
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']
filt_files = np.array( os.listdir(hs.data_dir+'ryan_ids/') )
filt_files = filt_files[np.argsort(filt_files)]
spot_filters = [pd.read_csv(hs.data_dir+'ryan_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]
for i in range(len(spot_filters)):
    spot_filters[i] = [prefixes[i]+bc for bc in spot_filters[i]]

datas = [hs.load(hs.data_dir+samples_[i], prefixes[i])[spot_filters[i],:]
         for i in range(len(samples_))]

# Adding in the species meta-information #
species_data = [pd.read_csv(hs.data_dir+'spot_meta/'+filt_file,
                            header=0, index_col=0, sep=' ')
                                                 for filt_file in species_metas]
for i, meta_data in enumerate(species_data):
    meta_data.index = [prefixes[i]+spot_bc for spot_bc in meta_data.index]
    datas[i].obs['species'] = meta_data.loc[:,'species'].astype('category')
    datas[i].obs['human_scores'] = meta_data.loc[:,'human_scores'].astype(float)
    datas[i].obs['mouse_scores'] = meta_data.loc[:,'mouse_scores'].astype(float)

# Storing the raw counts & normalising #
for data in datas:
    data.raw = data
    sc.pp.normalize_total(data)
    sc.pp.log1p(data)

################################################################################
                    # Making the plots #
################################################################################
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
        # Just using the original scores saved from Seurat #
        # counts = np.apply_along_axis(np.sum, 1, expr[:, gene_bool])
        # species_counts[species_prefix][prefixes[i]] = counts

colors = {'data': 'dodgerblue', 'mouse': 'mediumseagreen', 'mix': 'gold'}
spot_sizes = [15, 20, 15, 15]
for i, (data, prefix) in enumerate( zip(datas, prefixes) ):
    data.uns['species_colors'] = list(colors.values())
    data.uns['species_set'] = list(colors.keys())
    st.pl.cluster_plot(data, use_label='species', size=spot_sizes[i],
                       cell_alpha=.75, show_color_bar=False)
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{prefix}species.pdf', 300)

    for species_prefix in ['hg38-', 'mm10-']:
        st.pl.het_plot(data, use_het=f'{species_prefix}gene-counts',
                       size=spot_sizes[i], cell_alpha=.75,
                       show_color_bar=True, vmin=0, vmax=4000)
        vhs.dealWithPlot(True, True, True, out_plots,
                         f'{prefix}{species_prefix}gene_counts_spatial.pdf', 300)

    hs.plot_species_scatters(data.obs['species'].values,
                             species_gene_counts['hg38-'][prefix],
                             species_gene_counts['mm10-'][prefix], colors
                             )
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{prefix}species_gene_counts_scatter.pdf', 300)

    hs.plot_species_scatters(data.obs['species'].values,
                             data.obs['human_scores'].values,
                             data.obs['mouse_scores'].values, colors
                             )
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{prefix}species_counts_scatter.pdf', 300)















