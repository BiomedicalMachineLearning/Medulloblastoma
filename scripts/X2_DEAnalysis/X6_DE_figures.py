"""
(25/06/2021) -> Creates volcano plots/MDS plots to
                visualise the DE results for the data/mouse/mix, with
                emphasis on a few important genes we want to highlight.
                NOTE that this replaces the DE_figures in May_REDO/,
                but is still dependent on scripts/May_REDO/DE_figure_helpers.py

                INPUT:  data/DE_out/Pseudo_Limma*/de_results_*
                        data/DE_out/Pseudo_Limma*/data/*

                OUTPUT: figure_components/DE_figures/*
"""

################################################################################
                        # Environment setup #
################################################################################
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scripts.utils.helpers as hs
hs.setUp()

import scripts.X2_DEAnalysis.DE_figure_helpers as dhs

from sklearn.decomposition import PCA
from sklearn.preprocessing import scale

import scripts.utils.postprocessing.format.format as form

import seaborn as sb

out_plots = 'figure_components/DE_figures/'
dir_prefix = 'data/DE_out/'
data_dirs = [dir_prefix+'Pseudo_Limma_human/',
             dir_prefix+'Pseudo_Limma_mixHuman/', #NOTE, saved all genes in
                                                  # data/lcpms_all.txt for PCA
             dir_prefix+'Pseudo_Limma_mouse/']
species_ = np.array(['human', 'mix', 'mouse'])

################################################################################
                # Merging the mix results to save as one excel #
################################################################################
i, species = 1, 'mix'
mix_dir_prefix = dir_prefix+'Pseudo_Limma_mix/'
suffixes = ['Human', 'Mouse']

split_dirs = [mix_dir_prefix.strip('/')+suffix+'/' for suffix in suffixes]
dfs = []
sheet_names = []
for j, split_dir in enumerate(split_dirs):
    up_df_ = pd.read_excel(split_dir+
               f'de_results_PseudoLimma_allGenes_mix{suffixes[j]}.xlsx',
                      sheet_name=f'{species}.treatment_up', index_col=0)
    down_df_ = pd.read_excel(split_dir+
               f'de_results_PseudoLimma_allGenes_mix{suffixes[j]}.xlsx',
                    sheet_name=f'{species}.treatment_down', index_col=0)
    dfs.extend([up_df_, down_df_])
    sheet_names.extend([f'mix_{suffixes[j]}-genes_up',
                        f'mix_{suffixes[j]}-genes_down'])

form.writeDFsToExcelSheets(dir_prefix+'Table_S2_MixSpotsDE.xlsx',
                             dfs, sheet_names)

################################################################################
                    # Generating the volcano plots... #
################################################################################
sample_colors = {'A1_treated': 'orange', 'B1_treated': 'aquamarine',
                 'C1_untreated': 'hotpink', 'D1_untreated': 'orchid'}
de_colors = {'up': 'red', 'down': 'blue', 'not-de': 'k'}

for i, species in enumerate(species_):
    data_dir = data_dirs[i]
    sample_colors_ = dict(zip([sample+f'_{species}' for sample in sample_colors],
                              list(sample_colors.values())))
    if species == 'mix': # Since performed human/mouse gene DE separately.
        suffixes = ['Human', 'Mouse']
        split_dirs = [mix_dir_prefix.strip('/')+suffix+'/' for suffix in suffixes]
        up_df, down_df = None, None
        for j, split_dir in enumerate(split_dirs):
            up_df_ = pd.read_excel(split_dir+
                       f'de_results_PseudoLimma_allGenes_mix{suffixes[j]}.xlsx',
                              sheet_name=f'{species}.treatment_up', index_col=0)
            down_df_ = pd.read_excel(split_dir+
                       f'de_results_PseudoLimma_allGenes_mix{suffixes[j]}.xlsx',
                            sheet_name=f'{species}.treatment_down', index_col=0)
            if type(up_df)==type(None):
                up_df, down_df = up_df_, down_df_
            else:
                up_df = pd.concat([up_df, up_df_])
                down_df = pd.concat([down_df, down_df_])

            ### Making the diagnostic plots for human/mouse genes separately ###
            # Plotting the diagnostic plots #
            dhs.plot_diagnostics(split_dir, species+suffixes[j],
                                 sample_colors_, de_colors, out_plots)

    else: # For other species, all DE genes are in the one excel sheet.
        up_df = pd.read_excel(data_dir+
                              f'de_results_PseudoLimma_allGenes_{species}.xlsx',
                              sheet_name=f'{species}.treatment_up', index_col=0)
        down_df = pd.read_excel(data_dir+
                              f'de_results_PseudoLimma_allGenes_{species}.xlsx',
                            sheet_name=f'{species}.treatment_down', index_col=0)

    de_df = pd.concat([up_df, down_df])

    # Now plotting the volcano plot #
    adj_pvals = de_df.loc[:,'padj'].values
    log10_pvals = -np.log10( adj_pvals )
    logFC = de_df.loc[:,'log2FoldChange'].values

    sig_bool = adj_pvals < .05
    up_bool = np.logical_and(sig_bool, logFC>0)
    down_bool = np.logical_and(sig_bool, logFC<0)
    de_labels = np.array([' '*max([len(label) for label in ['up', 'down', 'non-de']])
                          ]*len(sig_bool))
    de_labels[up_bool] = 'up'
    de_labels[down_bool] = 'down'
    de_labels[sig_bool==False] = 'not-de'

    hs.plot_species_scatters(de_labels, log10_pvals, logFC, de_colors,
                             alpha=.4, linewidths=0,
                             out_path=out_plots+f'{species}_volcano_plot.pdf',
                             ymax=6, xlim=(-6,6))

    # Plotting the diagnostic plots #
    if species!='mix': #Already plotted diagnostics for mixHuman/mixMouse
        dhs.plot_diagnostics(data_dir, species,
                         sample_colors_, de_colors, out_plots)

################################################################################
                        # PCA analysis... #
################################################################################
# Now concatentating the lcpms together to create a PCA plot #
species_lcpms = []
for i, species in enumerate(species_):
    data_dir = data_dirs[i]
    lcpms_all = pd.read_csv(data_dir + 'data/lcpms_all.txt', sep='\t')
    species_lcpms.append(lcpms_all)

con_lcpms = pd.concat(species_lcpms, axis=1)
# Filtering to genes with a positive value somewhere #
gene_bool = [False]*con_lcpms.shape[0]
for i, gene in enumerate(con_lcpms.index):
    gene_bool[i] = len(np.where(con_lcpms.values[i,:]>0)[0]) > 3
print(sum(gene_bool))
human_genes = [gene for gene in con_lcpms.index.values[gene_bool]
               if 'hg38-' in gene]
mouse_genes = [gene for gene in con_lcpms.index.values[gene_bool]
               if 'mm10-' in gene]
print(len(human_genes))
print(len(mouse_genes))
# Looks good #

######## Performing the PCA #########
scaled_values = scale(con_lcpms.values.transpose(), axis=0)

pca = PCA(n_components=3)
pca_vals = pca.fit_transform(scaled_values)
print(pca.explained_variance_ratio_) # the proportion variation explained
print(sum(pca.explained_variance_ratio_)) # almost 90% explained by first two
#pca.components_ # has the matrix of feature weights #
# con_lcpms.index.values[np.argsort(pca.components_[1,:])] # ordered genes pc1
pca_loadings = pca.components_[2,:]
order = np.argsort(-pca_loadings)
plt.scatter(list(range(len(order))), pca_loadings[order])
plt.hlines(.02, 0, 10000)
plt.vlines(800, 0, .02)
plt.show()
top_genes = con_lcpms.index.values[order][0:800]
print(top_genes)

sb.clustermap(con_lcpms.loc[top_genes,:], z_score=0)
plt.show()

treat_colors = {'treated': 'red', 'untreated': 'blue'}
treat_labels = np.array([name.split('_')[1]
                         for name in con_lcpms.columns.values])
treat_label_colors = [treat_colors[label] for label in treat_labels]

species_colors = {'human': 'dodgerblue',
                  'mouse': 'mediumseagreen', 'mix': 'gold'}
species_labels = np.array([name.split('_')[-1]
                           for name in con_lcpms.columns.values])
species_label_colors = np.array([species_colors[label]
                                 for label in species_labels])

map = {'A1': 'treated', 'B1': 'treated', 'C1': 'untreated', 'D1': 'untreated'}
sample_labels = [name.split('_')[0] for name in con_lcpms.columns.values]
sample_label_colors = np.array([sample_colors[label+f'_{map[label]}']
                                for i, label in enumerate(sample_labels)])

treat_markers = {'treated': '*', 'untreated': '.'}

dhs.pca_plot(pca, pca_vals, treat_colors, treat_labels, treat_markers,
			 species_label_colors, out_plots+'pc1_pc2_pca.pdf',)

dhs.pca_plot(pca, pca_vals, treat_colors, treat_labels, treat_markers,
			 species_label_colors, out_plots+'pc1_pc3_pca.pdf', pc2=2,)






















