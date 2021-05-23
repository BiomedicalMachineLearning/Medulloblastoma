""" Loading in the data from the Pseudo_Limma approach for visualisation.

    OUTPUT:
        figure_components/PseudoLimma/DE_figures/*
"""

################################################################################
                    # Environment setup #
################################################################################
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scripts.May_REDO.helpers as hs
import scripts.May_REDO.DE_figure_helpers as dhs
from sklearn.preprocessing import scale

out_plots = 'figure_components/PseudoLimma/DE_figures/'

hs.setUp()

################################################################################
                    # Loading in the data #
################################################################################
data_dir = 'data/DE_out/Pseudo_TMM_Limma_Voom/data/'
lcpms = pd.read_csv(data_dir+'lcpms.txt', sep='\t')
rles = pd.read_csv(data_dir+'rles.txt', sep='\t')
gene_stats = pd.read_csv(data_dir+'gene_stats.txt', sep='\t')

colors = {'A1_treated_human': 'firebrick', 'B1_treated_human': 'royalblue',
          'C1_untreated_human': 'limegreen', 'D1_untreated_human': 'orchid'}

# The density plots #
dhs.density_plot(lcpms, colors, out_plots+'norm_density.pdf')
dhs.density_plot(rles, colors, out_plots+'rle_density.pdf',
                 vmin=-1.5, vmax=1.5)
dhs.density_plot(gene_stats.iloc[:,0:1],
                 {gene_stats.columns[0]: 'blue'},
                 out_plots+'logFC_density.pdf', vertical=True,
                 )

# The scatter plot #
de_colors = {'up': 'red', 'not-de': 'k', 'down': 'blue'}
avg_expr = lcpms.values.mean(axis=1)
hs.plot_species_scatters(gene_stats.loc[:,'de_status'].values,
                         gene_stats.loc[:,'logFC'].values, avg_expr,
                         de_colors, alpha=.4, s=10, linewidths=0,
                         out_path=out_plots+'logFC_avgExpr_scatter.pdf')

# The pie chart #
# Pie chart, where the slices will be ordered and plotted counter-clockwise:
dhs.pie_plot(gene_stats, out_plots+'DE_pie.pdf')

# Heatmap of the DE genes #
dhs.plot_heatmap(gene_stats, lcpms, out_plots+'DE_heatmap.pdf',
                 alpha=.7, cbar_pos=None, xticklabels=False)










""" Junk Code:

###### Making Venn diagrams of overlap with stable genes #########
stable_genes = pd.read_csv('data/third_party_data/stable_genes.txt', sep='\t')
human_stable = {'hg38-'+gene for gene in stable_genes.loc[:,'human'].values}
mouse_stable = {'mm10-'+gene for gene in stable_genes.loc[:,'mouse'].values
                         if type(gene)!=float}
                         
# Venn of overlap with stably expressed genes #
human_up = {gene for i, gene in enumerate(gene_stats.index)
            if gene.startswith('hg38-') and labels[i]=='up'}
human_down = {gene for i, gene in enumerate(gene_stats.index)
            if gene.startswith('hg38-') and labels[i]=='down'}
mouse_up = {gene for i, gene in enumerate(gene_stats.index)
            if gene.startswith('mm10-') and labels[i]=='up'}
mouse_down = {gene for i, gene in enumerate(gene_stats.index)
            if gene.startswith('mm10-') and labels[i]=='down'}

from matplotlib_venn import venn3

venn3([mouse_stable, mouse_up, mouse_down],
      ['Human stable genes', 'Human up', 'Human down'],
      set_colors=['palegreen', 'red', 'blue'])
plt.show()
"""











