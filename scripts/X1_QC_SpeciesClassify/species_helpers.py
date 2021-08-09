""" Helper functions for the species classify
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import scripts.utils.visualisation.helpers as vhs

def add_species_counts(data):
    """Adds the species counts to the anndata.
    """
    human_genes = np.array([gene for gene in data.var_names
                            if 'hg38-' in gene])
    mouse_genes = np.array([gene for gene in data.var_names
                            if 'mm10-' in gene])

    human_bool = data[:, human_genes].to_df().transpose().values > 0
    mouse_bool = data[:, mouse_genes].to_df().transpose().values > 0

    human_gene_counts = human_bool.sum(axis=0)
    mouse_gene_counts = mouse_bool.sum(axis=0)

    data.obs['human_gene_counts'] = human_gene_counts
    data.obs['mouse_gene_counts'] = mouse_gene_counts

    return data

def plot_species_scatters(species_labels, human_counts, mouse_counts, colors,
                          figsize=(3, 3), alpha=.5, s=10, linewidths=3,
                          out_path=None):
    """ Purpose is to plot the species scatter plots. """

    matplotlib.rcParams.update({'font.size': 8, 'font.weight': 'bold'})
    fig, ax = plt.subplots(figsize=figsize)
    for label in colors:
        label_bool = species_labels == label
        ax.scatter(mouse_counts[label_bool], human_counts[label_bool],
                    c=colors[label], alpha=alpha, s=s, linewidths=linewidths)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    if type(out_path)!=type(None):
        out_dir = '/'.join(out_path.split('/')[0:-1]) + '/'
        out_file = out_path.split('/')[-1]
        print(out_dir, out_file)
        vhs.dealWithPlot(True, True, True, out_dir, out_file, 300)














