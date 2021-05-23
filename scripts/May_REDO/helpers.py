"""
Helper functions for reading/writing the data.
"""

import scanpy as sc
import stlearn as st
import numpy, pandas
import numpy as np
import pandas as pd
import os, sys
import matplotlib
import matplotlib.pyplot as plt

import beautifulcells.visualisation.helpers as vhs

data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'
data_dir2 = '/Volumes/GML001-Q1851/Quan/Visium/Visium19_MB/'

def setUp():
    os.chdir('/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/')

def load(sample_dir, prefix):
    """Loads in the data and filters for QC and human/mix spots"""
    data = st.Read10X(sample_dir)
    data.var_names_make_unique()
    data.var_names = numpy.array([var_name.replace('_', '-')
                                  for var_name in data.var_names])
    data.obs_names = numpy.array([prefix+obs_name
                                  for obs_name in data.obs_names])

    return data

def getGenePercents(data, genes, n_expr=True):
    """n_expr=True if in terms of boolean numbers, or in-terms of total counts \
    per spot to calculate the percentage.
    """

    df = data.to_df()
    df_sub = df.loc[:, genes]

    if n_expr:
        n_genes = (df.values > 0).sum(axis=1)
        n_sub = (df_sub.values > 0).sum(axis=1)
        perc = n_sub/n_genes

    else:
        n_genes = (df.values).sum(axis=1)
        n_sub = (df_sub.values).sum(axis=1)
        perc = n_sub / n_genes

    return perc

def count_true(vals):
    return len(np.where(vals)[0])

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








