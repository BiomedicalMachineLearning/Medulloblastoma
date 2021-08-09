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

import scripts.utils.visualisation.helpers as vhs

file_dir = os.path.realpath(__file__)
project_dir = file_dir.split('scripts/')[0]
data_dir = project_dir+'data/'

def setUp():
    os.chdir(project_dir)

def load(sample_dir, prefix):
    """Loads in the data and filters for QC and data/mix spots"""
    data = st.Read10X(sample_dir)
    data.var_names_make_unique()
    data.var_names = numpy.array([var_name.replace('_', '-')
                                  for var_name in data.var_names])
    data.obs_names = numpy.array([prefix+obs_name
                                  for obs_name in data.obs_names])

    return data

def count_true(vals):
    return len(np.where(vals)[0])

def plot_species_scatters(species_labels, human_counts, mouse_counts, colors,
                          figsize=(3, 3), alpha=.5, s=10, linewidths=3,
                          out_path=None, ymax=None, xlim=None, **kwargs):
    """ Purpose is to plot the species scatter plots. """

    matplotlib.rcParams.update({'font.size': 8, 'font.weight': 'bold'})
    fig, ax = plt.subplots(figsize=figsize)
    for label in colors:
        label_bool = species_labels == label
        ax.scatter(mouse_counts[label_bool], human_counts[label_bool],
                    c=colors[label], alpha=alpha, s=s, linewidths=linewidths,
                   **kwargs)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    if type(ymax)!=type(None):
        ax.set_ylim((0, ymax))
    if type(xlim)!=type(None):
        ax.set_xlim(xlim)

    if type(out_path)!=type(None):
        out_dir = '/'.join(out_path.split('/')[0:-1]) + '/'
        out_file = out_path.split('/')[-1]
        print(out_dir, out_file)
        vhs.dealWithPlot(True, True, True, out_dir, out_file, 300)








