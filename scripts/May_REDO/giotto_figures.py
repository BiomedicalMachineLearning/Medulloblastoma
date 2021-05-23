""" Creating the giotto figures to display in the paper using stlearn
    visualisation.

    OUTPUT:

        figure_components/PseudoLimma/DE_figures/giotto/*
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
import scanpy as sc
import scripts.stlearn_helpers as st_hs
import stlearn as st

import beautifulcells.visualisation.helpers as vhs

out_plots = 'figure_components/PseudoLimma/DE_figures/'
score_dir = 'data/DE_out/Pseudo_TMM_Limma_Voom/giotto_out/'

hs.setUp()

################################################################################
                    # Loading in the data #
################################################################################
# AnnData files #
data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir+sample for sample in samples_]
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']

filt_files = np.array( os.listdir(data_dir+'ryan_ids/') )
filt_files = filt_files[np.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir+'ryan_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

datas = [st_hs.load(sample_dir, spot_filters[i], species_metas[i], subset=False)
         for i, sample_dir in enumerate(data_dirs)]

# The enrichment scores to add to these files #
df_names = os.listdir(score_dir)
giotto_scores = []
lims = [] #vmin, vmaxs
for df_name in df_names:
    score_df = pd.read_csv(score_dir+df_name, sep='\t')
    giotto_scores.extend(score_df.columns)
    lims.extend( [(min(score_df.values[:,i]), max(score_df.values[:,i]))
                  for i in range(score_df.shape[1])] )
    for i, data in enumerate(datas):
        data_bcs = np.array([bc for bc in score_df.index
                                if bc.startswith(samples[i])])
        bc_indices = [np.where(data.obs_names==bc.split('_')[-1])[0][0]
                      for bc in data_bcs]
        data_bcs = data_bcs[bc_indices]
        data_scores = score_df.loc[data_bcs,:]

        for j, col in enumerate(score_df.columns):
            data.obsm[col] = data_scores.values[:,j]

################################################################################
                    # Creating the enrichment plots #
################################################################################
lims[1] = (lims[1][0], 10)
for i in range(len(giotto_scores)):
    for j, data in enumerate(datas):
        vmin, vmax = lims[i]
        st.pl.het_plot(data, use_het=giotto_scores[i], vmin=vmin, vmax=vmax,
                       size=15)
        vhs.dealWithPlot(True, True, True, out_plots+'giotto/',
                         f'{samples[j]}_{giotto_scores[i]}.pdf', dpi=300)















