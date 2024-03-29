"""
Creates the spatial plots of the per-spot enrichment in a consistent way across
the samples for comparison.

        INPUT: * data/Visium8*/
               * data/filter_ids/*
               * data/giotto_out/SHH_gsea_v2_scores.txt

        OUTPUT: * figure_components/giotto_figures/*
"""

################################################################################
                    # Environment setup #
################################################################################

import os, sys
import numpy as np
import pandas as pd
import scripts.utils.helpers as hs
import scripts.utils.stlearn_helpers as st_hs
import scripts.utils.visualisation.helpers as vhs
import stlearn as st

score_dir = 'data/giotto_out/'
out_plots = 'figure_components/giotto_figures/'

hs.setUp()

################################################################################
                    # Loading in the data #
################################################################################
# AnnData files #
data_dir = hs.data_dir
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir+sample for sample in samples_]
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']

filt_files = np.array( os.listdir(data_dir+'filter_ids/') )
filt_files = filt_files[np.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir+'filter_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

datas = [st_hs.load(sample_dir, spot_filters[i], species_metas[i], subset=False)
         for i, sample_dir in enumerate(data_dirs)]

# The enrichment scores to add to these files #
df_names = ['SHH_gsea_v2_scores.txt']
species = '' if 'mouse' not in df_names[0] else 'mouse'
giotto_scores = []
lims = [] #vmin, vmaxs
upper_lim_scale = {'SHH.B enrich scores': 1,
                   'HALLMARK_E2F_TARGETS enrich scores': .1}
for df_name in df_names:
    score_df = pd.read_csv(score_dir+df_name, sep='\t')
    giotto_scores.extend(score_df.columns)
    scales = [.2 if score_name not in upper_lim_scale 
              else upper_lim_scale[score_name] 
              for score_name in score_df.columns]
    lims.extend( [(min(score_df.values[:,i]), max(score_df.values[:,i])*scales[i])
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
for i in range(len(giotto_scores)):
    vmin, vmax = lims[i]
    for j, data in enumerate(datas):
        st.pl.het_plot(data, use_het=giotto_scores[i], vmin=vmin, vmax=vmax,
                       size=15, cmap='plasma', cell_alpha=.65)
        vhs.dealWithPlot(True, True, True, out_plots,
                      f'{samples[j]}_{giotto_scores[i]}_{species}.pdf', dpi=300)















