"""
Loads in the data, performs SME normalisation, & saves output.

    INPUT: * data/Visium8*/
    OUTPUT: * data/scanpy_h5ads/*
"""


################################################################################
                    # Environment setup #
################################################################################

import os
import pandas as pd
import numpy as np

import scripts.utils.helpers as hs
hs.setUp()
import scripts.utils.stlearn_helpers as st_hs

################################################################################
                    # Loading data & normalising #
################################################################################
data_dir = hs.data_dir
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir+sample for sample in samples_]
species_metas = ['species_classify_v2/A1_species.txt',
                 'species_classify_v2/B1_species.txt',
                 'species_classify_v2/C1_species.txt',
                 'species_classify_v2/D1_species.txt']

filt_files = np.array( os.listdir(data_dir+'filter_ids/') )
filt_files = filt_files[np.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir+'filter_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

# Loading all the data that pass QC without subsetting #
datas = [st_hs.load_and_norm(sample_dir, spot_filters[i], species_metas[i],
                             subset=True, sme=True,
                             subset_labels='human_mix_mouse')
         for i, sample_dir in enumerate(data_dirs)]

# Saving the data #
[data.write_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}_all_species_SME.h5ad',
                 compression='gzip') for i, data in enumerate(datas)]


