""" Performs the per-spot LR analysis of the visium data or each sample, &
    subsequently outputs visualisations for the top few LR pairs across
    the samples !
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import stlearn as st
import matplotlib.pyplot as plt

import scripts.stlearn_helpers as st_hs

################################################################################
                    # Loading the data #
################################################################################
data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'
os.chdir(data_dir)
os.chdir('../')
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

datas = [st_hs.load_and_norm_LRAnalysis(
              sample_dir, spot_filters[i], species_metas[i],
            min_cell_prop=0.01, subset=False)
         for i, sample_dir in enumerate(data_dirs)]

# Filtering mouse genes #
# for i, data in enumerate(datas):
#     human_genes = [gene for gene in data.var_names if 'hg38-' in gene]
#     datas[i] = data[:,human_genes]

################################################################################
                # Performing cci analysis #
################################################################################
# Load the NATMI literature-curated database of LR pairs, data formatted #
lrs_ = st.tl.cci.load_lrs(['connectomeDB2020_lit'])
# Testing for mice-data interactions !!!!
lrs = ['hg38-'+lr.replace('_', '_mm10-') for lr in lrs_]
lrs.extend( ['hg38-'+lr.replace('_', '_hg38-') for lr in lrs_] )
lrs = np.array(lrs+['mm10-'+lr.replace('_', '_hg38-') for lr in lrs_])

lr_summaries = []
for data in datas:
    st.tl.cci.run(data, lrs,
                  use_label = None, #Need to add the label transfer results to object first, above code puts into 'label_transfer'
                  use_het = 'cell_het', #Slot for cell het. results in adata.obsm, only if use_label
                  min_spots = 6, #Filter out any LR pairs with no scores for less than 5 spots
                  distance=60, #distance=0 for within-spot mode
                  n_pairs=0, #Number of random pairs to generate
                  adj_method='fdr_bh', #MHT correction method
                  lr_mid_dist = 150, #Controls how LR pairs grouped when creating bg distribs, higher number results in less groups
                                    #Recommended to re-run a few times with different values to ensure results robust to this parameter.
                  min_expr=0,
                  )
    lr_summaries.append(data.uns['lr_summary'])

################################################################################
                # Visualising the results #
################################################################################
# Now looking at the LR pair with the highest number of sig. spots #
best_lr = datas[0].uns['lr_summary'].index.values[7]
best_lrs = ['hg38-TIMP1_hg38-CD63', 'mm10-C3_hg38-CD46']#'hg38-HLA-A_hg38-APLP2'
for i, data in enumerate(datas):
    # Binary LR coexpression plot for all spots #
    # st.pl.lr_plot(data, best_lr, inner_size_prop=0.1, outer_mode='binary',
    #               pt_scale=10,
    #               use_label=None, show_image=True,
    #               sig_spots=False)
    # plt.show()
    for best_lr in [lr for lr in best_lrs if lr in data.uns['lr_summary'].index]:
        # Significance scores for all spots #
        st.pl.lr_plot(data, best_lr, inner_size_prop=1, outer_mode=None,
                      pt_scale=10,
                      use_label='lr_scores', show_image=True,
                      sig_spots=False)
        plt.savefig(f'figure_components/lr_prelim/{samples[i]}_{best_lr}_sig.png')

        # Binary LR coexpression plot for significant spots #
        # st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode='binary',
        #               pt_scale=20,
        #               use_label=None, show_image=True,
        #               sig_spots=True)
        # plt.show()

        # Continuous LR coexpression for signficant spots #
        st.pl.lr_plot(data, best_lr,
                      inner_size_prop=0.1, middle_size_prop=.15, outer_size_prop=.4,
                      outer_mode='binary', pt_scale=50,
                      use_label=None, show_image=True,
                      sig_spots=False)
        plt.savefig(f'figure_components/lr_prelim/{samples[i]}_{best_lr}_expr.png')

    # Continous LR coexpression for significant spots with tissue_type information #
    # st.pl.lr_plot(data, best_lr,
    #               inner_size_prop=0.08, middle_size_prop=.3, outer_size_prop=.5,
    #               outer_mode='continuous', pt_scale=150,
    #               use_label='tissue_type', show_image=True,
    #               sig_spots=True)
    # plt.show()

# Now looking at the LR pair with the highest number of sig. spots #
best_lr = datas[0].uns['lr_summary'].index.values[7]
for data in datas:
    # LR enrichment scores
    data.obsm[f'{best_lr}_scores'] = data.uns['per_lr_results'][best_lr].loc[:,
                                                                 'lr_scores'].values
    # -log10(p_adj) of LR enrichment scores
    data.obsm[f'{best_lr}_log-p_adj'] = data.uns['per_lr_results'][best_lr].loc[:,
                                                             '-log10(p_adj)'].values
    # Significant LR enrichment scores
    data.obsm[f'{best_lr}_sig-scores'] = data.uns['per_lr_results'][best_lr].loc[:,
                                                             'lr_sig_scores'].values

    # Visualising these results #
    # st.pl.het_plot(data, use_het=f'{best_lr}_scores', cell_alpha=0.7)
    # plt.show()

    st.pl.het_plot(data, use_het=f'{best_lr}_sig-scores', cell_alpha=0.7)
    plt.show()

# Saving the results #
# [data.write_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}.h5ad',compression='gzip')
#  for i, data in enumerate(datas)]



