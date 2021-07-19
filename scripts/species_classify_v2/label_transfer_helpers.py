"""
Helper functions for performing scanorama label transfer.
"""

import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import cosine_distances

import scanorama

import stlearn as st
import scanpy as sc

import beautifulcells.visualisation.helpers as vhs

def scan_label_transfer_analysis(data, sample, species, count_ads, exper_ids,
                                 out_dir, out_plots, **kwargs):
    """ Performs the integration & label transfer using scanorama.
    """
    int_data = scan_integration(data, sample, species,
                                      count_ads, exper_ids, out_plots,
                                      knn=50, alpha=.1)
    scan_label_transfer(int_data, data, sample,
                              species, out_dir, out_plots)
    return int_data

def scan_integration(data, sample, species, count_ads, exper_ids,
                                                          out_plots, **kwargs):
    """Performs scanorama integration.
    """
    to_int_datas = [data] + count_ads
    # combine multiple ST array with scanorama normalization
    corrected_M = scanorama.correct_scanpy(to_int_datas, return_dimred=True,
                                           **kwargs)

    # Visualising labels in integrated space #
    int_data = corrected_M[0].concatenate(corrected_M[1:],
                                          batch_key="sample",
                                          uns_merge="first",
                                          batch_categories=[sample]+exper_ids,
                                          )

    # Adding the raw #
    merged_data = to_int_datas[0].concatenate(to_int_datas[1:],
                                          batch_key="sample",
                                          uns_merge="first",
                                          batch_categories=[sample]+exper_ids,
                                          )
    int_data.raw = merged_data

    # Getting neighbours... #
    st.pp.neighbors(int_data, n_neighbors=25, use_rep='X_scanorama')
    sc.tl.umap(int_data)

    # Plotting the integration results #
    if 'hov_leiden' in int_data.obs:
        prefix=''
        color = ['hov_leiden', 'exper_id', 'NHLH1', 'species']
    elif 'cell_labels' in int_data.obs:
        prefix='vlad'
        color = ['cell_labels_merged', 'species']
    else:
        prefix='integrated'
        color = ['glial', 'vasculature']
    sc.pl.umap(int_data,
               color=color, show=False)
    vhs.dealWithPlot(True, True, True,
                     out_plots,
                     f'scanorama_{species}_{sample}_{prefix}umaps.pdf', 300)

    return int_data

def scan_label_transfer(int_data, data, sample, species, out_dir, out_plots,
                        label_key='hov_leiden'):
    """ Performs the scanorama label transfer.
    """
    vis_bool = int_data.obs['sample'].values.astype(str) == sample
    vis_data = int_data[vis_bool, :].copy()
    data.obsm['X_scanorama'] = vis_data.obsm['X_scanorama']
    data.obsm['X_umap'] = vis_data.obsm['X_umap']
    sc_data = int_data[vis_bool == False, :]

    scan_label_transfer_(data, sc_data, label_key)

    leiden_clusters = list(np.unique(sc_data.obs[label_key]))  # ['0', '1', '2', '3', '4', '5', '6']
    sc.pl.umap(data,
               color=leiden_clusters + ['species'], show=False)
    vhs.dealWithPlot(True, True, True,
                     out_plots,
                  f'scanorama_{species}_{sample}_{label_key}_label_transfer_umaps.pdf', 300)
    sc.pl.spatial(data,
                  color=leiden_clusters + ['species'],
                  library_id=f'Visium8_{sample}_Hybrid', show=False)
    vhs.dealWithPlot(True, True, True,
                     out_plots,
               f'scanorama_{species}_{sample}_{label_key}_label_transfer_spatials.pdf', 300)

    # Saving the label transfer values #
    label_transfer_df = data.obs[leiden_clusters]
    label_transfer_df.to_csv(
        out_dir + f'{sample}_{species}_{label_key}_label_transfer.txt', sep='\t')


def scan_label_transfer_(vis_data, sc_data, ref_key):
    """Performs scanorama label transfer."""

    # vis_bool = int_data.obs[batch_key].values.astype(str) == query_label
    # vis_scan = int_data[vis_bool, :].obsm['X_scanorama']
    # sc_scan = int_data[vis_bool == False, :].obsm['X_scanorama']
    vis_scan = vis_data.obsm['X_scanorama']
    sc_scan = sc_data.obsm['X_scanorama']
    sc_clusters = sc_data.obs[ref_key]

    distances = 1 - cosine_distances(sc_scan, vis_scan)

    # def divide(vals, denom):
    #     return vals / denom
    #
    # def label_transfer(dist, labels):
    #     lab = pd.get_dummies(labels).to_numpy().T
    #     class_prob = lab @ dist
    #     norm = np.linalg.norm(class_prob, 2, axis=0)
    #     class_prob = class_prob / norm
    #     class_prob = (class_prob.T - class_prob.min(1))  # / class_prob.ptp(1)
    #     totals = class_prob.sum(axis=1)
    #     class_prob = np.apply_along_axis(divide, 0, class_prob, totals)
    #     return class_prob
    def label_transfer(dist, labels):
        """The original version."""
        lab = pd.get_dummies(labels).to_numpy().T
        class_prob = lab @ dist
        norm = np.linalg.norm(class_prob, 2, axis=0)
        class_prob = class_prob / norm
        class_prob = (class_prob.T - class_prob.min(1)) / class_prob.ptp(1)
        return class_prob

    class_probs = label_transfer(distances, sc_clusters)

    class_probs_df = pd.DataFrame(
        class_probs,
        columns=np.sort(sc_data.obs[ref_key].unique()),
        index=vis_data.obs_names
    )

    vis_data.obs = pd.concat([vis_data.obs, class_probs_df], axis=1)
    print("Added label transfer results to columns: ",
          np.sort(sc_data.obs[ref_key].unique()))

























