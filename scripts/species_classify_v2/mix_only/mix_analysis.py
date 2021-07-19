"""
mix_analysis.py (07/07/2021) -> Purpose of this script is to load only the
                         mix spots, but separate into human genes/mouse genes.
                         For the mix-humang, compare to Hovestadt reference.
                         For the mix-mouseg, compare to Vladoiu reference.
                         Going to do this for each sample independently, but then
                         add the scores from the comparison to the origin mix
                         anndata, & subsequently integrate the mix spots into a
                         common space. Thereafter, will examine the scores &
                         cluster, labelling the clusters according to the
                         cell type combinations.

                         INPUT: data/Visium8*/
                         OUTPUT: figure_components/species_classify_v2/mix_only/
                                 data/scanpy_h5ads/mix_int_data.h5ad
"""

################################################################################
                            # Environment setup #
################################################################################

import os
import scanpy as sc
import pandas as pd
import numpy
import numpy as np
import matplotlib.pyplot as plt

from anndata import AnnData
import scripts.helpers as hs

hs.setUp()
import scripts.stlearn_helpers as st_hs
import scripts.species_classify_v2.label_transfer_helpers as lt_hs
import scripts.species_classify_v2.mix_only.mix_helpers as mix_hs

out_dir = 'data/scanpy_h5ads/'
out_plots = 'figure_components/species_classify_v2/mix_only/'
hov_dir = 'data/third_party_data/Hovestadt2019_Nature_scRNA/'
vlad_dir = 'data/third_party_data/Vladoiu2019_Nature_scRNA/'

################################################################################
# Loading data & normalising #
################################################################################
# read in visium dataset downloaded from: support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Breast_Cancer_Block_A_Section_2
data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir + sample for sample in samples_]
species_metas = ['species_classify_v2/A1_species.txt',
                 'species_classify_v2/B1_species.txt',
                 'species_classify_v2/C1_species.txt',
                 'species_classify_v2/D1_species.txt']

filt_files = numpy.array(os.listdir(data_dir + 'ryan_ids/'))
filt_files = filt_files[numpy.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir + 'ryan_ids/' + filt_file,
                            header=None).values[:, 0] for filt_file in
                filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

# Loading only the the data that pass QC without subsetting #
datas = [st_hs.load_and_norm(sample_dir, spot_filters[i], species_metas[i],
                             subset=True, sme=False,
                             subset_labels='mix')
         for i, sample_dir in enumerate(data_dirs)]

###### Subsetting to just border spots ######
row_coords = [[[20, 60]], [[0, 100], [62, 100], [30, 50], [0, 15]], #treat
              [[0, 45], [40,77], [34,63], [20,30], [55,100]], #C1
              [[0,45],[27,45], [15,45], [38, 58], [58,75], [65,75]]] #D1
col_coords = [[[0, 50]], [[0, 35], [0, 200], [0, 50], [70, 200]], #treat
              [[30, 112], [20,66], [2,40], [110,200], [95,150]], #C1
              [[80,200], [30,200], [66,200], [26,76], [26,62], [10,26]]] #D1
for i in range(len(datas)):
    sc.pl.spatial(datas[i], color=['array_row', 'array_col'])
    keep_bool = np.array([True] * len(datas[i]))
    for j in range(len(row_coords[i])):
        row_bool = np.logical_and(
            datas[i].obs['array_row'].values > row_coords[i][j][0],
            datas[i].obs['array_row'].values < row_coords[i][j][1])
        col_bool = np.logical_and(
            datas[i].obs['array_col'].values > col_coords[i][j][0],
            datas[i].obs['array_col'].values < col_coords[i][j][1])
        outlier_bool = np.logical_and(row_bool, col_bool)
        keep_bool = np.logical_and(keep_bool, outlier_bool == False)

        sc.pl.spatial(datas[i][keep_bool, :], color=['array_row', 'array_col'])

    datas[i] = datas[i][keep_bool,:]

# Saving the data for faster reloading #
# [data.write_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}_all_species_SME.h5ad',
#                  compression='gzip') for i, data in enumerate(datas)]

# Reformatting each so just containing data genes/spots from the relevant species
# in the right formatting #
datas_human = []
datas_mouse = []
for i, data in enumerate(datas):
    datas_human.append(st_hs.species_split(data, species='mix-humang'))
    datas_mouse.append(st_hs.species_split(data, species='mix-mouseg'))

################################################################################
# Loading the Hovestadt reference data #
################################################################################
hov_data = sc.read_h5ad(hov_dir + 'hov_scrna.h5ad')  # This was already integrated
# with the labels.

count_files = [file_name for file_name in os.listdir(hov_dir)
               if file_name[0:3] == 'GSM']

count_ads = []
exper_ids = []
for file_name in count_files:
    count_df = pd.read_csv(hov_dir + file_name, sep='\t', index_col=0)

    exper_id = file_name.split('_')[1].replace('.txt', '')
    exper_ids.append(exper_id)
    exper_type = 'PDX' if 'RCMB' in file_name else 'patient'
    meta_df = pd.DataFrame([[exper_id] * count_df.shape[1],
                            [exper_type] * count_df.shape[1]]).transpose()
    meta_df.index = count_df.columns
    meta_df.columns = ['exper_id', 'exper_type']

    count_ad = AnnData(count_df.transpose(), obs=meta_df)
    bcs = [f'{bc}-{exper_id}' for bc in count_ad.obs_names]
    count_ad.obs['hov_leiden'] = hov_data[bcs, :].obs['leiden'].values.astype(
        str)
    count_ads.append(count_ad)

# Removing cluster 6 since is damaged/dying cells #
hov_data = hov_data[hov_data.obs['leiden'].values!='6',:]

for i, count_ad in enumerate(count_ads):
    count_ads[i] = count_ad[count_ad.obs['hov_leiden'].values!='6',:]

################################################################################
# Loading the Vladoiu reference data #
################################################################################
vlad_data = sc.read_h5ad(vlad_dir + 'Vladoiu2019_logcpms.h5ad')

# Subsetting to just the mature cells #
vlad_data_sub = vlad_data[vlad_data.obs['time'].values.astype(str) == 'P7', :]
cell_labels_sub = vlad_data_sub.obs['cell_labels'].values.astype(str)
print(np.unique(cell_labels_sub))
remove_cells = ['Mesenchymal stem cells-1',
                'Unipolar brush cell and GCP progenitor',
                'VZ progenitors']
cell_bool = [label not in remove_cells for label in cell_labels_sub]
mature_cells = np.unique(cell_labels_sub[cell_bool])

cell_labels = vlad_data.obs['cell_labels'].values.astype(str)
mature_cell_bool = [label in mature_cells for label in cell_labels]
vlad_data = vlad_data[mature_cell_bool, :]

# grouping cells #
cell_labels = vlad_data.obs['cell_labels'].values.astype(str)
cell_labels = mix_hs.vlaodiu_label_merge(cell_labels, np.unique(cell_labels))

vlad_data.obs['cell_labels_merged'] = cell_labels

# Just visualising the data #
sc.pp.highly_variable_genes(vlad_data,
                            min_disp=.4, min_mean=0.02, max_mean=3)
sc.pl.highly_variable_genes(vlad_data)
print(sum(vlad_data.var['highly_variable']))

sc.pp.pca(vlad_data, n_comps=100)
sc.pp.neighbors(vlad_data, n_neighbors=15)
sc.tl.umap(vlad_data)

sc.pl.umap(vlad_data, color=['cell_labels_merged', 'time', 'cell_labels', ])

vlad_data.write_h5ad(vlad_dir + 'Vladoiu2019_lateRef.h5ad', compression='gzip')

################################################################################
# Scanorama integration for each dataset independently for human ref. #
################################################################################
species = 'human'
datas_species = datas_human  # if species=='human' else datas_mouse
reload = True
params = [{'knn': 20, 'alpha': .05}, {'knn': 20, 'alpha': .05},
          {'knn': 20, 'alpha': .05}, {'knn': 20, 'alpha': .05}]
for i in range(len(datas_species)):
    if reload:
        datas_species[i] = st_hs.species_split(datas[i], species=species)
    int_data = lt_hs.scan_integration(datas_species[i], samples[i], species,
                                      count_ads, exper_ids, out_plots,
                                      knn=params[i]['knn'],
                                      alpha=params[i]['alpha'])
    lt_hs.scan_label_transfer(int_data, datas_species[i], samples[i],
                              species, out_dir, out_plots)

################################################################################
# Scanorama integration for each dataset independently for mouse ref. #
################################################################################
species = 'mouse'
datas_species = datas_mouse  # if species=='human' else datas_mouse
reload = True
params = [{'knn': 20, 'alpha': .05}, {'knn': 20, 'alpha': .05},
          {'knn': 20, 'alpha': .05}, {'knn': 20, 'alpha': .05}]
for i in range(len(datas_species)):
    if reload:
        datas_species[i] = st_hs.species_split(datas[i], species=species)
    int_data = lt_hs.scan_integration(datas_species[i], samples[i], species,
                                      [vlad_data], ['vladoiu'], out_plots,
                                      knn=params[i]['knn'],
                                      alpha=params[i]['alpha'])
    lt_hs.scan_label_transfer(int_data, datas_species[i], samples[i],
                              species, out_dir, out_plots,
                              label_key='cell_labels_merged')

################################################################################
# Pulling the information from the .obs into the original datas #
################################################################################
for i, data in enumerate(datas):
    orig_obs = data.obs.copy()
    h_obs = datas_human[i].obs
    h_obs = h_obs.loc[:, [col for col in h_obs.columns
                          if col not in orig_obs.columns]]
    m_obs = datas_mouse[i].obs
    m_obs = m_obs.loc[:, [col for col in m_obs.columns
                          if col not in orig_obs.columns and
                          col not in h_obs.columns]]
    obs_wHuman = pd.concat([orig_obs, h_obs], axis=1)
    obs_wHwM = pd.concat([obs_wHuman, m_obs], axis=1)
    data.obs = obs_wHwM

################################################################################
                            # Merging the datas #
################################################################################
int_data = lt_hs.scan_integration(datas[0], samples[0], 'mix',
                                  datas[1:], samples[1:], out_plots)

sc.tl.leiden(int_data, resolution=.3)

sc.pl.umap(int_data, color=['leiden', 'glial', 'vasculature', 'Granule cells',
                            'Purkinje cells', '0', '2', '4'])

samp = 'D1'
sc.pl.umap(int_data[int_data.obs['sample'].values == samp, :],
           color=['leiden', 'glial', 'vasculature', 'Purkinje cells',
                  'Granule cells', '0', '1', '2', '3', '4', '5'])

sc.pl.spatial(int_data[int_data.obs['sample'].values == samp, :],
              color=['leiden', 'glial', 'vasculature',
                     '0', '1', '3', '4'],
              library_id=f'Visium8_{samp}_Hybrid')

# Saving these results #
int_data.write_h5ad(f'{data_dir}scanpy_h5ads/integrated_border.h5ad',
                  compression='gzip')

################################################################################
        # Measuring the correlations between the two sets of scores #
################################################################################
from sklearn.metrics.pairwise import pairwise_distances
from scipy.stats import spearmanr

def spearman(vals1, vals2):
    rho, p = spearmanr(vals1, vals2)
    return rho

hov_set = np.unique(hov_data.obs['leiden'].values.astype(str))
vlad_set = np.unique(vlad_data.obs['cell_labels_merged'].values.astype(str))

all_set = list(hov_set) + list(vlad_set)
all_scores = int_data.obs.loc[:, all_set].values.transpose()

metric = spearman

distances = pairwise_distances(all_scores, all_scores, metric=metric)
dists_all = pd.DataFrame(distances, index=all_set, columns=all_set)

data_dists_comp = {}
data_dists_human = {}
data_dists_mouse = {}
for i in range(len(datas_human)):
    data_human = datas_human[i]
    data_mouse = datas_mouse[i]
    human_scores = data_human.obs.loc[:, hov_set].values.transpose()
    mouse_scores = data_mouse.obs.loc[:, vlad_set].values.transpose()

    distances = pairwise_distances(human_scores, mouse_scores, metric=metric)
    dists_comp = pd.DataFrame(distances, index=hov_set, columns=vlad_set)

    distances = pairwise_distances(human_scores, human_scores, metric=metric)
    dists_human = pd.DataFrame(distances, index=hov_set, columns=hov_set)

    distances = pairwise_distances(mouse_scores, mouse_scores, metric=metric)
    dists_mouse = pd.DataFrame(distances, index=vlad_set, columns=vlad_set)

    data_dists_comp[samples[i]] = dists_comp
    data_dists_human[samples[i]] = dists_human
    data_dists_mouse[samples[i]] = dists_mouse













