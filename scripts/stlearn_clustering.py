""" Created: 26th March 2021

Purpose is to run cluster within each of the samples, then compare these to
    pseudobulked reference profiles of developing cerebellum cell types
    in order to biologically annotate each of the medulloblastoma clusters.

"""

################################################################################
                    # Environment setup #
################################################################################

import os
import stlearn as st
import scanpy as sc
import pandas as pd
import pandas
import numpy
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.spatial.distance import cosine
from scipy.stats import zscore

from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering

import scripts.helpers as hs
hs.setUp()
import scripts.stlearn_helpers as st_hs
import beautifulcells.clustering.marker_labelling.marker_labelling \
                                                             as marker_labelling
import beautifulcells.normalisation.Discretize.Discretize as discretize

import beautifulcells.visualisation.heatmap.heatmap_helpers as heat_hs
import beautifulcells.visualisation.helpers as vhs

import beautifulcells.preprocessing.load_data.simple_pickle as simple_pickle

################################################################################
                    # Loading data & normalising #
################################################################################
# read in visium dataset downloaded from: support.10xgenomics.com/spatial-gene-expression/datasets/1.0.0/V1_Breast_Cancer_Block_A_Section_2
data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'
samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
data_dirs = [data_dir+sample for sample in samples_]
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']

filt_files = numpy.array( os.listdir(data_dir+'ryan_ids/') )
filt_files = filt_files[numpy.argsort(filt_files)]
spot_filters = [pd.read_csv(data_dir+'ryan_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]

samples = ['A1', 'B1', 'C1', 'D1']

datas = [st_hs.load_and_norm(sample_dir, spot_filters[i], species_metas[i])
         for i, sample_dir in enumerate(data_dirs)]

# Saving these to improve loading speed #
[data.write_h5ad(f'{data_dir}scanpy_h5ads/{samples[i]}.h5ad',compression='gzip')
 for i, data in enumerate(datas)]

################################################################################
            # Clustering each, with the goal of getting relatively
            # comparable clusters
################################################################################
##### Clustering ######
#resolutions = {0: 2.3, 1: 3, 2: 2.5, 3: 1.5} # from when including mouse spots
resolutions = {0: 1, 1: 1, 2: 1, 3: 1}
for i in range(0, 4):
    res = resolutions[i]
    st.tl.clustering.louvain(datas[i], resolution=res)
    st.pl.cluster_plot(datas[i], dpi=100, use_label='louvain')
    plt.show()
""" Cross referenced the above clusters with the SHHA/B/C scores, and in the 
    treated samples there is a tight correspondence !
"""

################################################################################
            ##### Assigning to equivalent clusters ######
################################################################################
####### Comparing the pseudo-bulked clusters between samples to get equivalent
####### cluster labels.
# Adding in the species information, subset to data/mix, cluster pseudobulk #
human_datas = []
dfs = []
obs_names = []
labels = []
for i, data in enumerate(datas):
    # Retrieving df and pseudo-bulking #
    df = data.to_df().transpose()
    human_genes = hs.getSpeciesGenes(df, 'hg38-')
    df = df.loc[human_genes, :]
    dfs.append( df )
    # The sample/cluster labels #
    sample = species_metas[i].split('_')[0]
    labels.extend( [f'{sample}_{cluster}'
                    for cluster in data.obs['louvain'].values] )
    obs_names.extend( df.columns.values )

labels = numpy.array(labels)
label_set = []
for sample in samples:
    sample_labels = labels[[sample in label for label in labels]]
    sample_labels = numpy.unique(sample_labels)
    clusters = numpy.array([int(label.split('_')[1]) for label in sample_labels])
    clusters_ordered = clusters[numpy.argsort(clusters)]
    label_set.extend([sample+'_'+str(cluster) for cluster in clusters_ordered])

count_datas = [st_hs.load(sample_dir, spot_filters[i], species_metas[i])
                                      for i, sample_dir in enumerate(data_dirs)]
count_dfs = [data.to_df().transpose() for data in count_datas]

# Concatenating the datas & pseudobulking #
data_df = pandas.concat(count_dfs, axis=1)
nan_bool = numpy.apply_along_axis(numpy.isnan, 1, data_df.values)
nan_bool = numpy.apply_along_axis(numpy.any, 1, nan_bool)
data_df = data_df.loc[nan_bool==False, :]

labels = numpy.array(labels)
# Saving the clusters #
cluster_meta = pandas.DataFrame(labels, index=data_df.columns,
                                columns=['clusters'])
cluster_meta.to_csv(data_dir+'spot_meta/clusters_human_mix.txt', sep='\t')

human_genes = hs.getSpeciesGenes(data_df, 'hg38-')
data_df_human = data_df.loc[human_genes, :]
data_df_human.index = [gene.replace('hg38-', '')
                       for gene in data_df_human.index]
data_df_human = heat_hs.formatData(data_df_human, labels, label_set, 'sum')[0]

## Loading in the pseudobulked data from the developing mouse cerebellum
## dataset to measure similarity against that as a reference
ref = pandas.read_csv('data/third_party_data/Vladoiu2019_Nature_scRNA/'
                      'Vladoiu2019_pseudo_counts.txt', sep='\t', index_col=0)
# Converting gene names to data format #
ref_genes = [gene.upper() for gene in ref.index]
ref.index = ref_genes

overlap = [gene for gene in data_df_human.index
           if gene in ref_genes and gene!='PISD']

ref = ref.loc[overlap, :]
data_df_human = data_df_human.loc[overlap, :]

# normalising #
ref_norm = st_hs.norm_df(ref, scale=True)
human_norm = st_hs.norm_df(data_df_human, scale=True)

# Getting cosine similarity between each tumor cluster & mouse cell types
sims = numpy.zeros((data_df_human.shape[1], ref.shape[1]))
for i in range(sims.shape[0]):
    tumor_expr = human_norm.values[:, i]
    for j in range(sims.shape[1]):
        cell_expr = ref_norm.values[:, j]
        sims[i, j] = 1-cosine(tumor_expr, cell_expr)
sims = pandas.DataFrame(sims, index=data_df_human.columns,
                        columns=ref.columns)
sims.to_csv(data_dir+'spot_meta/clusters_human_mix_vladiouSims.txt', sep='\t')
"""
sims = pandas.read_csv(data_dir+'spot_meta/clusters_human_mix_vladiouSims.txt',
                      sep='\t', index_col=0)
"""

######### Based on the above, assigning cluster labels #########################
# for each cluster, getting top 4 scores and their values #
cluster_scores = {}
for i, cluster in enumerate(sims.index.values):
    scores = sims.values[i,:]
    order = numpy.argsort(-scores)
    top_ = sims.columns.values[order[0:6]]
    top_scores = sims.values[i, order[0:6]]
    cluster_scores[cluster] = {}
    for j, top in enumerate(top_):
        cluster_scores[cluster][top] = top_scores[j]

cluster_to_type = {'A1_0': 'neuron-like',
                   'A1_1': 'differentiated-neural-like',
                   'A1_2': 'mouse-interface',
                   'A1_3': 'neuron-stem-like',
                   'A1_4': 'neuron-like',
                   'A1_5': 'glia-stem-like',
                   'A1_6': 'proliferative-stem-like', # End A
                   'B1_0': 'neuron-like',
                   'B1_1': 'glia-stem-like',
                   'B1_2': 'differentiated-neural-like',
                   'B1_3': 'mouse-interface',
                   'B1_4': 'meninges/vasculogenic-mimicry',
                   'B1_5': 'mouse-interface',
                   'B1_6': 'proliferative-stem-like',
                   'B1_7': 'differentiated-neural-like',
                   'B1_8': 'meninges/vasculogenic-mimicry', # End B
                   'C1_0': 'proliferative-stem-like',
                   'C1_1': 'mesenchymal-stem-like',
                   'C1_2': 'vasculogenic-stem-like',
                   'C1_3': 'differentiated-neural-like',
                   'C1_4': 'vasculogenic-stem-like',
                   'C1_5': 'mouse-interface', #NOTE: Very similiar to vasculature mimicry
                   'C1_6': 'vasculogenic-mimicry',
                   'C1_7': 'proliferative-stem-like',
                   'C1_8': 'mouse-interface', #NOTE: Very similiar to vasculature mimicry, unlike treated mouse interface
                   'D1_0': 'proliferative-stem-like',
                   'D1_1': 'mesenchymal-stem-like',
                   'D1_2': 'proliferative-stem-like',
                   'D1_3': 'vasculogenic-mimicry',
                   'D1_4': 'mouse-interface', #NOTE: Very similiar to vasculature mimicry..
                   'D1_5': 'mouse-interface',
                   'D1_6': 'mouse-interface',
                   'D1_7': 'vasculogenic-mimicry',
                   'D1_8': 'differentiated-neural-like',
                   'D1_9': 'mouse-interface'
                   }

# If updating using the pre-saved file #
"""
human_meta = pandas.read_csv(data_dir+'spot_meta/clusters_human_mix.txt',
                                                        sep='\t', index_col=0)
labels = human_meta.values[:,0]
barcodes = human_meta.index
new_labels = human_meta.values[:, 1]
"""

barcodes = data_df.columns.values
barcodes = [cluster.split('_')[0]+'-'+barcodes[i]
            for i, cluster in enumerate(labels)]

new_labels = [cluster_to_type[label] for label in labels]

human_meta = pandas.DataFrame([labels, new_labels]).transpose()
human_meta.columns = ['cluster', 'tissue_type']
human_meta.index = barcodes
human_meta.to_csv(data_dir+'spot_meta/clusters_human_mix.txt', sep='\t')

# Adding in the respective annotations for visualisation #
# Fast loading to just plot the clusters #
"""
datas = [st_hs.load_and_norm(sample_dir, spot_filters[i], species_metas[i],
                             sme=False)
         for i, sample_dir in enumerate(data_dirs)]
"""
import stlearn.plotting.classes as cl

# Defining colors #
colors = vhs.getColors(new_labels)
colors['proliferative-stem-like'] = 'springgreen'
colors['differentiated-neural-like'] = 'magenta'
colors['vasculogenic-mimicry'] = 'crimson'
colors['vasculogenic-stem-like'] = 'pink'
#simple_pickle.saveAsPickle(data_dir+'spot_meta/colors.pkl', colors)
#colors = simple_pickle.loadPickle(data_dir+'spot_meta/colors.pkl')

sizes=[6,6,7,7]
for i, data in enumerate(datas):
    obs_names = [samples[i]+'-'+obs for obs in data.obs_names]
    tissue_labels = human_meta.loc[obs_names,:].values[:,1]
    cluster_labels = human_meta.loc[obs_names,:].values[:,0]
    data.obs['tissue_type'] = tissue_labels
    data.obs['tissue_type'] = data.obs['tissue_type'].astype('category')
    data.obs['cluster'] = cluster_labels
    data.obs['cluster'] = data.obs['cluster'].astype('category')

    # Adding colors #
    data.uns['tissue_type_colors'] = []
    for j, cluster_tuple in enumerate(data.obs.groupby('tissue_type')):
        cluster = cluster_tuple[0]
        data.uns['tissue_type_colors'].append( colors[cluster] )

    x = cl.ClusterPlot(data, dpi=300, use_label='tissue_type',
                                        figsize=(10, 7), size=sizes[i])
    vhs.dealWithPlot(True, True, True,
                     data_dir+'../figure_components/stlearn_clustering/',
                     samples[i]+'_tissue_types.pdf', 300)

    st.pl.cluster_plot(data, dpi=300, use_label='cluster', figsize=(7.8, 6),
                       size=10)
    vhs.dealWithPlot(True, True, True,
                     data_dir + '../figure_components/stlearn_clustering/',
                     samples[i] + '_cluster.pdf', 300)

#### Make a table summarising the top scores for each. ######
score_table = numpy.empty((sims.shape[0], 7), object)
for i, cluster in enumerate(sims.index.values):
    cluster_label = cluster_to_type[cluster]
    score_table[i, 0] = cluster_label

    scores = sims.values[i,:]
    order = numpy.argsort(-scores)
    top_ = sims.columns.values[order[0:6]]
    top_scores = sims.values[i, order[0:6]]
    for j, top in enumerate(top_):
        score_table[i, j+1] = f'{top}: {top_scores[j].round(3)}'

score_df = pandas.DataFrame(score_table, index=sims.index,
                            columns=['tissue_type']+
                                    [f'Top cell type {i}' for i in range(1,7)])
score_df.to_excel(data_dir+'supps/cluster_to_tissues.xlsx',
                  "cluster-celltype-cosines")

# References:
# * granule cells:
#   https://www.frontiersin.org/articles/10.3389/fncir.2020.611841/full
# * EMT & vascular mimicry
#   https://pubmed.ncbi.nlm.nih.gov/26558950/

""" NOTE: * barcodes not random between visium slides !!!!
"""











################################################################################
                        # Junk Code #
################################################################################
"""
Ref of some CCIs in MDB:

Immunotherapy in Medullablastoma, contains L-R pairs:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7182450/

NOTE: NKG2D encoded by KLRK1 -
        is a receptor NK cells important in MDB immuno-resistance.

        -> notably expressed by both the tumour cells & mouse

    Important ligands of the above, which are present on MDB cell surface:
    * ULBP2 (found in the tumor genes, not in mouse)
    * MICA (found only in the data)
"""
z_scores = pandas.DataFrame(zscore(sims.values, 1),
                               index=data_df_human.columns, columns=ref.columns)

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(distance_threshold=0, n_clusters=None,
                                affinity='euclidean', linkage='ward')
model = model.fit(sims)

plt.title('Hierarchical Clustering Dendrogram')
# plot the top three levels of the dendrogram
plot_dendrogram(model, truncate_mode='level', #p=0,
                labels=human_norm.columns.values)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()


######## Criteria based, whereby call DE genes using t-tests and use these to
######## generate interpretable criteria for different tissue types
### DE gene testing ###
sc.tl.rank_genes_groups(datas[0], 'louvain', method='wilcoxon')
de_genes = pd.DataFrame(datas[0].uns['rank_genes_groups']['names'])

st.pl.gene_plot(datas[0], gene_symbols='hg38-CD63')
plt.show()

### binary-based labelling ###
marker_df = marker_labelling.readMarkerFile(data_dir +
                                            'cluster_criteria/mdb_clusters.txt')

params = {'classifier_params': {'C': .1, 'penalty': 'l1', 'solver': 'saga',
                                'tol': 1e-4, 'max_iter': 300},
          'n_top': 300, 'scale': True}
binary_params= {'method': 'quantile_non_zero', 'cutoff': .7}

marker_labeller = marker_labelling.MarkerLabeller(label_transfer_params=params,
                                                  binary_params=binary_params)

expr = datas[0].to_df().transpose()
broad_labels = marker_labeller.fit_predict(expr, marker_df)
if type(marker_labeller._classifier) != type(None):
    marker_labeller._classifier.printConfusionMatrix()

datas[0].obs['binary_labels'] = broad_labels
datas[0].obs['binary_labels'] = datas[0].obs['binary_labels'].astype('category')

st.pl.cluster_plot(datas[0], dpi=100, use_label='binary_labels')
plt.show()

######## Attempts at measuring similarity between clusters/clustering ##########
# Adding in the species information, subset to data/mix, cluster pseudobulk #
human_datas = []
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                 'C1_untreated_species.txt', 'D1_untreated_species.txt']
dfs = []
labels = []
for i, species_meta in enumerate(species_metas):
    # Retrieving the saved meta data #
    meta = pandas.read_csv(data_dir+'spot_meta/'+species_meta, sep=' ',
                           index_col=0, header=0)
    species = meta.values[:, 0].astype(str)
    human_mix = numpy.logical_or(species=='data', species=='mix')

    # Subsetting to data/mix and checking with plot #
    data_ordered = datas[i][meta.index.values, :]
    data_human_mix = data_ordered[human_mix, :]
    data_human_mix.obs['species'] = species[human_mix]
    sc.pl.spatial(data_human_mix, color='species')
    human_datas.append( data_human_mix )

    # Retrieving df and pseudo-bulking #
    df = data_human_mix.to_df().transpose()
    human_genes = hs.getSpeciesGenes(df, 'hg38-')
    df = df.loc[human_genes, :]
    dfs.append( df )
    # The sample/cluster labels #
    sample = species_meta.split('_')[0]
    labels.extend( [f'{sample}_{cluster}'
                    for cluster in data_human_mix.obs['louvain'].values] )

labels = numpy.array(labels)

# Concatenating the datas & pseudobulking #
data_df = pandas.concat(dfs, axis=1)
nan_bool = numpy.apply_along_axis(numpy.isnan, 1, data_df.values)
nan_bool = numpy.apply_along_axis(numpy.any, 1, nan_bool)
data_df = data_df.loc[nan_bool==False, :]

label_set = []
samples = ['A1', 'B1', 'C1', 'D1']
for sample in samples:
    sample_labels = labels[[sample in label for label in labels]]
    sample_labels = numpy.unique(sample_labels)
    clusters = numpy.array([int(label.split('_')[1]) for label in sample_labels])
    clusters_ordered = clusters[numpy.argsort(clusters)]
    label_set.extend([sample+'_'+str(cluster) for cluster in clusters_ordered])

data_df_med = heat_hs.formatData(data_df, labels, label_set, 'median')[0]

# Removing any genes which have the consistent minimum across clusters #
# for sample in samples:
#     sample_bool = [sample in cluster for cluster in data_df_med.columns]
#     sample_expr = data_df_med.values[:,sample_bool]
#     min_vals = numpy.apply_along_axis(min, 1, sample_expr)
#     min_bool =

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(distance_threshold=0, n_clusters=None,
                                affinity='cosine', linkage='average')
model = model.fit(data_df_med.values.transpose())

plt.title('Hierarchical Clustering Dendrogram')
# plot the top three levels of the dendrogram
plot_dendrogram(model, truncate_mode='level', #p=0,
                labels=data_df_med.columns.values)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.show()

# TODO color dendrogram based on colors of clusters
# TODO add clusters to annData and visualise the results to evaluate if
#       clusters look consistent across samples.

# Getting pair-wise cosine similarity #
dists = pairwise_distances(data_df_med.values.transpose(), metric='cosine')
sims = 1-dists
sims = pandas.DataFrame(sims,
                        index=data_df_med.columns, columns=data_df_med.columns)








