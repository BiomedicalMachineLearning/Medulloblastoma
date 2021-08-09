""" Helper functions for stlearn_clustering.py predominantly.
"""

import os
import scanpy as sc
import stlearn as st
import numpy, pandas
import numpy as np
import pandas as pd

from anndata import AnnData

file_dir = os.path.realpath(__file__)
project_dir = file_dir.split('scripts/')[0]
data_dir = project_dir+'data/'

def species_subset(data, species_meta, subset_labels):
    """Subsets to the species"""
    # Retrieving the saved meta data #
    meta = pandas.read_csv(data_dir+'spot_meta/'+species_meta, sep='\t', #NOTE was sep=' ' for old species classifications.
                           index_col=0, header=0)
    species = meta.values[:, 0].astype(str)
    subset_labels = subset_labels.split('_')
    subset_bool = species == subset_labels[0]
    for label in subset_labels[1:]:
        subset_bool = numpy.logical_or(subset_bool, species==label)

    # Subsetting to data/mix and checking with plot #
    data_ordered = data[meta.index.values, :]
    data_human_mix = data_ordered[subset_bool, :]
    data_human_mix.obs['species'] = species[subset_bool]
    sc.pl.spatial(data_human_mix, color='species')

    return data_human_mix

def load(sample_dir, spot_filter, species_meta, subset=True,
         subset_labels='human_mix'):
    """Loads in the data and filters for QC and data/mix spots"""
    data = st.Read10X(sample_dir)[spot_filter, :]
    data.var_names_make_unique()
    data.var_names = numpy.array([var_name.replace('_', '-')
                                  for var_name in data.var_names])
    # corrects errors down-stream #
    # data.uns["spatial"] = data.uns["spatial"]["Visium8_A1_Hybrid"]

    ### Subsetting to data/mix spots ###
    if subset:
        data = species_subset(data, species_meta, subset_labels)

    return data

def load_and_norm(sample_dir, spot_filter, species_meta,
                  min_cells=3, n_pca=50, n_neighbors=25, scale=False,
                  sme=True, subset=True, subset_labels='human_mix'):
    """ Loads in the 10X visium data, filters, and normalises.

    Args:
        data_dir (str): Directory containing the visium data.
        spot_filter :
        min_cells (int): min no. cells gene needs to be expressed in.

    Returns:
        AnnData: STSME normalised data.
    """
    data = load(sample_dir, spot_filter, species_meta,
                subset=subset, subset_labels=subset_labels)

    st.pp.filter_genes(data, min_cells=min_cells)
    st.pp.normalize_total(data)
    st.pp.log1p(data)

    #### Adding morphology #####
    if sme:
        st.pp.tiling(data, out_path="./temp_tiling", crop_size=40)
        st.pp.extract_feature(data)
        # apply stSME to normalise log transformed data
        st.em.run_pca(data, n_comps=n_pca)
        st.spatial.SME.SME_normalize(data, use_data="raw")
        data.X = data.obsm['raw_SME_normalized']

    if scale:
        st.pp.scale(data)
    st.em.run_pca(data, n_comps=n_pca)
    st.pp.neighbors(data, n_neighbors=n_neighbors, use_rep='X_pca')

    return data

def just_norm(data,
              min_cells=3, n_pca=50, n_neighbors=25, scale=False, sme=True):
    st.pp.filter_genes(data, min_cells=min_cells)
    st.pp.normalize_total(data)
    st.pp.log1p(data)

    #### Adding morphology #####
    if sme:
        st.pp.tiling(data, out_path="./temp_tiling", crop_size=40)
        st.pp.extract_feature(data)
        # apply stSME to normalise log transformed data
        st.em.run_pca(data, n_comps=n_pca)
        st.spatial.SME.SME_normalize(data, use_data="raw")
        data.X = data.obsm['raw_SME_normalized']

    if scale:
        st.pp.scale(data)
    st.em.run_pca(data, n_comps=n_pca)
    st.pp.neighbors(data, n_neighbors=n_neighbors, use_rep='X_pca')

def subset_to_human_spots(datas):
    """ Subsets the anndatas to data spots.

    :return:
    """
    human_datas = []
    species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']
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

    return human_datas

def norm_df(df, scale=False, log1p=True):
    """Normalises data using scanpy."""
    ref_ad = sc.AnnData(df.values.transpose())
    sc.pp.normalize_total(ref_ad)
    if log1p:
        sc.pp.log1p(ref_ad)
    if scale:
        st.pp.scale(ref_ad)
    ref_norm = pandas.DataFrame(ref_ad.X.transpose(),
                                index=df.index, columns=df.columns)

    return ref_norm

def species_split(data, species='human', remove_genes=False,
                  out_format=None):
    """ Splits anndata into data/mix spots for data genes.
        remove_genes indicates whether to remove mitochondrial/ribosomal genes.
        out_format='human' or 'mouse' to format genes either way, default is same
        as input species.
    """
    species_labels = data.obs['species'].values.astype(str)
    mix_bool = species_labels == 'mix'
    human_bool = species_labels == 'human'
    mouse_bool = species_labels == 'mouse'

    h_prefix, m_prefix = 'hg38-', 'mm10-'
    prefix = h_prefix if 'human' in species else m_prefix
    if species == 'human':
        species_bool = np.logical_or(human_bool, mix_bool)
    elif species == 'mouse':
        species_bool = np.logical_or(mouse_bool, mix_bool)
    elif species == 'human-only':
        species_bool = human_bool
    elif species == 'mouse-only':
        species_bool = mouse_bool
    elif species == 'mix-only':
        species_bool = mix_bool
        prefix = ''
    elif species == 'mix-humang':
        species_bool = mix_bool
    elif species == 'mix-mouseg':
        species_bool = mix_bool

    # NOT implimented.... could be.
    if remove_genes and 'human' in species:
        remove_prefixes = ['MT-', 'RPL', 'RPS']
    elif remove_genes and 'mouse' in species:
        remove_prefixes = ['mt-', 'Rpl', 'Rps']
    elif remove_genes and 'mix-only':
        remove_prefixes = ['MT-', 'RPL', 'RPS', 'mt-', 'Rpl', 'Rps']
    else:
        remove_prefixes = []

    data_species = data[species_bool, :]

    species_genes = [gene for gene in data.var_names if prefix in gene]
    species_df = data_species[:,species_genes].to_df()
    split_cond = prefix if prefix!='' else ' '
    species_genes_formatted = np.array([gene.split(split_cond)[-1]
                                        for gene in species_genes])
    species_df.columns = species_genes_formatted

    if out_format=='human' and prefix!='hg38-':
        genes = [gene.upper() for gene in species_df.columns]
        species_df.columns = genes
    elif out_format=='mouse' and prefix!='mm10-':
        genes = [gene[0].upper()+gene[1:].lower() for gene in species_df.columns]
        species_df.columns = genes

    species_ad = AnnData(species_df, obs=data_species.obs, uns=data_species.uns,
                                 obsm=data_species.obsm, obsp=data_species.obsp)

    return species_ad



