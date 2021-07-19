""" Helper functions for stlearn_clustering.py predominantly.
"""

import scanpy as sc
import stlearn as st
import numpy, pandas
import numpy as np
import pandas as pd

data_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/data/'

def species_subset(data, species_meta):
    """Subsets to the species"""
    # Retrieving the saved meta data #
    meta = pandas.read_csv(data_dir+'spot_meta/'+species_meta, sep=' ',
                           index_col=0, header=0)
    species = meta.values[:, 0].astype(str)
    human_mix = numpy.logical_or(species=='data', species=='mix')

    # Subsetting to data/mix and checking with plot #
    data_ordered = data[meta.index.values, :]
    data_human_mix = data_ordered[human_mix, :]
    data_human_mix.obs['species'] = species[human_mix]
    sc.pl.spatial(data_human_mix, color='species')

    return data_human_mix

def load(sample_dir, spot_filter, species_meta, subset=True):
    """Loads in the data and filters for QC and data/mix spots"""
    data = st.Read10X(sample_dir)[spot_filter, :]
    data.var_names_make_unique()
    data.var_names = numpy.array([var_name.replace('_', '-')
                                  for var_name in data.var_names])
    # corrects errors down-stream #
    # data.uns["spatial"] = data.uns["spatial"]["Visium8_A1_Hybrid"]

    ### Subsetting to data/mix spots ###
    if subset:
        data = species_subset(data, species_meta)

    return data

def load_and_norm(sample_dir, spot_filter, species_meta,
                  min_cells=3, n_pca=50, n_neighbors=25, scale=False,
                  sme=True, subset=True):
    """ Loads in the 10X visium data, filters, and normalises.

    Args:
        data_dir (str): Directory containing the visium data.
        spot_filter :
        min_cells (int): min no. cells gene needs to be expressed in.

    Returns:
        AnnData: STSME normalised data.
    """
    data = load(sample_dir, spot_filter, species_meta, subset=subset)

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
    st.pp.neighbors(data,
                    n_neighbors=n_neighbors, use_rep='X_pca')

    return data

def load_and_norm_LRAnalysis(sample_dir, spot_filter, species_meta,
                             min_cell_prop=0.025, subset=True):
    """Loads & normalises the data in appropriate way for LR analysis."""
    data = load(sample_dir, spot_filter, species_meta, subset=subset)
    st.pp.filter_genes(data, min_cells=int(min_cell_prop * data.n_obs))
    st.pp.normalize_total(data)

    # Adding in tissue_type info #
    # sample = species_meta.split('_')[0]
    # meta = pd.read_csv(data_dir + 'spot_meta/clusters_human_mix.txt',
    #                    sep='\t', index_col=0)
    # obs_names = np.array(
    #     [name for name in data.obs_names if f'{sample}-{name}' in meta.index])
    # data = data[obs_names, :]
    # obs_names = np.array([f'{sample}-{name}' for name in data.obs_names])
    # tissue_type = meta.loc[obs_names, 'tissue_type'].astype('category')
    # tissue_type.index = data.obs_names
    # data.obs['tissue_type'] = tissue_type
    return data

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

def norm_df(df, scale=False):
    """Normalises data using scanpy."""
    ref_ad = sc.AnnData(df.values.transpose())
    sc.pp.normalize_total(ref_ad)
    sc.pp.log1p(ref_ad)
    if scale:
        st.pp.scale(ref_ad)
    ref_norm = pandas.DataFrame(ref_ad.X.transpose(),
                                index=df.index, columns=df.columns)

    return ref_norm

