""" Created: 9th April 2021

Uses the stlearn clusters, pseudobulks, normalises, & then perform TRIAGE
analysis on these to identify biologically important genes for each cluster to
help justify annotations.

OUTPUT: data/DE_out/stlearn_cluster_DE/triage_by_cluster.xlsx
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

import scripts.helpers as hs
hs.setUp()
import scripts.stlearn_helpers as st_hs

import beautifulcells.tools.triage.triage_analysis as triage_analysis
import beautifulcells.postprocessing.format.format as format

################################################################################
                    # Loading data & gene filtering #
################################################################################
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

# Loading in the un-normalised counts, filtered to human/mix spots.
datas = [st_hs.load(sample_dir, spot_filters[i], species_metas[i])
         for i, sample_dir in enumerate(data_dirs)]
# Removing lowly expressed genes
[st.pp.filter_genes(data, min_cells=5) for data in datas]

################################################################################
                    # Pseudobulk & normalise clusters #
################################################################################
# Loading the cluster labels #
meta = pandas.read_csv(data_dir+'spot_meta/clusters_human_mix.txt',
                                                          sep='\t', index_col=0)
tissue = meta.loc[:,'tissue_type'].values

sheet_names = []
triages = []
for i, data in enumerate(datas):
    sample = samples[i]
    sample_indices = [i for i in range(meta.shape[0])
                      if sample in meta.index.values[i]]
    ordered_obs = meta.index.values[sample_indices]
    obs_names = numpy.array([f'{sample}-{name}' for name in data.obs_names])
    counts = data.to_df().transpose()
    counts.columns = obs_names
    counts = counts.loc[:,ordered_obs]
    print(numpy.all(counts.columns.values==ordered_obs))

    # Subsetting to human genes #
    human_genes = hs.getSpeciesGenes(counts, 'hg38-')
    counts = counts.loc[human_genes, :]
    counts.index = hs.getSpeciesGenes(counts, 'hg38-', remove_prefix=True)

    # Creating the pseudobulk #
    sample_tissues = tissue[sample_indices]
    tissue_types = numpy.unique(sample_tissues)
    tissue_types = tissue_types[numpy.argsort(tissue_types)]
    pseudobulk = numpy.zeros( (counts.shape[0], len(tissue_types)) )
    for j, tissue_type in enumerate(tissue_types):
        pseudobulk[:,j] = numpy.apply_along_axis(sum, 1,
                                   counts.values[:,sample_tissues==tissue_type])
    pseudobulk = pandas.DataFrame(pseudobulk,
                                  index=counts.index, columns=tissue_types)

    # Normalising & Scaling #
    pseudobulk = st_hs.norm_df(pseudobulk, scale=True)

    # Performing TRIAGE analysis #
    triage = triage_analysis.Triage(organism='human', average=False,
                                                                 binary_data='')
    triaged = triage.fit_transform(pseudobulk)
    triaged_top = triaged.iloc[0:300,:]
    triaged_top.columns = [f'{name}_TOP' for name in triaged_top.columns]
    triaged_bottom = triaged.iloc[-300::,:]
    triaged_bottom.columns = [f'{name}_BOTTOM' for name in triaged_bottom.columns]
    triaged_bottom = triaged_bottom.iloc[::-1,:]
    triaged_bottom.index = triaged_top.index
    new_triaged = pandas.concat((triaged_top, triaged_bottom), axis=1)

    sheet_names.append(f'{sample}_TRIAGED')
    triages.append( new_triaged )

format.writeDFsToExcelSheets(data_dir+
                             'DE_out/stlearn_cluster_DE/triage_by_cluster.xlsx',
                             triages, sheet_names)











