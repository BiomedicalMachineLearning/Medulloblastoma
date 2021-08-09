"""
Performs the QC to filter spots & looks at the distribution of the counts/genes.

             INPUT: data/Visium8_{sample_name}_Hybrid/*
             OUTPUT: * data/filter_ids/{sample_name}_filtered.txt ->
                              Contains the barcodes of the spots that passed QC.
                     * figure_components/QC_figure/*
"""

################################################################################
                    # Environment setup #
################################################################################
import os, sys
import stlearn as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scripts.X1_QC_SpeciesClassify.QC_helpers as QC_hs
import scripts.utils.helpers as hs
import matplotlib
import scripts.utils.visualisation.helpers as vhs
import scripts.utils.preprocessing.QC.cell_QC as cell_QC

out_plots = 'figure_components/QC_figure/'
hs.setUp()

samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
samples_2 = ['A/', 'B/', 'C/', 'D/']

################################################################################
                    # Loading in the data #
################################################################################
prefixes = ['A1-', 'B1-', 'C1-', 'D1-']
datas = [hs.load(hs.data_dir+samples_[i], prefixes[i]) for i in
         range(len(samples_))]

################################################################################
                # Performing QC spot filtering #
################################################################################
dfs = [data.to_df() for data in datas]
df_all = pd.concat(dfs).transpose()
df_all.values[np.isnan(df_all.values)] = 0 #Removing nans

count_gene_ratios, umis, gene_counts = cell_QC.getQCMetrics(df_all)

matplotlib.rcParams.update({'font.size': 8, 'font.weight': 'bold'})
cutoffs = [0, 0, 200]
sf = f'MBVisium_'
passed = cell_QC.plotQCHistograms(count_gene_ratios, umis, gene_counts, cutoffs,
                                  folder=out_plots, savePlot=True,
                                  suffix=sf,
                                  linewidth=2, figsize=(5, 3), showPlot=True
                                  )

stats = cell_QC.getQCStats(count_gene_ratios, umis, gene_counts, cutoffs,
                          saveFile=True, folder='data/spot_meta/QC/', suffix=sf)

################################################################################
                    # Performing spot QC #
################################################################################
# Identifying ribosomal & mitochondrial genes #
genes = datas[0].var_names
ribo = [gene for gene in genes if \
        '-RPL' in gene or '-RPS' in gene or \
        '-Rpl' in gene or '-Rps' in gene]
mito = [gene for gene in genes if '-MT-' in gene or '-mt-' in gene]
good_genes = [gene for gene in genes if gene not in ribo and gene not in mito]

# Getting percentage of ribo/mito genes expressed per spot #
ribo_gene_perc = [QC_hs.getGenePercents(data, ribo, n_expr=True)
                                                              for data in datas]
ribo_count_perc = [QC_hs.getGenePercents(data, ribo, n_expr=False)
                                                              for data in datas]

mito_gene_perc = [QC_hs.getGenePercents(data, mito, n_expr=True)
                                                              for data in datas]
mito_count_perc = [QC_hs.getGenePercents(data, mito, n_expr=False)
                                                              for data in datas]

overlap_matrix = np.zeros((len(datas), len(datas)))
for i, data in enumerate( datas ):
    print(samples_[i])
    print(f"Total genes detected: {len(data.var_names)}")
    for j, data2 in enumerate( datas ):
        overlap_matrix[i,j] = len(set(data.var_names).intersection(set(data2.var_names)))

# Extremely high gene recovery overlap between samples, which is great !
"""
Visium8_A1_Hybrid_treated/
Total genes detected: 64591

Visium8_B1_Hybrid_treated/
Total genes detected: 64591

Visium8_C1_Hybrid_untreated/
Total genes detected: 64591

Visium8_D1_Hybrid_untreated/
Total genes detected: 64591
"""

################################################################################
            # Visualising total counts after filtering genes #
################################################################################
mins, maxs = [], []
sample_counts = []
for i, data in enumerate(datas):
    df = data.to_df()
    tot_counts = df.values.sum(axis=1)
    min_counts = min(tot_counts)
    max_counts = max(tot_counts)
    mins.append( min_counts )
    maxs.append( max_counts )

    data.obsm['total_counts'] = tot_counts
    sample_counts.append(tot_counts)

min_counts = 0#min(mins)
max_counts = 8000#np.min(maxs)

r_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('rcmap',
                                    [(0, 0, 0)]+[(.6, 0, 0)]*1+[(.7, 0, 0)]*1+\
                                    [(.8, 0, 0)]*1+[(.9, 0, 0)]*4+[(.95, 0, 0)])

for i, data in enumerate(datas):
    st.pl.het_plot(data, cmap='Spectral_r',#r_cmap,
                   size=12,
                   use_het='total_counts', vmin=min_counts, vmax=max_counts)
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{samples_[i].split("_")[1]}_total_counts.pdf', dpi=300)

################################################################################
       # Performing QC spot filtering for each sample indepedently #
################################################################################
dfs_t = [df.transpose() for df in dfs]

# To make sure on the same x-axis scale #
mins = [min(vals) for vals in [count_gene_ratios, umis, gene_counts]]
maxs = [max(vals) for vals in [count_gene_ratios, umis, gene_counts]]
valRanges = [(mins[i], maxs[i]) for i in range(len(mins))]

sample_stats = []
for i, df in enumerate(dfs_t):
    prefix = prefixes[i]

    count_gene_ratios, umis, gene_counts = cell_QC.getQCMetrics(df)

    matplotlib.rcParams.update({'font.size': 8, 'font.weight': 'bold'})
    cutoffs = [0, 0, 200]
    sf = f'MBVisium_{prefix}'
    passed = cell_QC.plotQCHistograms(count_gene_ratios, umis, gene_counts,
                                      cutoffs,
                                      folder=out_plots, savePlot=True,suffix=sf,
                                      linewidth=2, figsize=(5, 3), showPlot=True,
                                      ratioValRange=valRanges[0],
                                      umiValRange=valRanges[1],
                                      geneValRange=valRanges[2]
                                      )

    stats = cell_QC.getQCStats(count_gene_ratios, umis, gene_counts, cutoffs,
                               saveFile=True, folder='data/spot_meta/QC/',
                               suffix=sf)
    sample_stats.append( stats )

cols = ['umis_median', 'umis_mean', 'gene_counts_median', 'gene_counts_mean',
        'umi_gene_ratio_median', 'umi_gene_ratio_mean']
for i in range(len(sample_stats)):
    print(prefixes[i])
    print(sample_stats[i].loc['after',cols])
    print('\n')

################################################################################
                # Visualising total counts after QCd #
################################################################################
passed_bcs = df_all.columns.values[passed]
datas_qcd = []
for data in datas:
    data_bcs = [bc for bc in passed_bcs if bc in data.obs_names]
    datas_qcd.append( data[data_bcs, :] )

for i, data in enumerate(datas_qcd):
    st.pl.het_plot(data, cmap='Spectral_r',#r_cmap,
                   size=12,
                   use_het='total_counts', vmin=min_counts, vmax=max_counts)
    vhs.dealWithPlot(True, True, True, out_plots,
                     f'{samples_[i].split("_")[1]}_total_counts_QCd.pdf',
                     dpi=300)






