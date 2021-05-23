"""
Performs the QC looks at the distribution of the counts/genes & also the
ribosomal/mt gene captured.

OUTPUT:

        * data/spot_meta/QC/{sample_name}_QCd_barcodes.txt ->
                              Contains the barcodes of the spots that passed QC.

        * figure_components/QC_out/ -> Contains the plots from the QC.
"""

################################################################################
                    # Environment setup #
################################################################################
import os, sys
import stlearn as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scripts.May_REDO.helpers as hs
import matplotlib
import beautifulcells.visualisation.helpers as vhs
import beautifulcells.preprocessing.QC.cell_QC as cell_QC

out_plots = 'figure_components/QC_out/'

hs.setUp()

samples_ = ['Visium8_A1_Hybrid_treated/', 'Visium8_B1_Hybrid_treated/',
            'Visium8_C1_Hybrid_untreated/', 'Visium8_D1_Hybrid_untreated/']
samples_2 = ['A/', 'B/', 'C/', 'D/']

################################################################################
                    # Loading in the data #
################################################################################
prefixes = ['A1-', 'B1-', 'C1-', 'D1-']

# Filtering the spots #
species_metas = ['A1_treated_species.txt', 'B1_treated_species.txt',
                     'C1_untreated_species.txt', 'D1_untreated_species.txt']
filt_files = np.array( os.listdir(hs.data_dir+'ryan_ids/') )
filt_files = filt_files[np.argsort(filt_files)]
spot_filters = [pd.read_csv(hs.data_dir+'ryan_ids/'+filt_file,
                           header=None).values[:,0] for filt_file in filt_files]
for i in range(len(spot_filters)):
    spot_filters[i] = [prefixes[i]+bc for bc in spot_filters[i]]

datas = [hs.load(hs.data_dir+samples_[i], prefixes[i])[spot_filters[i],:] for i in
         range(len(samples_))]
# Loading the data from the second model #
#datas.extend( [hs.load(hs.data_dir2+samples_2[i]) for i in range(len(samples_2))] )
# No counts in the tumour areas...

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
ribo_gene_perc = [hs.getGenePercents(data, ribo, n_expr=True)
                                                              for data in datas]
ribo_count_perc = [hs.getGenePercents(data, ribo, n_expr=False)
                                                              for data in datas]

mito_gene_perc = [hs.getGenePercents(data, mito, n_expr=True)
                                                              for data in datas]
mito_count_perc = [hs.getGenePercents(data, mito, n_expr=False)
                                                              for data in datas]

# for i in range(len(datas)):
#     datas[i].obsm['ribo_gene_perc'] = ribo_gene_perc[i]
#     st.pl.het_plot(datas[i], use_het='ribo_gene_perc')
#     plt.show()
#
#     datas[i].obsm['mito_gene_perc'] = mito_gene_perc[i]
#     st.pl.het_plot(datas[i], use_het='mito_gene_perc')
#     plt.show()

# Removing all of the ribosomal/mitochondrial genes before filtering further #
#datas = [data[:,good_genes] for data in datas]

# Removing genes with extremely low counts #
#[st.pp.filter_genes(data, min_cells=4) for data in datas]

overlap_matrix = np.zeros((len(datas), len(datas)))
human_matrix = np.zeros((len(datas), len(datas)))
mouse_matrix = np.zeros((len(datas), len(datas)))
for i, data in enumerate( datas ):
    human_genes = [gene for gene in data.var_names if 'hg38-' in gene]
    mouse_genes = [gene for gene in data.var_names if 'mm10-' in gene]
    print(samples_[i])
    print(f"Total genes detected: {len(data.var_names)}")
    print(f"Human genes detected: {len(human_genes)}")
    print(f"Mouse genes detected: {len(mouse_genes)}")
    for j, data2 in enumerate( datas ):
        overlap_matrix[i,j] = len(set(data.var_names).intersection(set(data2.var_names)))
        human_genes2 = [gene for gene in data2.var_names if 'hg38-' in gene]
        mouse_genes2 = [gene for gene in data2.var_names if 'mm10-' in gene]
        human_matrix[i, j] = len(set(human_genes).intersection(set(human_genes2)))
        mouse_matrix[i, j] = len(set(mouse_genes).intersection(set(mouse_genes2)))
# Extremely high gene recovery overlap between samples, which is great !
"""
Visium8_A1_Hybrid_treated/
Total genes detected: 64591
Human genes detected: 33538
Mouse genes detected: 31053

Visium8_B1_Hybrid_treated/
Total genes detected: 64591
Human genes detected: 33538
Mouse genes detected: 31053

Visium8_C1_Hybrid_untreated/
Total genes detected: 64591
Human genes detected: 33538
Mouse genes detected: 31053

Visium8_D1_Hybrid_untreated/
Total genes detected: 64591
Human genes detected: 33538
Mouse genes detected: 31053
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






