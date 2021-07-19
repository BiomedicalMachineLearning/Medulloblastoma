""" The purpose of this script is to perform broad GO term analysis
    using gseapy in order to figure out broad changes between treatment
    and control.
"""

import os, sys
import numpy, pandas

import scripts.helpers as hs
hs.setUp()

import gseapy

import beautifulcells.postprocessing.format.format as format

################################################################################
                        # Reading in the data #
################################################################################
gene_de_names = ['mouse.treatment_up', 'mouse.treatment_down',
                 'data.treatment_up', 'data.treatment_down',
                 'mix.treatment_up', 'mix.treatment_down']
de_results = pandas.read_excel('data/supps/de_results_MDB.xlsx', gene_de_names,
                             engine='openpyxl', index_col=0)
# Getting pre-ranked lists to perform GSEA
# going from down to up-regulated, ranked by logFC #
gene_treats = numpy.unique([gene_name.split('_')[0]
                            for gene_name in gene_de_names])
preranked = {}
for treat in gene_treats:
    up_results = de_results[f'{treat}_up']
    down_results = de_results[f'{treat}_down']

    # Filtering by species #
    if 'data' in treat or 'mouse' in treat:
        prefix = 'hg38-' if 'data' in treat else 'mm10-'
        up_results = up_results.loc[hs.getSpeciesGenes(up_results, prefix), :]
        down_results = down_results.loc[
                                    hs.getSpeciesGenes(down_results, prefix), :]
    else:
        continue

    # Ranking by log2FoldChange #
    up_results = hs.reOrderDF(up_results, 'log2FoldChange', False)
    down_results = hs.reOrderDF(down_results, 'log2FoldChange', False)
    pre_ranked_i = pandas.concat((up_results, down_results))
    pre_ranked_i = pre_ranked_i.loc[:,['log2FoldChange']]
    pre_ranked_i.index = hs.getSpeciesGenes(pre_ranked_i, prefix, True, True)
    preranked[treat] = pre_ranked_i

"""
For the data spots, just using the data DE genes, for the mouse spots,
just using mouse DE genes. Ignoring mixed for now, complicated & not important.
"""

#### Getting supported library names #####
human_libs = gseapy.get_library_name('Human')
mouse_libs = gseapy.get_library_name('Mouse')

gene_sets = [#'CCLE_Proteomics_2020', 'Cancer_Cell_Line_Encyclopedia',
             #   'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018',
             #   'GO_Biological_Process_2018',
              'MSigDB_Hallmark_2020', 'MSigDB_Oncogenic_Signatures',
              'TRRUST_Transcription_Factors_2019',
            'KEGG_2019_{0}', 'WikiPathways_2019_{0}'
              ]

######### Performing the GSEA #######
# NOTE: I edited the source code of GSEApy to get this to work, namely had to
#       change the url string query.
gsea_results = {}
for treat in ['mouse.treatment', 'data.treatment']:
    df = preranked[treat]
    species = treat.split('.')[0]
    species = species[0].upper()+species[1::]
    for gene_set in gene_sets:
        gene_set = gene_set.format(species)
        pre_res = gseapy.prerank(rnk=df, gene_sets=gene_set,
                     processes=2,
                     permutation_num=100, # reduce number to speed up testing
                     outdir=f'data/gsea_out/{treat}/{gene_set}',
                                 format='png', seed=6)
        gsea_results[f'{treat}/{gene_set}'] = pre_res

######### Determining the significant results, and moving to a separate
########  summary folder
sig_dfs = {'mouse.treatment': [], 'data.treatment': []}
db_gene_set_labels = {'mouse.treatment': [], 'data.treatment': []}
for folder in gsea_results:
    results = gsea_results[folder].results
    treat = folder.split('/')[0]
    db = folder.split('/')[1]
    results_df = pandas.read_csv(f'data/gsea_out/{treat}/{db}/'
                                 f'gseapy.prerank.gene_sets.report.csv',
                                 index_col=0)
    sig_gene_sets = [] # Storing which gene sets significant
    for gene_set in results:
        fdr = results[gene_set]['fdr']
        if fdr < .05:
            sig_gene_sets.append(gene_set)
            set_name = gene_set.replace(" ", "\ ").\
                        replace("(", "\(").replace(")", "\)")
            out_folder = f'data/gsea_out/interesting_results/{treat}/{db}/'
            os.system(f'mkdir {out_folder}')
            os.system(f'cp data/gsea_out/{folder}/{set_name}.prerank.png '
                      f'{out_folder}')

    db_gene_set_labels[treat] = db_gene_set_labels[treat] + \
                                [db]*len(sig_gene_sets)
    sig_df = results_df.loc[sig_gene_sets, :]
    sig_dfs[treat].append(sig_df)

sig_dfs_treat = []
for treat in sig_dfs:
    sig_df_final = pandas.concat(sig_dfs[treat])
    db_df = pandas.DataFrame(db_gene_set_labels[treat],
                             columns=['database'], index=sig_df_final.index)
    sig_df_final = pandas.concat([db_df, sig_df_final], axis=1)

    sig_df_final = hs.reOrderDF(sig_df_final, 'nes', False)
    sig_dfs_treat.append( sig_df_final)

sheet_names = list(sig_dfs.keys())

format.writeDFsToExcelSheets('data/gsea_out/interesting_results/'
                             'gsea_summary.xlsx', sig_dfs_treat, sheet_names)






