""" Loading in the data from GSEA_analyse to visualise by clustering the
    genes found in each gene set & their logFC scores.

    OUTPUT:
        figure_components/PseudoLimma/DE_figures/*
"""

################################################################################
                    # Environment setup #
################################################################################
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scripts.May_REDO.helpers as hs
import scripts.May_REDO.DE_figure_helpers as dhs
from sklearn.preprocessing import scale

out_plots = 'figure_components/PseudoLimma/DE_figures/'

hs.setUp()

################################################################################
                    # Loading in the data #
################################################################################
de_dir = 'data/DE_out/Pseudo_TMM_Limma_Voom/data/'
gsea_dir = 'data/DE_out/Pseudo_TMM_Limma_Voom/gsea_out/interesting_results/'
gene_stats = pd.read_csv(de_dir+'gene_stats.txt', sep='\t')
gsea_results = pd.read_excel(gsea_dir+'gsea_summary_human.treatment.xlsx',
                           engine='openpyxl')
gsea_results.index = gsea_results.loc[:,'Term']

# Getting the total number of unique genes #
interesting_terms = np.array(['Epithelial Mesenchymal Transition',
                     'Cell adhesion molecules (CAMs)', 'Myogenesis',
                      'cGMP-PKG signaling pathway', 'KRAS.600 UP.V1 UP',
                      'cAMP signaling pathway', 'Coagulation',
                      'G1 to S cell cycle control WP45', 'ATM human',
#'ATM Signaling Network in Development and Disease  WP3878',
#'DNA IR-Double Strand Breaks (DSBs) and cellular response via ATM WP3959',
'E2F3 human', 'Integrated Cancer Pathway WP1971', 'E2F4 human',
'GCNP SHH UP EARLY.V1 UP', 'Myc Targets V1', 'GCNP SHH UP LATE.V1 UP',
'PRC2 EZH2 UP.V1 UP', 'E2F1 UP.V1 UP', 'Cell cycle', 'E2F1 human',
'VEGF A UP.V1 DN', 'E2F Targets', 'G2-M Checkpoint',
    ])
nes_values = gsea_results.loc[interesting_terms,'nes'].values
interesting_terms = interesting_terms[np.argsort(-nes_values)]
terms = gsea_results.loc[:,'Term'].values
term_indices = [np.where(terms==term)[0][0] for term in interesting_terms]
ledge_genes = gsea_results.loc[:,'ledge_genes'].values[term_indices]
genes = []
terms_to_genes = {}
for i, gene_str in enumerate(ledge_genes):
    gene_ = gene_str.split(';')
    genes.extend( gene_ )

    terms_to_genes[interesting_terms[i]] = gene_
genes = np.unique(genes)
# Only 288 of the genes found in gene sets !

# Creating a matrix to visualise the results #
gsea_matrix = np.zeros((len(interesting_terms), len(genes)), dtype=np.float64)
for i, term in enumerate(interesting_terms):
    for gene in terms_to_genes[term]:
        j = np.where(genes==gene)[0][0]
        if 'ORF' in gene:
            gene = gene.replace('ORF', 'orf')

        gsea_matrix[i,j] = gene_stats.loc[f'hg38-{gene}','logFC']
gsea_df = pd.DataFrame(gsea_matrix, index=interesting_terms, columns=genes)

################################################################################
                    # Creating the Visualisation #
################################################################################
dhs.plot_gseaClusterMap(gsea_df, out_plots+'gsea_map.pdf')
dhs.plot_gseaDotPlot(interesting_terms, gsea_results,
                     out_plots+'gsea_dotplot.pdf', )









