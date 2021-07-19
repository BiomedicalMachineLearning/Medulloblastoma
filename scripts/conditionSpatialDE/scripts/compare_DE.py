""" The purpose of this script is to read-in and compare the results from the
    DE analysis from various different analyses. Assumes that stored in excel
    sheets, with the equivalent gene sets have the same gene names and the
    first column indicating the gene names, and a column indicating the padj.
"""

import os, sys
import numpy, pandas
from collections import defaultdict
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

import scripts.conditionSpatialDE.scripts.helpers as hs
hs.setUp()

import gseapy

import beautifulcells.postprocessing.format.format as format
import beautifulcells.visualisation.helpers as vhs

################################################################################
                        # Reading in the data #
################################################################################
gene_de_names = ['data.treatment_up', 'data.treatment_down']
de_files = ['data/DE_out/Pseudo_TMM_Limma_Voom/de_results_PseudoLimma_withoutLogFC.xlsx',
            'data/supps/de_results_MDB.xlsx']
de_names = ['PseudoLimma', 'DESeq2']
de_dfs = [pandas.read_excel(de_file, gene_de_names,
                        engine='openpyxl', index_col=0) for de_file in de_files]

# Getting pre-ranked lists to perform GSEA
# going from down to up-regulated, ranked by logFC #
de_species = numpy.unique([name.split('.')[0] for name in gene_de_names])
for de_specie in de_species:
    for gene_de_name in gene_de_names:
        for species in ['hg38-', 'mm10-']:
            gene_sets = []
            for de_df in de_dfs:
                species_de = hs.getSpeciesDERanked(de_df, [gene_de_name])[species]
                gene_sets.append(set(species_de))

            venn2(gene_sets, set_labels=[f'{name}_{species}' for name in de_names])
            vhs.dealWithPlot(True, True, False, 'data/DE_out/',
                             f'DE_compare_venn_{de_specie}_{species}.pdf', 300)

















