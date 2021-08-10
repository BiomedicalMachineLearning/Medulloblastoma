"""
    Loads in the fetal brain data, & UMAP from seurat
    to display the data in a figure.
    Ref. paper: https://www.nature.com/articles/s41586-020-2157-4
    Download link: https://db.cngb.org/HCL/gallery.html?tissue=Fetal-Brain3

      INPUT: * data/third_party_data/HCL2020/Fetal-Brain3_dge.txt
             * data/third_party_data/HCL2020/Fetal-Brain3_Anno.txt

      OUTPUT: * data/third_party_data/HCL2020/FetalBrain3.h5ad
              * figure_components/HumanAnnot_figures/*
"""

################################################################################
                        # Environment setup #
################################################################################
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad

import scripts.utils.helpers as hs
hs.setUp()

import scripts.utils.visualisation.helpers as vhs

data_dir = 'data/third_party_data/HCL2020/'
out_dir = 'data/third_party_data/HCL2020/'
out_plots = 'figure_components/HumanAnnot_figures/'

################################################################################
                        # Loading in the data #
################################################################################
hcl_counts = pd.read_csv(data_dir+'Fetal-Brain3_dge.txt', index_col=0)
hcl_meta = pd.read_csv(data_dir+'Fetal-Brain3_Anno.csv', index_col=0)

hcl = ad.AnnData(hcl_counts.transpose(), obs=hcl_meta)

################################################################################
                    # Dimensionality reduction #
################################################################################
##### Using saved UMAP from Seurat FetalBrain3_SingleR.R
##### and will load here
hcl_umap = pd.read_csv(data_dir+'Fetal-Brain3_seuratUMAP.txt',
                       index_col=0, sep='\t')
hcl.obsm['X_umap'] = hcl_umap.values

sc.pl.umap(hcl, color='CT')

################################################################################
               # Saving UMAP for visualisation in figure #
################################################################################
fig, ax = plt.subplots(figsize=(12, 8))
sc.pl.umap(hcl, color='CT', ax=ax, show=False)
vhs.dealWithPlot(True, True, True,
                 out_plots, 'FetalBrain3_celltype_umap.pdf', 300)
