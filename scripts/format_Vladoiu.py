""" The purpose of this script is to format the Vladoiu developing Cerebellum
    scRNA-seq data for fast loading.

    INPUT: data/third_part_data/Vladoiu2019_Nature_scRNA/
                                GSE118068_RAW/*

    OUTPUT: data/third_part_data/Vladoiu2019_Nature_scRNA/
                                counts.hdf & cell_meta.hdf.
"""

import os, sys
import numpy, pandas

import beautifulcells.preprocessing.load_data.load_10X as load_10X

################################################################################
                        # Reading in the data #
################################################################################
data_dir = '/30days/uqbbalde/MedullaBlastoma/'
raw_dir = data_dir

# Getting file prefixes and grouping into file-type (matrix, gene, barcode) #
prefixes = numpy.unique(['_'.join(file_name.split('_')[0:2])
                        for file_name in os.listdir(raw_dir)
                         if not file_name.startswith('.')
			and file_name.endswith('.gz')])
suffixes = ['_matrix.mtx.gz', '_genes.tsv.gz', '_barcodes.tsv.gz']

# Loading in the count matrices as dataframes #
all_counts = None
time_labels = []
for i, prefix in enumerate(prefixes):
    # Loading in the 10X data #
    paths = [raw_dir+prefix+suffix for suffix in suffixes]
    count_df = load_10X.read_10X(paths[0], paths[1], paths[2])
    n_cells = count_df.shape[1]
    time = prefix.split('_')[1] if 'E18' not in prefix else 'E18.5'
    count_df.columns = [f'{time}_{name.replace("-1", "")}' for name in count_df.columns]
    if type(all_counts) == type(None):
        all_counts = count_df
    else:
        all_counts = pandas.concat([all_counts, count_df], axis=1)
        del count_df

    print(all_counts.shape)
    # Saving the time point labels #
    time_labels += [time] * n_cells

# Formatting the cell meta data #
cell_meta = pandas.read_csv(data_dir+'cluster_annotations.csv',
                                      sep=',', index_col=0)

#subsetted to just annotated
all_counts = all_counts.loc[:, cell_meta.index]

cell_meta['time'] = [name.split('_')[0] for name in cell_meta.index]

# Converting the cluster labels to cell types #
cluster_to_label = {0: 'Excitatory cerebellar nuclei neurons',
		   1: 'Embryonic and postnatal GCPs-1',
		   2: 'Neural stem cells', 3: 'Unipolar brush cell and GCP progenitor',
		   4: 'Unipolar brush cells', 5: 'GABA interneurons',
		   6: 'Brainstem progenitors', 7: 'Granule cells',
		   8: 'VZ progenitors', 9: 'Unipolar brush cell precursors',
		   9: 'Unipolar brush cell precursors', 
		   10: 'Differentiating Purkinje cells',
		   11: 'Gliogenic progenitors-1', 
		   12: 'Upper rhombic lip progenitors',
		   13: 'Mesenchymal stem cells-1',
		   14: 'Purkinje cells', 15: 'Postnatal GCPs-2',
		   16: 'Post mitotic NTZ neurons', 
		   17: 'Roof plate-like stem cells',
		   18: 'Proliferating VZ progenitors',
		   19: 'Oligodendrocyte precursor cells', 
		   20: 'Gliogenic progenitors-2',
		   21: 'Astrocyte/Bergmann glia precursors',
		   22: 'Endothelial cells', 
		   23: 'Postnatal excitatory cerebellar nuclei neurons',
		   24: 'GABA interneuron precursors',
		   25: 'Pericytes', 26: 'Early proliferating VZ progenitors',
		   27: 'Mesenchymal stem cells-2', 
		   28: 'Microglia', 29: 'Meninges', 30: 'Red blood cells'}

labels = [cluster_to_label[cluster] for cluster in cell_meta['Cluster identity']]

cell_meta['cell_labels'] = labels

###### Removing any genes with litle expression ######
new_counts = pandas.DataFrame(all_counts.values,
                              index=all_counts.index, columns=all_counts.columns)

def getFreq(x):
	return len(numpy.where(x>0)[0])

gene_freqs = numpy.apply_along_axis(getFreq, 1, new_counts.values)
gene_bool = gene_freqs>50
print(sum(gene_bool))

filt_counts = new_counts.loc[gene_bool, :]

# For some reason, there are cases where a gene has a duplicate name.
# Going to filter these rare cases out, and double check downstream that
# the gene expressions match up to those expected from the paper.
index_freqs = [len(numpy.where(filt_counts.index==gene)[0]) for gene in filt_counts.index]
genes = filt_counts.index.values[numpy.array(index_freqs)==1]

filt_counts = filt_counts.loc[genes, :]

###### Saving the output ######
data_out = '/30days/uqbbalde/MedullaBlastoma/Vladoiu2019_Nature_scRNA/'
cell_meta.to_csv(data_out+'Vladoiu2019_cell_meta.txt', sep='\t')

filt_counts.to_hdf(data_out+'Vladoiu2019_cell_counts.hdf', 'counts')






















