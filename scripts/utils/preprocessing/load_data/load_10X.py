""" Contains functions for reading 10X data from CellRanger.
"""

import scipy.io
import gzip
import csv

import pandas

def getReader(file_path):
    """Gets an appropriate file reader, whether file gzipped or normal text."""
    if '.gz' in file_path:
        reader = csv.reader(gzip.open(file_path, mode='rt'), delimiter="\t")
    else:
        reader = csv.reader(open(file_path, 'r'), delimiter="\t")

    return reader

def read_10X(matrix_path, feature_path, barcode_path):
    """ Reads 10X data from the files typically outputted from CellRanger, as
        indicated here:
        https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices

    Args:
        matrix_path (str): Path to the matrix.

        feature_path (str): Path to the feature (genes) information.

        barcode_path (str): Path to the barcode (cells) information.

    Returns:
        pandas.DataFrame: Dataframe of the counts, but each column is stored \
                                                               as a SparseArray.
    """
    counts_1_10X = scipy.io.mmread(matrix_path)

    feature_reader = getReader(feature_path)
    barcode_reader = getReader(barcode_path)

    gene_names = [row[1] for row in feature_reader]
    cell_names = [row[0] for row in barcode_reader]

    counts_frame = pandas.DataFrame.sparse.from_spmatrix(counts_1_10X,
                                           index=gene_names, columns=cell_names)

    return counts_frame





