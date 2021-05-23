""" The purpose of this script is to impliment broadly useful helper functions.
"""

import os, sys
import numpy, pandas

work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/MedullaBlastoma/'

def setUp():
    os.chdir(work_dir)

################################################################################
        # Functions for dealing with multi-species count matrices #
################################################################################
def getSpeciesGenes(gene_matrix_or_array, species_prefix,
                    remove_prefix=False, to_upper=False):
    """ Retrieves the genes belonging to the respective species, assuming that \
    the gene name is prefixed with the inputted species prefix.

    Args:
        gene_matrix_or_array (pandas.DataFrame or numpy.array<str>): In the \
                            case of the former, assumes the gene names are the \
                            index.

        species_prefix (str): Either 'mm10-' or 'hg38-'.

        remove_prefix (bool): Whether to remove the indicates species prefix \
                                from the gene names.

        to_upper (bool): Whether to put the gene names to upper case or not.

    Returns:
         numpy.array<str>: Just the species names.
    """
    if type(gene_matrix_or_array) == pandas.DataFrame:
        genes = gene_matrix_or_array.index.values
    else:
        genes = gene_matrix_or_array

    species_genes = []
    for i, gene in enumerate(genes):
        if gene.startswith(species_prefix):
            if remove_prefix:
                gene = gene.replace(species_prefix, '')
            if to_upper:
                gene = gene.upper()
            species_genes.append( gene )

    return species_genes

def reOrderDF(dataframe, col_name, small_to_large=True):
    """Re-orders the df based on the indicated column and direction."""
    col_values = dataframe.loc[:,col_name].values
    if not small_to_large:
        col_values = -col_values
    order = numpy.argsort(col_values)

    dataframe = dataframe.iloc[order, :]
    return dataframe







