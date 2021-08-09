""" Miscallenuous but commonly performed calculations.
"""

import numpy, pandas
from tabulate import tabulate
import scripts.utils.utils.inputHelpers as input

import importlib
importlib.reload(input)

def normToOne(values):
	return values/sum(values)

def greaterBool(values, cutoffs):
    return values > cutoffs

def countNonZero(values, cutoff=0):
	return len(numpy.where(values>cutoff)[0])

def binarise(data, cutoff=0):
	""" Takes inputted expressed data and binarises.

	Args:
		data (pandas.DataFrame): Contains int or float data.

	Returns:
		numpy.array<int>: 0 for vals <= cutoff, 1 for vals greater.
	"""
	dataVals = data.values
	boolean = dataVals > cutoff
	binary = boolean.astype(int)
	return binary

def getOrderedLabelSet(labels):
	""" Gets the set of labels associated with labels ordered from most to \
		least frequent.
	"""
	labelSet = list(set(labels))
	freqs = numpy.array([len(numpy.where(labels == label)[0])
						 for label in labelSet])
	order = numpy.array(list((-freqs).argsort()))
	labelSetOrdered = list(numpy.array(labelSet)[order])

	return labelSetOrdered

def filterLowFreqExpr(data, min_cell_expr):
	"""Filters genes which are expressed in less than a given number of cells

	Args:
		data (pandas.DataFrame): genes*cells data frame.

		min_expr_freq (float): Minimum number of cells to express gene for \
									gene to be considered expressed.
	"""
	expr_counts = numpy.apply_along_axis(countNonZero, 1, data.values)
	gene_bool = expr_counts > min_cell_expr
	data = data.loc[gene_bool, :]
	return data

def getWhichPresent(list1, list2):
	"""Returns the items in list 1 that can also be found in list2."""
	present = [item for item in list1
					if item in list2]
	if len(present) != len(list1):
		print("Some genes not present, excluding from selection.")

	return list1

def convertIndicesToBool(indices, length):
	"""Converts a list of indices into a list of bool of specified length,
	with the indices indicated set to True, and the rest False.
	"""
	initBool = numpy.array( [False] * length )
	initBool[indices] = True
	return initBool

def getAvgByLabel(expr, labels, label_set=None):
	""" Gets the average expression of each gene for the different labels of \
		the inputted expression data.

	Args:
		expr (pandas.DataFrame): Genes*Cells data frame.

		labels (numpy.array<str>): Length of samples, indicates a discrete \
									label for each cell.

	Returns:
		pandas.DataFrame: Genes*label_set dataframe, where the number of genes \
							is the same as that inputted, but the label set is \
							the allowed values in the labels.
	"""
	label_set = input.returnDefaultIfNone(label_set, numpy.unique(labels))

	avg_expr = numpy.zeros( (expr.shape[0], len(label_set)) )
	for i, label in enumerate(label_set):
		label_indices = numpy.where(labels == label)[0]
		label_avg_exprs = numpy.apply_along_axis(numpy.mean, 1,
												  expr.values[:, label_indices])
		avg_expr[:, i] = label_avg_exprs

	avg_expr = pandas.DataFrame(avg_expr, index=expr.index, columns=label_set)
	return avg_expr

def printDataFrame(dataframe):
	"""Prints the dataframe in readable format"""
	print( tabulate(dataframe, headers='keys', tablefmt='psql') )