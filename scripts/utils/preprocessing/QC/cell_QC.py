""" Contains functions which allow for cell QC, mostly to do with removal of
		cells which are below total count, nGenes expressed, and
		nGenes expressed/total count thresholds.

	Example usage:

	### QC stats ###
	count_gene_ratios, umis, gene_counts = cell_QC.getQCMetrics(counts)

	matplotlib.rcParams.update({'font.size': 8, 'font.weight': 'bold'})
	cutoffs = [0, 1200, 1000]
	sf = f'10X_{i}'
	passed = cell_QC.plotQCHistograms(count_gene_ratios, umis, gene_counts, cutoffs,
						   folder='figure_components/', savePlot=False, suffix=sf,
									  linewidth=2, figsize=(5, 3), showPlot=False
									  )

	stats = cell_QC.getQCStats(count_gene_ratios, umis, gene_counts, cutoffs,
										saveFile=True, folder='data/', suffix=sf)
"""

import numpy, pandas
import matplotlib.pyplot as plt
import scripts.utils.visualisation.histogram.histogram as hist

import importlib
importlib.reload(hist)

def getQCMetrics(counts_df):
	""" Gets quality control statistics from cells, such as the umis/gene ratio,
						  number of umis per cell, the number of genes per cell.

	Args:
		counts_df (pandas.DataFrame): genes*cells dataframe of counts.

	Returns:
		 numpy.array<float>, numpy.array<int>, numpy.array<int>: \
			umis/gene ratio, number of umis per cell, \
			the number of genes per cell.
	"""
	umis = numpy.apply_along_axis(sum, 0, counts_df.values)
	gene_counts = numpy.apply_along_axis(getGeneCounts, 0, counts_df.values)
	count_gene_ratios = umis/gene_counts
	return count_gene_ratios, umis, gene_counts

def getGeneCounts(cellExprs):
	""" Returns the number of non-zero entires in the inputted \
		array.
	"""
	return len(numpy.where(cellExprs>0)[0])

def getQualityCells(count_gene_ratios, cell_counts, gene_counts, cutoffs):
	""" Retrieves a boolean array indicating which cells pass the QC cutoffs.

	Args:
		count_gene_ratios:
		cell_counts:
		gene_counts:
		cutoffs:

	Returns:
		numpy.array<bool>: True if passed, false otherwise.
	"""
	metrics = [count_gene_ratios, cell_counts, gene_counts]
	pass_bool = numpy.array([True] * len(gene_counts))
	for i, cutoff in enumerate(cutoffs):
		if cutoff==0:
			continue
		else:
			metric_bool = metrics[i] > cutoff
			pass_bool = numpy.logical_and(pass_bool, metric_bool)

	return pass_bool

def plotQCHistograms(count_gene_ratios, cell_counts, gene_counts, cutoffs,
					  folder='./', suffix='', savePlot=False,
					 linewidth=2, figsize=(5, 3), obs='Cell',
                     side_by_side=False,
                     umiValRange=None, geneValRange=None, ratioValRange=None,
                     **kwargs):
	""" Plots a series of histograms associated with the QC. \
		This include histograms of umis, genes expressed, and the \
		umis/gene ratio.

	Args:
		count_gene_ratios (numpy.array<float>): Ratio of umis/gene.

		cell_counts (numpy.array<int>): Total umis per cell.

		gene_counts (numpy.array<int>): Genes expressed per cell.

		cutoffs (list<float, int, int>): Cutoffs used for each QC metric for \
							ratios, counts, gene_counts, respectively. Use \
							None if don't want to use that metric for cutoff.

		folder (str): Where to save the plots to. Default is current directory.

		suffix (str): The suffix to use for each file ending, example \
					counts_{suffix}.pdf.

		savePlot (bool): True to save the plots, false otherwise.

		linewidth (float): Indicates the thickness of the cutoff bars.

		kwargs (dict): Extra arguments parsed to hist.visualise

	Returns:
		numpy.array<int>: Indices of cells which passed the quality control.
	"""
	### Applying the cutoffs ###
	pass_bool = getQualityCells(count_gene_ratios, cell_counts,
								gene_counts, cutoffs)
	ratios_pass = count_gene_ratios[pass_bool]
	counts_pass = cell_counts[pass_bool]
	genes_pass = gene_counts[pass_bool]

	### Getting the figs & axes ###
	if side_by_side:
		fig, axs = plt.subplots(1, 3, figsize=figsize)
		fig2, axs2 = plt.subplots(1, 3, figsize=figsize)
	else:
		fig, axs = None, [None]*3
		fig2, axs2 = None, [None]*3

	### The Counts/Gene ratio histogram ###
	hist_obj = hist.Histogram()
	hist_obj.fit(count_gene_ratios, cutoff=cutoffs[0] if cutoffs[0]!=0 else None,
                 valRange=ratioValRange)
	hist_obj.visualise(count_gene_ratios, color='gold',
					   x_label=f'UMIs/Gene Per {obs}', y_label='frequency',
					   folder=folder, savePlot=savePlot, figsize=figsize,
				  plotName=f'umiPerGene_hist_{suffix}.pdf',
					   linewidth=linewidth, fig=fig, axes=axs[0], **kwargs)
	hist_obj.visualise(ratios_pass, color='gold',
						x_label=f'UMIs/Gene Per {obs}', y_label='frequency',
						folder=folder, savePlot=savePlot, figsize=figsize,
			 plotName=f'umiPerGene_hist_pass_{suffix}.pdf',
					   linewidth=linewidth, fig=fig2, axes=axs2[0], **kwargs)

	### The Counts per Cell Histogram ###
	hist_obj = hist.Histogram()
	hist_obj.fit(cell_counts, cutoff=cutoffs[1] if cutoffs[1]!=0 else None,
                 valRange=umiValRange)
	hist_obj.visualise(cell_counts,
					   x_label=f'UMIs Per {obs}', y_label='frequency',
					   folder=folder, savePlot=savePlot, figsize=figsize,
					   plotName=f'umi_hist_{suffix}.pdf',
					   linewidth=linewidth, fig=fig, axes=axs[1], **kwargs)
	hist_obj.visualise(counts_pass,
					   x_label=f'UMIs Per {obs}', y_label='frequency',
					   folder=folder, savePlot=savePlot, figsize=figsize,
					plotName=f'umi_hist_pass_{suffix}.pdf',
					   linewidth=linewidth, fig=fig2, axes=axs2[1], **kwargs)

	### Genes Detected per Cell Histogram ###
	hist_obj = hist.Histogram()
	hist_obj.fit(gene_counts, cutoff=cutoffs[2] if cutoffs[2]!=0 else None,
                 valRange=geneValRange)
	hist_obj.visualise(gene_counts, color='mediumseagreen',
					   x_label=f'Genes Per {obs}', y_label='frequency',
					   folder=folder, savePlot=savePlot, figsize=figsize,
				 plotName=f'genePerCell_hist_{suffix}.pdf',
					   linewidth=linewidth, fig=fig, axes=axs[2], **kwargs)
	hist_obj.visualise(genes_pass, color='mediumseagreen',
					   x_label=f'Genes Per {obs}', y_label='frequency',
					   folder=folder, savePlot=savePlot, figsize=figsize,
			plotName=f'genePerCell_hist_pass_{suffix}.pdf',
					   linewidth=linewidth, fig=fig2, axes=axs2[2], **kwargs)

	return pass_bool

def getQCStats_helper(values):
	return [min(values), numpy.median(values), numpy.mean(values), max(values)]

def getQCStats(count_gene_ratios, umis, gene_counts, cutoffs,
			   saveFile=False, folder='./', suffix=''):
	""" Gets quality control statistics.

	Args:
		count_gene_ratios (numpy.array<float>): Ratio of umis/gene.

	   cell_counts (numpy.array<int>): Total umis per cell.

	   gene_counts (numpy.array<int>): Genes expressed per cell.

		cutoffs (list<float, int, int>): Cutoffs used for each QC metric for \
							ratios, counts, gene_counts, respectively. Use \
							None if don't want to use that metric for cutoff.

		saveFile (bool): Whether to save the statistics to a file or not.

		folder (str): Folder to save the file to.

		suffix (str): The suffix to add to the file ending. \
						e.g. f'QC_stats_{suffix}.txt'
	"""
	metrics = {'umi_gene_ratio': count_gene_ratios,
			   'umis': umis, 'gene_counts': gene_counts}
	stats = numpy.zeros((2, 3*4))  # rows are before and after, cols are metrics.
	pass_bool = getQualityCells(count_gene_ratios, umis, gene_counts, cutoffs)
	j = 0
	stat_metrics = ['minimum','median','mean','maximum']
	columns = []
	for i, metric in enumerate(metrics):
		before = metrics[metric]
		after = before[pass_bool]
		stats_before = getQCStats_helper(before)
		stats_after = getQCStats_helper(after)
		stats[0,j:j+4] = stats_before
		stats[1, j:j+4] = stats_after

		columns.extend( [f'{metric}_{stat}' for stat in stat_metrics] )

		j+=4

	### Adding in the number of cells detected ###
	cutoff_info = numpy.array( [[0]*len(cutoffs), cutoffs] )
	cutoff_cols = [f'{metric}_cutoff' for cutoff in list(metrics.keys())]
	ncells = numpy.array([len(umis), len(numpy.where(pass_bool)[0])])
	stats = numpy.concatenate((
							  ncells.reshape(-1,1), cutoff_info, stats), axis=1)
	columns = ['nCells']+cutoff_cols+columns

	stats = pandas.DataFrame(stats,
							 index=['before', 'after'], columns=columns)

	if saveFile:
		stats.to_csv(f'{folder}QC_stats_{suffix}.txt', sep='\t')

	return stats





