""" Creates plots in matplotlib from reticulated input from R in order to
	make nice plots with good control.
"""

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sb
from sklearn.preprocessing import scale

def density_plot(vals, colors, out_path, figsize=(3, 3),
				 vmin=None, vmax=None, **kwargs):
	""" Plots the density of inputted values.

	Args:
		vals (pd.DataFrame): Keys name of values.
		colors (dict<str, str>): Keys name of values, values are hexcode or named color.
		out_path (str): output path for the figure.
	"""
	matplotlib.rcParams.update({'font.size': 8, 'font.weight': 'bold'})

	vmin = vmin if type(vmin) != type(None) else min(vals.values.ravel())
	vmax = vmax if type(vmax) != type(None) else max(vals.values.ravel())

	fig, ax = plt.subplots(figsize=figsize)
	for label in colors:
		values = vals.loc[:,label].values
		values = values[values>vmin]
		values = values[values<vmax]
		val_norms = values
		# neg_values = values < 0
		# abs_vals = np.abs(values)
		# val_norms = abs_vals/np.sum(abs_vals)
		# val_norms[neg_values] = -val_norms[neg_values]
		# ax.hist(values, color = colors[label],
		# 		histtype='step', cumulative=False, density=True, bins=100)

		sb.distplot(val_norms, hist=False, kde=True,
					#ax=ax,
					color=colors[label], **kwargs)

		# sb.displot(val_norms, kind='kde')

		# df = pd.DataFrame(values)
		# df.plot(df.plot(kind = 'density'), ax=ax)

	ax.spines["right"].set_visible(False)
	ax.spines["top"].set_visible(False)

	save(out_path)

def pie_plot(gene_stats, out_path, figsize=(1,1)):

	labels = gene_stats.loc[:, 'de_status'].values
	label_set = ['up', 'down']
	sizes = [len(np.where(labels == label)[0]) for label in label_set]

	alpha = .7
	fig1, ax1 = plt.subplots(figsize=figsize)
	ax1.pie(sizes,  # labels=label_set,
			colors=[(1, 0, 0, alpha), (0, 0, 1, alpha)],
			shadow=False, startangle=90, )
	ax1.axis(
		'equal')  # Equal aspect ratio ensures that pie is drawn as a circle.


	save(out_path)

def plot_heatmap(gene_stats, lcpms, out_path, alpha=.7, font_size=12, **kwargs):
	""" Plots heatmap of the DE gene expression.
	"""

	matplotlib.rcParams.update({'font.size': font_size, 'font.weight': 'bold'})

	de_genes = gene_stats.index.values[
		gene_stats.loc[:, 'de_status'].values != 'not-de']
	lcpm_de = lcpms.loc[de_genes, :]
	lcpm_z = pd.DataFrame(
		np.apply_along_axis(scale, 1, lcpm_de.values),
		index=de_genes, columns=lcpms.columns)
	sb.clustermap(lcpm_z, row_cluster=True, col_cluster=False,
				   cmap='bwr', col_colors=[(1, 0, 0, alpha), (1, 0, 0, alpha),
										   (0, 0, 1, alpha), (0, 0, 1, alpha)],
				  **kwargs)

	save(out_path)

def plot_gseaClusterMap(gsea_df, out_path,
                        font_size=10, figsize=(8,5), **kwargs):
	""" Plots the clustermap of the GSEA results to visualise gene overlap. """
	matplotlib.rcParams.update({'font.size': font_size, 'font.weight': 'bold'})
	sb.clustermap(gsea_df, row_cluster=False, col_cluster=True,
				  xticklabels=False, cbar_pos=None, figsize=figsize,
					   cmap='bwr', #col_colors=[(1, 0, 0, alpha), (1, 0, 0, alpha),
									#		   (0, 0, 1, alpha), (0, 0, 1, alpha)],
					  **kwargs
				  )
	save(out_path)

def plot_gseaDotPlot(interesting_terms, gsea_results, out_path,
					 font_size=6, figsize=(5,3)):
	""" Plots the dotplots of the nes scores & number of genes in each gene set.
	"""
	matplotlib.rcParams.update({'font.size': font_size, 'font.weight': 'bold'})
	fig, ax = plt.subplots(figsize=figsize)
	plt.setp(ax, yticks=list(range(len(interesting_terms))),
			 yticklabels=interesting_terms[::-1],
			 xlabel='Normalised Enrichment Score (NES)')
	# for i, term in enumerate(interesting_terms[::-1]):
	#     nes = gsea_results.loc[term,'nes']
	#     n_genes = gsea_results.loc[term,'geneset_size']
	#     fdr = gsea_results.loc[term,'fdr']
	nes = gsea_results.loc[interesting_terms[::-1], 'nes'].values
	n_genes = np.array([len(gsea_results.loc[term, 'genes'].split(';'))
						for term in interesting_terms[::-1]])
	ys = list(range(len(interesting_terms)))
	scatter = ax.scatter(nes, ys, s=n_genes,
						 c='red', vmin=0, vmax=1)

	# produce a legend with a cross section of sizes from the scatter
	handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
	legend2 = ax.legend(handles, labels, bbox_to_anchor=(1.01, .2),
						loc="lower left", title="Genes detected")
	plt.tight_layout()
	save(out_path)


def save(out_path):
	out_dir = '/'.join(out_path.split('/')[0:-1]) + '/'
	out_file = out_path.split('/')[-1]
	print(out_dir, out_file)
	dealWithPlot(True, True, True, out_dir, out_file, 300)

def dealWithPlot(savePlot, showPlot, closePlot, folder, plotName, dpi,
				 tightLayout=True):
	""" Deals with the current matplotlib.pyplot.
	"""

	if tightLayout:
		plt.tight_layout()

	if savePlot:
		plt.savefig(folder+plotName, dpi=dpi,
					format=plotName.split('.')[-1])

	if showPlot:
		plt.show()

	if closePlot:
		plt.close()


