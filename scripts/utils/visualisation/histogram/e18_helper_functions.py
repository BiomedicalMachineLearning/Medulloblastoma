""" The purpose of this script is to provide helper functions to 
	e18QC.py, mostly the plotting and saving of histgrams for the
	QC/Biological filtering supplementary figure.
"""

import matplotlib.pyplot as plt
import numpy

def getGeneCounts(cellExprs):
	""" Returns the number of non-zero entires in the inputted \
		array.
	"""
	return len(numpy.where(cellExprs>0)[0])


def plotHistogram(data, title, cutoff, cutoffY,
		  color='blue', fileName=None, 
		  bins=200, showPlot=True, 
		  plotCutoff=True, yMax=None, **kwargs):
	""" Plots the histogram with a vertical line at a specified 
		position with a given height. Point is to show a QC
		cutoff. 

	Args:
	   data (numpy.array<int or float>): Array containing the data. 

	   title (str): The title to give the outputted plot.

	   cutoff (int): The cutoff to use.

	   cutoffY (int): The height of the vertical line to show the \ 
			   cutoff to use.

	   color (str): The colour of the histogram. 

	   fileName (str or None): None indicates not to save the \
			           file, whereas a string indicates \
				   the name to give to the saved \ 
					figure.
	
	   bins (int): The number of bins to put the data into \
			to show the histogram.

	   showPlot (bool): Whether to show the plot or clear.

	   yMax (int): The limit of the height of the y-axis. 

	   **kwargs: Extra input to matplotlib.pypot.hist
	"""
	plt.hist(data, bins=bins, color=color, **kwargs)
	mean, median = numpy.mean(data), numpy.median(data)
	print(title+' mean:', mean)
	print(title+' median:', median)
	if plotCutoff:
		plt.vlines(cutoff, 0, cutoffY)
	
	if type(yMax) != type(None):
		plt.ylim(0, yMax)	

	plt.vlines(mean, 0, cutoffY*.75, 'indianred')
	plt.vlines(median, 0, cutoffY*.75, 'darkorange')
	plt.title(title)
	#plt.axes(frameon=False)
	plt.savefig(fileName, dpi=300)
	if showPlot:
		plt.show()
	else:
		plt.close()

def plotQCHistograms(cellCounts, geneCounts, countGeneRatios,
		     ranges=None, yMaxs=None, suffix=''):
	""" Plots a series of histograms associated with the QC. \
		This include histograms of umis, genes expressed, and the \
		umis/gene ratio.

	Args:
	   cellCounts (numpy.array<int>): Total umis per cell.

	   geneCounts (numpy.array<int>): Genes expressed per cell.

	   countGeneRatios (numpy.array<float>): Ratio of umis/gene. 

	   ranges (list<tuple(int, int)>): For each of the above, the min and \
				max values to display on the histogram.

	   yMaxs (list<int>): The maximum height of the histograms for \
				each of the above. 

	   suffix (str): The suffix to use for each file ending.

	Returns:
	   list<tuple(int, int)>, list<int>: \
			The ranges used for plotting each histogram. \
			And the heights of the histogram. 
	"""
	if type(ranges) == type(None):
		ranges = []
		for array in [cellCounts, geneCounts, countGeneRatios]:
			ranges.append( (min(array), max(array)) )	

	if type(yMaxs) == type(None):
		yMaxs = [None, None, None]

	### The Counts per Cell Histogram ###
	countsCutoff, countsCutoffY = 600, 1000
	plotHistogram(cellCounts, 'UMIs Per Cell',
	                  countsCutoff, countsCutoffY,
        	          color='dodgerblue',
			  yMax=yMaxs[0],
                	  fileName=f'figures/QC_cellUMIs_Histogram{suffix}.pdf',
                	  showPlot=True, range=ranges[0])

	### Genes Detected per Cell Histogram ###
	genesCutoff, genesCutoffY = 750, 500
	plotHistogram(geneCounts, 'Genes Per Cell',
                  	genesCutoff, genesCutoffY,
                  	color='mediumseagreen',
			yMax=yMaxs[1],
                 	fileName=f'figures/QC_cellGenes_Histogram{suffix}.pdf',
       			plotCutoff=False,
	         	showPlot=True, range=ranges[1])

	### The Counts/Gene ratio histogram ###
	ratioCutoff, ratioCutoffY = 1.2, 250
	plotHistogram(countGeneRatios, 'UMIs/Gene Per Cell',
        	          ratioCutoff, ratioCutoffY,
                	  color='gold',
			  yMax=yMaxs[2],
                 	  fileName=f'figures/QC_UMIsPerGenes_Histogram{suffix}.pdf',
                	  showPlot=True, range=ranges[2])

	return ranges


 	















