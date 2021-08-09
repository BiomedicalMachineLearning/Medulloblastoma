""" Impliments a histogram that looks nice.
"""

import numpy

import scripts.utils.utils.inputHelpers as input
import scripts.utils.interfaces.visualisation as interface

import scripts.utils.visualisation.helpers as vhs

class Histogram(interface.VisualInterface):

	# Defined on object creation #
	bins = 200
	color = 'dodgerblue' #can override when visualising
	showMean = True
	showMedian = True,
	meanColor = 'indianred'
	medianColor = 'darkorange'
	lineHeightRatio = .75
	param_dict = None

	# Defined by fit #
	valRange = None # range of values to plot
	histHeight = None #Total height of the y-axis
	lineHeight = None # the actual height of the verticle lines to plot
	cutoff = None # an extra black verticle line to plot for cutoffs
	title = '' #title of the histogram

	def __init__(self, bins=200, color = 'dodgerblue',
				 showMean=True, showMedian=True,
				 meanColor='indianred', medianColor='darkorange',
				 lineHeightRatio=.75, **param_dict):
		""" Initialises the histogram object for visualisation.

		Args:
			bins (int): Number of bins for histogram.

			color (str): Color of the bins.

			showMean (bool): Whether to plot a verticle line indicating mean.

			showMedian (bool): Whether to plot verticle line indicating median.

			meanColor (str): Color of mean line.

			medianColor (str): Color of median line.

			lineHeightRatio (float): Height of the verticle line, given as fraction \
								of total height of bars.

			param_dict (dict): Extra arguments parsed to matplotlib hist.
		"""
		self.bins = bins
		self.color = color
		self.showMean = showMean
		self.showMedian = showMedian
		self.meanColor = meanColor
		self.medianColor = medianColor
		self.lineHeightRatio = lineHeightRatio
		self.param_dict = input.returnBlankIfNone(param_dict, dict)

	def fit(self, data, cutoff=None, title='', valRange=None):
		""" Fits the histogram to the data.

		Args:
			data (list-like<float>): A sequence of values to plot histogram for.
		"""
		if type(valRange)==type(None):
			self.valRange = (min(data), max(data))
		else:
			self.valRange = valRange
		histVals = numpy.histogram(data, bins=self.bins, range=self.valRange)[0]
		self.histHeight = max(histVals)
		self.histHeight += self.histHeight*0.05
		self.cutoff = cutoff
		self.title = title

	def visualise(self, data, title='', x_label='', y_label='',
				  color=None, fig=None, axes=None, figsize=(10,10), 
				  linewidth=3,
				  savePlot=False, showPlot=True, closePlot=False,
				  folder='./', plotName=None, dpi=300, **kwargs):
		""" Actually creates the visualisation.

		Args:
			data (list-like<float>): A sequence of values to plot histogram for.
		"""
		fig, axes = vhs.getFigAxes(fig, axes, figsize, **kwargs)
		color = input.returnDefaultIfNone(color, self.color)

		# Plotting the histogram
		axes.hist(data, bins=self.bins, color=color, range=self.valRange,
				  **self.param_dict)
		axes.set_ylim(0, self.histHeight)
		axes.set_xlabel( x_label )
		axes.set_ylabel( y_label )
		axes.spines['top'].set_visible(False)
		axes.spines['right'].set_visible(False)

		if self.showMean:
			mean = numpy.mean(data)
			axes.axvline(mean, 0, self.lineHeightRatio,
						 color=self.meanColor, linewidth=linewidth)

		if self.showMedian:
			median = numpy.median(data)
			axes.axvline(median, 0, self.lineHeightRatio,
						 color=self.medianColor, linewidth=linewidth)

		if type(self.cutoff) != type(None):
			axes.axvline(self.cutoff, 0, self.lineHeightRatio,
						 color='black', linewidth=linewidth*(4/3))

		axes.set_title(title)

		# Dealing with the plot #
		vhs.dealWithPlot(savePlot, showPlot, closePlot, folder, plotName, dpi)

	def fit_visualise(self, data, fitArgs=None, visArgs=None, **kwargs):
		""" Performs fit and visualise sequentially.

		Args:
			data (list-like<float>): A sequence of values to plot histogram for.
			fitArgs (dict): Named arguments parsed to fit().
			visArgs (dict): Named arguments parsed to vis().
			kwargs (dict): Named arguments parsed to vis() which flow to \
							matplotlib hist() call.
		"""
		fitArgs = input.returnBlankIfNone(fitArgs, dict)
		visArgs = input.returnBlankIfNone(visArgs, dict)
		self.fit(data, **fitArgs)
		self.visualise(data, **visArgs, **kwargs)









