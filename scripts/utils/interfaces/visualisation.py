""" Interface for visualisations.
"""

from abc import ABC, abstractmethod

class VisualInterface(ABC):

	@abstractmethod
	def fit(self, data):
		"""Fit the visualisation to the data, learn axis, colors, etc."""
		pass

	@abstractmethod
	def visualise(self, data):
		"""Actually create the visualisation and handle it."""
		pass

	@abstractmethod
	def fit_visualise(self, data):
		"""Fits the visualisation to the data and also creates it."""
		pass




