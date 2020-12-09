"""Classes with commonly used methods."""

from __future__ import division

import sys
import time
import six
import pprint
import pdb
import os

from pvder import templates,specifications
from pvder.defaults import logConfig
from pvder.logutil import LogUtil


class Utilities(object):
	""" Utility class for common methods."""
	pp = pprint.PrettyPrinter(indent=4)

	def name_instance(self,identifier=''):
		"""Generate a name for the model instance."""
		
		self.ID = self.count  #ID is same as current instance count
		model_name = type(self).__name__

		if model_name in templates.DER_design_template:
			self.name = model_name + '_' + str(self.ID) #Object name			
		elif model_name == 'DynamicSimulation':
			self.name = 'sim_'+str(self.ID)  #Object name
		elif model_name == 'SimulationEvents':
			self.name = 'events_'+str(self.ID)  #Object name
		elif model_name == 'SimulationResults':
			self.name = 'results_'+str(self.ID)  #Object name
		else:
			raise ValueError('{} is not a valid instance name.'.format(model_name))
		self.name  = str(identifier) + '-' +self.name  #Add additional identifier to object name if it was provided
	
	def initialize_logger(self,logging_level):
		"""Initialize loggers for different classes."""		
		self.verbosity = logging_level
		
		logConfig['logLevel']= specifications.logging_levels_dict[self.verbosity]		
		LogUtil.set_logger(logConfig['logLevel'],logConfig['mode'])
		
#	@property
#	def verbosity(self):
#		return self.__verbosity

#	@verbosity.setter
#	def verbosity(self,verbosity):
#		"""Method to set verbosity of logging on terminal. Different classes may have different levels of verbosity."""
#		self.__verbosity=10
#		return self.__verbosity


