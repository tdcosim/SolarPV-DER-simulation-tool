"""Classes with commonly used methods."""

from __future__ import division

import sys
import time
import six
import logging
import pprint

from pvder import templates,specifications

class Logging(object):
    """ Utility class for common methods."""
    
    pp = pprint.PrettyPrinter(indent=4)
    
    def name_instance(self,identifier=''):
        #Provide a name to the instance.
        
        self.ID = self.count  #ID is same as current instance count
        
        model_name = type(self).__name__
        
        if model_name in templates.DER_design_template:
            self.name = model_name + '_' + str(self.ID) #Object name            
        elif type(self).__name__ == 'DynamicSimulation':
            self.name = 'sim_'+str(self.ID)  #Object name
        elif type(self).__name__ == 'SimulationEvents':
            self.name = 'events_'+str(self.ID)  #Object name
        elif type(self).__name__ == 'SimulationResults':
            self.name = 'results_'+str(self.ID)  #Object name
        else:
            raise ValueError('{} is not a valid instance name.'.format(type(self).__name__))
            
        self.name  = str(identifier) + '-' +self.name  #Add additional identifier to object name if it was provided
    
    def initialize_logger(self,logging_level):
        """Initialize loggers for different classes."""        
        
        self.logger = logging.getLogger(type(self).__name__)       
        
        self.verbosity = logging_level
    
    @property
    def verbosity(self):
        return self.__verbosity
    
    @verbosity.setter
    def verbosity(self,verbosity):
        """Method to set verbosity of logging on terminal. Different classes may have different levels of verbosity."""
        
        if isinstance(verbosity,specifications.string_type):
            if verbosity not in specifications.logging_levels:
                raise ValueError('{} is not a valid logging level!'.format(verbosity))                
            self.__verbosity = verbosity.upper()
        elif isinstance(verbosity,int):
            if verbosity not in specifications.logging_levels_integer:
                raise ValueError('{} is not a valid logging level!'.format(verbosity))
            self.__verbosity = logging.getLevelName(verbosity)
        else:
            raise ValueError('{} is an valid type for verbosity!'.format(type(verbosity)))
        
        #Set logging level - {DEBUG,INFO,WARNING,ERROR.CRITICAL}
        if self.__verbosity == 'DEBUG':
            self.logger.setLevel(logging.DEBUG)            
        
        elif self.__verbosity == 'INFO':
            self.logger.setLevel(logging.INFO)            
            
        elif self.__verbosity == 'WARNING':
            self.logger.setLevel(logging.WARNING)
        
        elif self.__verbosity == 'ERROR':
            self.logger.setLevel(logging.ERROR)
        
        elif self.__verbosity == 'CRITICAL':
            self.logger.setLevel(logging.CRITICAL)
            
        self.logger.debug('{}:Logging level is set to:{}'.format(self.name,self.__verbosity))
                    
        return self.__verbosity
