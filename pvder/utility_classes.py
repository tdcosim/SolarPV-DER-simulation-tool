"""Classes with commonly used methods."""

from __future__ import division

import sys
import time
import six
import logging

class Logging(object):
    """ Utility class for common methods."""
    
    logging_levels = ['DEBUG','INFO','WARNING']#,'ERROR'
    
    def name_instance(self,identifier):
        #Provide a name to the instance.
        
        self.ID = self.count  #ID is same as current instance count
        
        if type(self).__name__ == 'SolarPV_DER_SinglePhase':
            self.name = 'PVDER-1ph_'+str(self.ID)  #Object name
        elif type(self).__name__ == 'SolarPV_DER_ThreePhase':
            self.name = 'PVDER-3ph_'+str(self.ID)  #Object name
        elif type(self).__name__ == 'DynamicSimulation':
            self.name = 'sim_'+str(self.ID)  #Object name
        elif type(self).__name__ == 'SimulationEvents':
            self.name = 'events_'+str(self.ID)  #Object name
        elif type(self).__name__ == 'SimulationResults':
            self.name = 'results_'+str(self.ID)  #Object name
        else:
            raise ValueError('{} is not a valid instance name.'.format(type(self).__name__))
            
        if identifier is not None:
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
        
        assert verbosity in self.logging_levels, '{} is not a valid logging level!'.format(verbosity)
        
        self.__verbosity = verbosity
        
        #Set logging level - {DEBUG,INFO,WARNING,ERROR}
        if verbosity == 'DEBUG':
            self.logger.setLevel(logging.DEBUG)            
        
        elif verbosity == 'INFO':
            self.logger.setLevel(logging.INFO)            
            
        elif verbosity == 'WARNING':
            self.logger.setLevel(logging.WARNING)
            
        self.logger.debug('{}:Logging level is set to:{}'.format(self.name,verbosity))
                    
        return self.__verbosity
