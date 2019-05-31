"""Classes with commonly used methods."""

from __future__ import division

import sys
import time
import six
import logging

class Logging():
    """ Utility class for common methods."""
    
    logging_levels = ['DEBUG','INFO']#,'WARNING','ERROR'
    
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
            logging.getLogger().setLevel(logging.DEBUG)
        
        elif verbosity == 'INFO':
            logging.getLogger().setLevel(logging.INFO)
            
        elif verbosity == 'WARNING':
            logging.getLogger().setLevel(logging.INFO)
            
        print('{}:Logging level is set to:{}'.format(self.name,verbosity))
                    
        return self.__verbosity
