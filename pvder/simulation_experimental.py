# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:48:42 2020

@author: splathottam
"""

import numpy as np

class SimulationExperimental():
    """Class for experimental methods and features."""
    
    def initialize_y0_t(self):
        """Initialize y0_t using dictionary keywords."""        
        
        y0_dict = dict((key+'_t', 0.0) for key in self.PV_model.DER_design_template['initial_states'].keys())
        
        if not len(self.PV_model.y0)==self.PV_model.DER_design_template['basic_specs']['n_ODE']:
            raise ValueError('Expected {} states, but only found {} states!'.format(self.PV_model.DER_design_template['basic_specs']['n_ODE'],len(self.PV_model.y0)))
            
        for i,y0 in enumerate(y0_dict):
            y0_dict[y0] = np.array([self.PV_model.y0[i]])
        
        self.__dict__.update(y0_dict)
    
    
