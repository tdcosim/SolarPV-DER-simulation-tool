# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:20:55 2020

@author: splathottam
"""

DER_basic_spec = {'SolarPV_DER_SinglePhase':{'n_phases':1,'n_ODE':11,'t_stable':0.5,'m_steady_state':0.96},
                  'SolarPV_DER_ThreePhase':{'n_phases':3,'n_ODE':23,'t_stable':0.5,'m_steady_state':0.96},
                  'SolarPV_DER_ThreePhaseBalanced':{'n_phases':3,'n_ODE':11,'t_stable':0.5,'m_steady_state':0.96},
                  'SolarPV_DER_SinglePhaseConstantVdc':{'n_phases':1,'n_ODE':10,'t_stable':0.5,'m_steady_state':0.96}
                 }

solver_spec = {'SLSQP':{'ftol': 1e-10, 'disp': True, 'maxiter':10000},
               'nelder-mead':{'xtol': 1e-8, 'disp': True, 'maxiter':10000}}