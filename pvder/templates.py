# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:19:16 2020

@author: splathottam
"""
import six
from pvder.grid_components import Grid,BaseValues
 
DER_design_template = {"basic_specs":["n_phases","Sinsol"],
                       "module_parameters":['Np','Ns','Vdcmpp0','Vdcmpp_min','Vdcmpp_max'],
                       "inverter_ratings":["Srated","Vdcrated","Ioverload","Vrmsrated"],
                       "circuit_parameters":["Rf_actual","Lf_actual","C_actual","R1_actual","X1_actual"],
                       "controller_gains":["Kp_GCC","Ki_GCC","Kp_DC","Ki_DC","Kp_Q","Ki_Q","wp"],
                       "steadystate_values":["maR0","maI0","iaR0","iaI0"],
                       "initial_values":["iaR0","iaI0","xaR0","xaI0","uaR0","uaI0","xDC0","xQ0","xPLL0","wte0"]
                      }
if six.PY3:
    string_type = (str)
elif six.PY2:
    string_type = (str,unicode)
           
DER_argument_template = {'derId':{'default_value':None,'type':string_type},
                         'powerRating':{'default_value':None,'type':(int,float)},
                         'VrmsRating':{'default_value':None,'type':(int,float)},'Vdcrated':{'default_value':None,'type':(int,float)},
                         'gridModel':{'default_value':None,'type':Grid},
                         'verbosity':{'default_value':'INFO','type':string_type},'identifier':{'default_value':'','type':string_type},
                         'derConfig':{'default_value':{},'type':dict},
                         'gridVoltagePhaseA':{'default_value':None,'type':complex},
                         'gridVoltagePhaseB':{'default_value':None,'type':complex},
                         'gridVoltagePhaseC':{'default_value':None,'type':complex},
                             'gridFrequency':{'default_value':None,'type':(int,float)},
                             'standAlone':{'default_value':True,'type':bool},'steadyStateInitialization':{'default_value':True,'type':bool},
                             'allowUnbalancedM':{'default_value':False,'type':bool},
                             'ia0':{'default_value':None,'type':complex},'xa0':{'default_value':None,'type':complex},
                             'ua0':{'default_value':None,'type':complex},'xDC0':{'default_value':None,'type':complex},
                             'xQ0':{'default_value':None,'type':complex}}

#Voltage and frequency ride through settings from IEEE 1557-2018 Category III (Table 16, page 48) 
#V1 to V2 - zone 2,V1 < - zone 1 
    
RT_config_template = {'LVRT':{'0':['V_threshold','t_threshold','mode','t_start','threshold_breach']},
                      'HVRT':{'0':['V_threshold','t_threshold','mode','t_start','threshold_breach']},
                      'OUTPUT_CESSATION_DELAY':'','OUTPUT_RESTORE_DELAY':'','RESTORE_Vdc':'',
                      'F_LF1':'','F_LF2':'','t_LF1_limit':'','t_LF2_limit':'',
                      'F_HF1':'','F_HF2':'',
                      't_HF1_limit':'','t_HF2_limit':'',
                      'FRT_INSTANTANEOUS_TRIP':''
                     }