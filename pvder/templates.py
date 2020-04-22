# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:19:16 2020

@author: splathottam
"""

parameter_properties = {"Rf":{'base':'impendance','type':(int,float)},
                        "Lf":{'base':'inductance','type':(int,float)},
                        "C":{'base':'capacitance','type':(int,float)},
                        "R1":{'base':'impendance','type':(int,float)},
                        "X1":{'base':'impendance','type':(int,float)},
                        "Zf":{'base':'impendance','type':(complex)},
                        "Srated":{'base':'power','type':(int,float)}
                       }

phase_properties =    {"i":{'base':'current','type':(complex,float),'description':'Current'},
                       "x":{'base':'','type':(complex,float),'description':'Controller state'},
                       "u":{'base':'','type':(complex,float),'description':'Controller state'}
                       }
                       
state_properties = {"iaR":{'physical_type':'real'},
                    "iaI":{'physical_type':'imag'},
                    "ibR":{'physical_type':'real'},
                    "ibI":{'physical_type':'imag'},
                    "icR":{'physical_type':'real'},
                    "icI":{'physical_type':'imag'},
                    
                    "xaR":{'physical_type':'real'},
                    "xaI":{'physical_type':'imag'},
                    "xbR":{'physical_type':'real'},
                    "xbI":{'physical_type':'imag'},
                    "xcR":{'physical_type':'real'},
                    "xcI":{'physical_type':'imag'},
                    
                    "uaR":{'physical_type':'real'},
                    "uaI":{'physical_type':'imag'},
                    "ubR":{'physical_type':'real'},
                    "ubI":{'physical_type':'imag'},
                    "ucR":{'physical_type':'real'},
                    "ucI":{'physical_type':'imag'}, 
                    
                    "Vdc":{'physical_type':'real'},
                    
                    "xDC":{'physical_type':'real'},
                    "xP":{'physical_type':'real'},
                    "xQ":{'physical_type':'real'},
                    
                    "xPLL":{'physical_type':'real'},
                    "wte":{'physical_type':'real'}
                   }

controller_properties = {"current_controller":{"gains":["Kp_GCC","Ki_GCC","wp"],"description":"Current controller"},
                         "dc_link_voltage_controller": {"gains":["Kp_DC","Ki_DC"],"description":"DC link voltage controller"},
                         "active_power_controller": {"gains":["Kp_P","Ki_P"],"description":"Active power controller"},
                         "reactive_power_controller": {"gains":["Kp_Q","Ki_Q"],"description":"Reactive power controller"},
                         "pll_controller": {"gains":["Kp_PLL","Ki_PLL"],"description":"PLL controller"}
                         }

DER_design_template = {"SolarPV_DER_SinglePhase":
                       {"basic_specs":{'phases':('a'),'n_phases':1,'n_ODE':11},
                        "basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0},
                       "module_parameters":{"Np":2,"Ns":1000,"Vdcmpp0":750.0,"Vdcmpp_min": 650.0,"Vdcmpp_max": 800.0},
                       "inverter_ratings":{"Srated":10e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0},
                       "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"C_actual":30.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
                       "controller_gains":{"Kp_GCC":12000.0,"Ki_GCC":4000.0,"Kp_DC":-4.0,"Ki_DC":-20.0,"Kp_Q":0.4,"Ki_Q":20.0,"wp": 20e4},
                       "steadystate_values":{"maR0":0.7,"maI0":0.0,"iaR0":0.5,"iaI0":0.01},
                       "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
                                         "Vdc":550.0,
                                         "xDC":0.0,"xQ":0.0,
                                         "xPLL":0.0,"wte":6.28}
                       },
                       
                       "SolarPV_DER_ThreePhase":
                       {"basic_specs":{'phases':('a','b','c'),'n_phases':3,'n_ODE':23,'unbalanced':True},
                        "basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0},
                        "module_parameters":{"Np":11,"Ns":735,"Vdcmpp0":550.0,"Vdcmpp_min": 525.0,"Vdcmpp_max": 650.0}, 
                        "inverter_ratings":{"Srated":50e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0},
                       "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"C_actual":300.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
                       "controller_gains":{"Kp_GCC":6000.0,"Ki_GCC":2000.0,"Kp_DC":-2.0,"Ki_DC":-10.0,"Kp_Q":0.2,"Ki_Q":10.0,"wp": 20e4},
                       "steadystate_values":{"maR0":0.89,"maI0":0.0,"iaR0":1.0,"iaI0":0.001},
                       "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
                                         "ibR":0,"ibI":0.0,"xbR":0.0,"xbI":0.0,"ubR":0.0,"ubI":0.0,
                                         "icR":0,"icI":0.0,"xcR":0.0,"xcI":0.0,"ucR":0.0,"ucI":0.0,
                                         "Vdc":550.0,"xDC":0.0,"xQ":0.0,"xPLL":0.0,"wte":6.28}
                       },
                       
                       "SolarPVDERThreePhaseBalanced":
                       {"basic_specs":{'phases':('a','b','c'),'n_phases':3,'n_ODE':11,'unbalanced':False},
                        "basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0},
                       "module_parameters":{"Np":11,"Ns":735,"Vdcmpp0":550.0,"Vdcmpp_min": 525.0,"Vdcmpp_max": 650.0}, 
                       "inverter_ratings":{"Srated":50e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0},
                       "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"C_actual":300.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
                       "controller_gains":{"Kp_GCC":6000.0,"Ki_GCC":2000.0,"Kp_DC":-2.0,"Ki_DC":-10.0,"Kp_Q":0.2,"Ki_Q":10.0,"wp": 20e4},
                       "steadystate_values":{"maR0":0.89,"maI0":0.0,"iaR0":1.0,"iaI0":0.001},
                       "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
                                         "Vdc":550.0,
                                         "xDC":0.0,"xQ":0.0,
                                         "xPLL":0.0,"wte":6.28}
                        },
                       
                       "SolarPVDER_SinglePhaseConstantVdc":
                       {"basic_specs":{'phases':('a'),'n_phases':1,'n_ODE':10},
                        "basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0,'use_Pref':False},
                       "module_parameters":{"Np":2,"Ns":1000,"Vdcmpp0":750.0,"Vdcmpp_min": 650.0,"Vdcmpp_max": 800.0},
                       "inverter_ratings":{"Srated":10e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0},
                       "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
                       "controller_gains":{"Kp_GCC":12000.0,"Ki_GCC":4000.0,"Kp_P":5000.0,"Ki_P":500.0,"Kp_Q":0.4,"Ki_Q":20.0,"wp": 20e4},
                       "steadystate_values":{"maR0":0.7,"maI0":0.0,"iaR0":0.5,"iaI0":0.01},
                       "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
                                         "xP":0.0,"xQ":0.0,
                                         "xPLL":0.0,"wte":6.28}
                       }
                       }

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
