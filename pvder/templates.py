# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 10:19:16 2020

@author: splathottam
"""

single_phase_models = ["SolarPVDERSinglePhase","SolarPVDERSinglePhaseConstantVdc"]
three_phase_models = ["SolarPVDERThreePhase","SolarPVDERThreePhaseConstantVdc","SolarPVDERThreePhaseBalanced","SolarPVDERThreePhaseNumba"]

constant_Vdc_models = ["SolarPVDERSinglePhaseConstantVdc","SolarPVDERThreePhaseConstantVdc"]

model_types = ['SinglePhase','SinglePhaseConstantVdc',
			   'ThreePhaseUnbalanced','ThreePhaseUnbalancedConstantVdc','ThreePhaseBalanced','ThreePhaseUnbalancedNumba']

DER_design_template = {"SolarPVDERSinglePhase":
					   {"parent_config":"",						
						"basic_specs":{'phases':('a'),'n_phases':1,'n_ODE':11,"model_type":"SolarPVDERSinglePhase"},
						"basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0,"current_gradient_limiter":False},
					   "module_parameters":{"Np":2,"Ns":1000,"Vdcmpp0":750.0,"Vdcmpp_min": 650.0,"Vdcmpp_max": 800.0},
					   "inverter_ratings":{"Srated":10e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0,"Iramp_max_gradient_real":1.0,"Iramp_max_gradient_imag":1.0},
					   "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"C_actual":30.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
					   "controller_gains":{"Kp_GCC":12000.0,"Ki_GCC":4000.0,"Kp_DC":-4.0,"Ki_DC":-20.0,"Kp_Q":0.4,"Ki_Q":20.0,"wp": 20e4},
					   "steadystate_values":{"maR0":0.7,"maI0":0.0,"iaR0":0.5,"iaI0":0.01},
					   "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
										 "Vdc":550.0,
										 "xDC":0.0,"xQ":0.0,
										 "xPLL":0.0,"wte":6.28}
					   
					   },
					   
					   "SolarPVDERSinglePhaseConstantVdc":
					   {"parent_config":"",
						"basic_specs":{'phases':('a'),'n_phases':1,'n_ODE':10,"model_type":"SolarPVDERSinglePhaseConstantVdc"},
						"basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0,'use_Pref':False,"current_gradient_limiter":False},
					   "module_parameters":{"Np":2,"Ns":1000,"Vdcmpp0":750.0,"Vdcmpp_min": 650.0,"Vdcmpp_max": 800.0},
					   "inverter_ratings":{"Srated":10e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0,"Iramp_max_gradient_real":1.0,"Iramp_max_gradient_imag":1.0},
					   "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
					   "controller_gains":{"Kp_GCC":12000.0,"Ki_GCC":4000.0,"Kp_P":5000.0,"Ki_P":500.0,"Kp_Q":0.4,"Ki_Q":20.0,"wp": 20e4},
					   "steadystate_values":{"maR0":0.7,"maI0":0.0,"iaR0":0.5,"iaI0":0.01},
					   "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
										 "xP":0.0,"xQ":0.0,
										 "xPLL":0.0,"wte":6.28}
					   },
					   
					   "SolarPVDERThreePhase":
					   {"parent_config":"",						
						"basic_specs":{'phases':('a','b','c'),'n_phases':3,'n_ODE':23,'unbalanced':True,"model_type":"SolarPVDERThreePhase"},
						"basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0,"current_gradient_limiter":False},
						"module_parameters":{"Np":11,"Ns":735,"Vdcmpp0":550.0,"Vdcmpp_min": 525.0,"Vdcmpp_max": 650.0}, 
						"inverter_ratings":{"Srated":50e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0,"Iramp_max_gradient_real":1.0,"Iramp_max_gradient_imag":1.0},
					   "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"C_actual":300.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
					   "controller_gains":{"Kp_GCC":6000.0,"Ki_GCC":2000.0,"Kp_DC":-2.0,"Ki_DC":-10.0,"Kp_Q":0.2,"Ki_Q":10.0,"wp": 20e4},
					   "steadystate_values":{"maR0":0.89,"maI0":0.0,"iaR0":1.0,"iaI0":0.001},
					   "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
										 "ibR":0,"ibI":0.0,"xbR":0.0,"xbI":0.0,"ubR":0.0,"ubI":0.0,
										 "icR":0,"icI":0.0,"xcR":0.0,"xcI":0.0,"ucR":0.0,"ucI":0.0,
										 "Vdc":550.0,"xDC":0.0,"xQ":0.0,"xPLL":0.0,"wte":6.28}					  
					   
					   },
					   
					   "SolarPVDERThreePhaseNumba":
					   {"parent_config":"",						
						"basic_specs":{'phases':('a','b','c'),'n_phases':3,'n_ODE':23,'unbalanced':True,"model_type":"SolarPVDERThreePhase"},
						"basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0,"current_gradient_limiter":False},
						"module_parameters":{"Np":11,"Ns":735,"Vdcmpp0":550.0,"Vdcmpp_min": 525.0,"Vdcmpp_max": 650.0}, 
						"inverter_ratings":{"Srated":50e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0,"Iramp_max_gradient_real":1.0,"Iramp_max_gradient_imag":1.0},
					   "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"C_actual":300.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
					   "controller_gains":{"Kp_GCC":6000.0,"Ki_GCC":2000.0,"Kp_DC":-2.0,"Ki_DC":-10.0,"Kp_Q":0.2,"Ki_Q":10.0,"wp": 20e4},
					   "steadystate_values":{"maR0":0.89,"maI0":0.0,"iaR0":1.0,"iaI0":0.001},
					   "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
										 "ibR":0,"ibI":0.0,"xbR":0.0,"xbI":0.0,"ubR":0.0,"ubI":0.0,
										 "icR":0,"icI":0.0,"xcR":0.0,"xcI":0.0,"ucR":0.0,"ucI":0.0,
										 "Vdc":550.0,"xDC":0.0,"xQ":0.0,"xPLL":0.0,"wte":6.28}					  
					   
					   },
					   
					   
					   "SolarPVDERThreePhaseBalanced":
					   {"parent_config":"",						
						"basic_specs":{'phases':('a','b','c'),'n_phases':3,'n_ODE':11,'unbalanced':False,"model_type":"SolarPVDERThreePhaseBalanced"},
						"basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0,"current_gradient_limiter":False},
					   "module_parameters":{"Np":11,"Ns":735,"Vdcmpp0":550.0,"Vdcmpp_min": 525.0,"Vdcmpp_max": 650.0}, 
					   "inverter_ratings":{"Srated":50e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0,"Iramp_max_gradient_real":1.0,"Iramp_max_gradient_imag":1.0},
					   "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"C_actual":300.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
					   "controller_gains":{"Kp_GCC":6000.0,"Ki_GCC":2000.0,"Kp_DC":-2.0,"Ki_DC":-10.0,"Kp_Q":0.2,"Ki_Q":10.0,"wp": 20e4},
					   "steadystate_values":{"maR0":0.89,"maI0":0.0,"iaR0":1.0,"iaI0":0.001},
					   "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
										 "Vdc":550.0,
										 "xDC":0.0,"xQ":0.0,
										 "xPLL":0.0,"wte":6.28}
						},
											  
					   "SolarPVDERThreePhaseConstantVdc":
					   {"parent_config":"",						
						"basic_specs":{'phases':('a','b','c'),'n_phases':3,'n_ODE':22,'unbalanced':True,"model_type":"SolarPVDERThreePhaseConstantVdc"},
						"basic_options":{'t_stable':0.5,'m_steady_state':0.96,"Sinsol":100.0,"current_gradient_limiter":False},
						"module_parameters":{"Np":11,"Ns":735,"Vdcmpp0":550.0,"Vdcmpp_min": 525.0,"Vdcmpp_max": 650.0}, 
						"inverter_ratings":{"Srated":50e3,"Vdcrated":550.0,"Ioverload":1.3,"Vrmsrated":177.0,"Iramp_max_gradient_real":1.0,"Iramp_max_gradient_imag":1.0},
					   "circuit_parameters":{"Rf_actual":0.002,"Lf_actual":25.0e-6,"R1_actual":0.0019,"X1_actual":0.0561},
					   "controller_gains":{"Kp_GCC":6000.0,"Ki_GCC":2000.0,"Kp_P":500.0,"Ki_P":50.0,"Kp_Q":0.2,"Ki_Q":10.0,"wp": 20e4},
					   "steadystate_values":{"maR0":0.89,"maI0":0.0,"iaR0":1.0,"iaI0":0.001},
					   "initial_states":{"iaR":0,"iaI":0.0,"xaR":0.0,"xaI":0.0,"uaR":0.0,"uaI":0.0,
										 "ibR":0,"ibI":0.0,"xbR":0.0,"xbI":0.0,"ubR":0.0,"ubI":0.0,
										 "icR":0,"icI":0.0,"xcR":0.0,"xcI":0.0,"ucR":0.0,"ucI":0.0,
										 "xP":0.0,"xQ":0.0,
										 "xPLL":0.0,"wte":6.28}
					   }
					   }
	

VRT_config_template = {'LVRT':{'config_id':'',
								  'config':{'1':{'V_threshold':0.88,
												 't_threshold':5.0,
                                                 't_min_ridethrough':2,
												 'mode':'mandatory_operation',
												 't_start':0.0,
												 'threshold_breach':False},
											'2':{'V_threshold':0.7,
												 't_threshold':3.5,
                                                 't_min_ridethrough':1,
												 'mode':'mandatory_operation', #'momentary_cessation'
												 't_start':0.0,
												 'threshold_breach':False},
											'3':{'V_threshold':0.5,
												 't_threshold':1.0,
                                                 't_min_ridethrough':0.5,
												 'mode':'momentary_cessation',
												 't_start':0.0,
												 'threshold_breach':False},
								  }
						},
						'HVRT':{'config_id':'',
								  'config':{'1':{'V_threshold':1.06,
												't_threshold':1.0,
                                                't_min_ridethrough':0.5,
												'mode':'mandatory_operation',
												't_start':0.0,
												'threshold_breach':False},
										   '2':{'V_threshold':1.12,
												't_threshold':0.5,
                                                't_min_ridethrough':0.5,
												'mode':'mandatory_operation',
												't_start':0.0,
												'threshold_breach':False},
							  }
						 },
						'VRT_delays':{'config_id':'',
									 'config':{'output_cessation_delay':0.01,
											  'output_restore_delay':1.75,
											  'restore_Vdc':False}
					   }}

FRT_config_template =  {'LFRT':{'parent_config':'',
						   'config':{'1':{'F_LF':57.0,
										  't_LF_limit':1/60,
										  't_LFstart':0.0},
									 '2':{'F_LF':58.8,
										  't_LF_limit':299.0,
										  't_LFstart':0.0}
									}},
						'HFRT':{'parent_config':'',
						   'config':{'1':{'F_HF':61.2,
										  't_HF_limit':299.0,
										  't_HFstart':0.0},
									 '2':{'F_HF':62.0,
										  't_HF_limit':1/60,
										  't_HFstart':0.0}
									}},
						'FRT_delays':{'parent_config':'',
								 'config':{'FRT_INSTANTANEOUS_TRIP':False}}
						}						   
					  
#Voltage and frequency ride through settings from IEEE 1557-2018 Category III (Table 16, page 48) 
#V1 to V2 - zone 2,V1 < - zone 1  SolarPVDER_ThreePhaseConstantVdc
	
RT_config_template_old = {'LVRT':{'0':['V_threshold','t_threshold','mode','t_start','threshold_breach']},
					  'HVRT':{'0':['V_threshold','t_threshold','mode','t_start','threshold_breach']},
					  'OUTPUT_CESSATION_DELAY':'','OUTPUT_RESTORE_DELAY':'','RESTORE_Vdc':'',
					  'F_LF1':'','F_LF2':'','t_LF1_limit':'','t_LF2_limit':'',
					  'F_HF1':'','F_HF2':'',
					  't_HF1_limit':'','t_HF2_limit':'',
					  'FRT_INSTANTANEOUS_TRIP':''
					 }
