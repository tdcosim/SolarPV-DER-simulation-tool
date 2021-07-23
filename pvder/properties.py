# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 12:43:26 2020

@author: splathottam
"""

import six

if six.PY3:
	string_type = (str)
elif six.PY2:
	string_type = (str,unicode)

parameter_properties = {"Rf":{'base':'impendance','type':(int,float)},
						"Lf":{'base':'inductance','type':(int,float)},
						"C":{'base':'capacitance','type':(int,float)},
						"R1":{'base':'impendance','type':(int,float)},
						"X1":{'base':'impendance','type':(int,float)},
						"Zf":{'base':'impendance','type':(complex)},
						
						"Rf_actual":{'base':'impendance','type':(int,float)},
						"Lf_actual":{'base':'impendance','type':(int,float)},
						"C_actual":{'base':'impendance','type':(int,float)},
						"R1_actual":{'base':'impendance','type':(int,float)},
						"X1_actual":{'base':'impendance','type':(int,float)},
						
						"Np":{'type':(int,float)},
						"Ns":{'type':(int,float)},
						"Vdcmpp0":{'base':'voltage','type':(int,float)},
						"Vdcmpp_min":{'base':'voltage','type':(int,float)},
						"Vdcmpp_max":{'base':'voltage','type':(int,float)},
						
						"Srated":{'base':'power','type':(int,float)},
						"Vdcrated":{'base':'voltage','type':(int,float)},
						"Ioverload":{'type':(int,float)},
						"Vrmsrated":{'base':'voltage','type':(int,float)},
						"Iramp_max_gradient_real":{'type':(int,float)},
						"Iramp_max_gradient_imag":{'type':(int,float)},
						
						"Kp_GCC":{'type':(int,float)},
						"Ki_GCC":{'type':(int,float)},
						"Kp_DC":{'type':(int,float)},
						"Ki_DC":{'type':(int,float)},
						"Kp_P":{'type':(int,float)},
						"Ki_P":{'type':(int,float)},						
						"Kp_Q":{'type':(int,float)},
						"Ki_Q":{'type':(int,float)},
						"wp": {'type':(int,float)},
							   
						"maR0":{'type':(int,float)},
						"maI0":{'type':(int,float)},
						"iaR0":{'type':(int,float)},
						"iaI0":{'type':(int,float)}, 
						
						"iaR":{'type':(int,float)},
						"iaI":{'type':(int,float)},
						"ibR":{'type':(int,float)},
						"ibI":{'type':(int,float)},
						"icR":{'type':(int,float)},
						"icI":{'type':(int,float)},
						
						"xaR":{'type':(int,float)},
						"xaI":{'type':(int,float)},
						"xbR":{'type':(int,float)},
						"xbI":{'type':(int,float)},
						"xcR":{'type':(int,float)},
						"xcI":{'type':(int,float)},
						
						"uaR":{'type':(int,float)},
						"uaI":{'type':(int,float)},
						"ubR":{'type':(int,float)},
						"ubI":{'type':(int,float)},
						"ucR":{'type':(int,float)},
						"ucI":{'type':(int,float)},
						
						"Vdc":{'type':(int,float)},
						"xDC":{'type':(int,float)},
						"xP":{'type':(int,float)},
						"xQ":{'type':(int,float)},
						"xPLL":{'type':(int,float)},
						"wte":{'type':(int,float)},
						
						"n_phases":{'type':(int)},
						"n_ODE":{'type':(int)},
						"model_type":{'type':string_type},
						"phases":{'type':(tuple)},
						"unbalanced":{'type':(bool)},
						
						't_stable':{'type':(int,float)},
						'm_steady_state':{'type':(float)},
						"Sinsol":{'type':(int,float)},
						'use_Pref':{'type':(bool)},
						'current_gradient_limiter':{'type':(bool)}
						
						
					   }

phase_properties =	{"i":{'base':'current','type':(complex,float),'description':'Current'},
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