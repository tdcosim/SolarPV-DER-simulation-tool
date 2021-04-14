# -*- coding: utf-8 -*-
from __future__ import division
"""
Created on Wed Mar 25 09:36:27 2020

@author: splathottam
"""

"""PV-DER base class."""

import os
import math
import pdb

import numpy as np
import scipy
from scipy.optimize import fsolve, minimize

from pvder.DER_check_and_initialize import PVDER_SetupUtilities
from pvder.DER_features import PVDER_SmartFeatures
from pvder.DER_utilities import PVDER_ModelUtilities
from pvder.grid_components import BaseValues
from pvder import utility_functions
from pvder import defaults,templates,specifications
from pvder.logutil import LogUtil


class SolarPVDER(PVDER_SetupUtilities,PVDER_SmartFeatures,PVDER_ModelUtilities,BaseValues):
	"""
	Class for describing a Solar Photo-voltaic Distributed Energy Resource consisting of panel, converters, and
	control systems.
	
	Attributes:
		  count (int): Number of instances of `SolarPV_DER_SinglePhase`.
		  n_ODE (int): Number of ODE's.
	
	"""
	
	#PLL controller parameters
	Kp_PLL = 180 #1800
	Ki_PLL = 320 #32000	
   
	winv = we = 2.0*math.pi*60.0 #Frequency of fundamental waveform
	fswitching  = 10e3 #Inverter switching frequency (not used by model)

	def setup_DER(self,events,configFile,**kwargs):
		"""Setup pvder instance"""
		try:
			self.events = events
			self.DER_model_type = type(self).__name__ 
			
			self.parameter_ID = self.get_DER_id(**kwargs)
			self.create_template()
			
			DER_config,config_dict = self.get_DER_config(configFile,self.parameter_ID)
			
			DER_arguments = self.get_DER_arguments(DER_config,**kwargs)   
			
			self.name_instance(DER_arguments['identifier']) #Generate a name for the instance  
			self.initialize_logger(DER_arguments['verbosity'])  #Set logging level - {DEBUG,INFO,WARNING,ERROR}
			
			self.DER_parent_ID = self.get_DER_parent_id(DER_config)
			DER_parent_config = self.get_DER_parent_config(configFile,self.DER_parent_ID)
			
			for DER_component in self.DER_design_template:
				if DER_component not in DER_config:
					DER_config.update({DER_component:{}})
				if DER_component not in DER_parent_config:
					DER_parent_config.update({DER_component:{}})
			
			self.check_model_type(DER_config,DER_parent_config)
			self.update_DER_config(DER_config,DER_parent_config,DER_arguments,self.parameter_ID,self.DER_parent_ID)
			self.check_basic_specs()
			
			self.RT_config = {}
			self.update_RT_config(config_dict) #Checks and updates RT_config if any entries are missing
			self.check_RT_config()
			
			return DER_arguments
		except:
			LogUtil.exception_handler()

	def get_DER_parent_id(self,DER_config):
		"""Check if user has specified a parent configuration."""
		try:
			if 'parent_config' in DER_config:
				DER_parent_id = DER_config['parent_config']
			else:
				DER_parent_id = ''
			
			return DER_parent_id
		except:
			LogUtil.exception_handler()
	
	def get_DER_parent_config(self,configFile,derParentId):
		"""Check if user has specified a parent configuration."""		
		try:
			if derParentId:
				LogUtil.logger.log(20,'{}:Reading parent DER config:{} for DER config:{}'.format(self.name,derParentId,self.parameter_ID)) 
				DER_parent_config,_ = self.get_DER_config(configFile,derParentId)	
			else:
				DER_parent_config = {}

			return DER_parent_config
		except:
			LogUtil.exception_handler()
	
	def get_DER_config(self,configFile,derId):
		"""Check DER ID in config file."""
		try:
			config_dict = self.read_config(configFile) #Read configuration dictionary
			
			available_ids = list(config_dict.keys())
			if derId in available_ids:
				LogUtil.logger.info('DER configuration with ID:{} was found in {}'.format(derId,configFile))
			else:
				LogUtil.logger.error('DER configuration with ID:{} could not be found in {}! - Available IDs are:{}'.format(derId,configFile,available_ids))
			return config_dict[derId],config_dict
		except:
			LogUtil.exception_handler()
	
	def create_template(self):
		"""Create templates for DER model."""
		try:
			self.DER_design_template = templates.DER_design_template[self.DER_model_type]
			self.DER_config =  dict((key, {}) for key in self.DER_design_template.keys())
		except:
			LogUtil.exception_handler()
		
	def get_DER_id(self,**kwargs):
		"""Create a parameter ID from inverter rated power output.		
		
		Args:
			 powerRating (float): The inverter power rating in kVA.
			 derId (str): User specified parameter ID (can be None). 

		Returns:
			 str: Parameter id
		"""	   
		try:
			assert 'derId' in kwargs or 'powerRating' in kwargs,\
			'Neither DER parameter ID nor DER power rating was provided!'
			if 'derId' in kwargs:
				DER_id = kwargs['derId']   
			elif 'powerRating' in kwargs:
				DER_id = str(int(kwargs['powerRating']/1e3)) 
			return DER_id
		except:
			LogUtil.exception_handler()
	
	def get_DER_arguments(self,DER_config,**kwargs):
		"""Initialize flags"""
		try:
			DER_arguments = {}		
			found_arguments = kwargs.keys() #Arguments which were passed
			used_arguments =[]
		   
			for key, value in specifications.DER_argument_spec .items():
				
				if key in kwargs.keys():				
					if isinstance(kwargs[key],specifications.DER_argument_spec [key]['type']):
					   DER_arguments.update({key:kwargs[key]})
					   used_arguments.append(key)
					else:
						LogUtil.logger.log(40,'Found {} to have type:{} - Valid type:{}'.format(key,type(kwargs[key]),specifications.DER_argument_spec [key]['type']))
				elif key in DER_config: #Check if key available in config file
					if isinstance(DER_config[key],specifications.DER_argument_spec[key]['type']):
						DER_arguments.update({key:DER_config[key]})
						LogUtil.logger.info('Updated DER argument {} from DER_config'.format(key))	   
					else:
						LogUtil.logger.info('Found {} to have type:{} - Valid type:{}'.format(key,type(DER_config[key]),specifications.DER_argument_spec [key]['type']))
				elif key in specifications.DER_argument_spec:
					if specifications.DER_argument_spec [key]['default_value'] is not None:
						DER_arguments.update({key:specifications.DER_argument_spec [key]['default_value']})

			LogUtil.logger.log(10,'Used arguments:{} and Invalid arguments:{}'.format(\
			used_arguments,list(set(found_arguments).difference(set(used_arguments)))))
			
			if 'powerRating' in found_arguments:
				DER_arguments.update({'Srated':kwargs['powerRating']})
				
			if 'VrmsRating' in found_arguments:
				DER_arguments.update({'Vrmsrated':kwargs['VrmsRating']})  
						
			return DER_arguments
		except:
			LogUtil.exception_handler()
	
	def initialize_DER(self,DER_arguments):
		"""Initialize flags"""
		try:
			self.initialize_basic_specs()		
			self.initialize_basic_options()
			
			self.initialize_flags(DER_arguments)
			self.attach_grid_model(DER_arguments)
			
			self.initialize_grid_measurements(DER_arguments)
			self.initialize_DER_model() #DER model parameters
			self.RT_initialize() #VRT and FRT settings  
			
			self.initialize_jacobian()
			
			self.update_Qref(t=0.0) #Reference
			
			self.initialize_states(DER_arguments) #initialize_states
			
			self.initialize_derived_quantities()
			self.initialize_Volt_VAR() #Volt-VAR settings		
			
			self.reset_reference_counters()
			
			if self.standAlone:
				self._vag_previous = self.grid_model.vag
			self._va_previous = self.va
		except:
			LogUtil.exception_handler()
	
	def initialize_flags(self,DER_arguments):
		"""Initialize flags"""
		try:
			self.standAlone = DER_arguments['standAlone']
			self.steady_state_initialization = DER_arguments['steadyStateInitialization']
			self.allow_unbalanced_m = DER_arguments['allowUnbalancedM']
		except:
			LogUtil.exception_handler()
	
	def check_model_type(self,DER_config,DER_parent_config):
		"""Check basic specs in configuration."""
		try:
			if 'model_type' in DER_config['basic_specs']:
				DER_model_type_in_config = DER_config['basic_specs']['model_type']
			elif 'model_type' in DER_parent_config['basic_specs']:
				DER_model_type_in_config = DER_parent_config['basic_specs']['model_type']
			else:
				LogUtil.logger.log(40,'{}:Model type was not found for parameter ID {}!'.format(self.name,self.parameter_ID))
			
			if not self.DER_model_type == DER_model_type_in_config:
				LogUtil.logger.log(40,'{}:DER configuration with ID:{} is defined to be used with model type {} but is being used with model {}!'.format(self.name,self.parameter_ID,DER_model_type_in_config,self.DER_model_type))
		except:
			LogUtil.exception_handler()

	def check_basic_specs(self):
		"""Check basic specs in DER config."""
		try:
			if self.DER_model_type in templates.DER_design_template:
				
				n_phases = self.DER_config['basic_specs']['n_phases']
				if not n_phases == templates.DER_design_template[self.DER_model_type]['basic_specs']['n_phases']:
					LogUtil.logger.log(40,'{}:DER configuration with ID:{} has {} phases which is invalid for {} DER model!'.format(self.name,self.parameter_ID,n_phases,self.DER_model_type))
				
				if not n_phases == len(templates.DER_design_template[self.DER_model_type]['basic_specs']['phases']):
					LogUtil.logger.log(40,'{}:DER configuration with ID:{} has {} phases buf following phases were found {}!'.format(self.name,self.parameter_ID,n_phases,len(templates.DER_design_template[self.DER_model_type]['basic_specs']['phases'])))
				
				n_ODE = self.DER_config['basic_specs']['n_ODE']
				if not n_ODE == templates.DER_design_template[self.DER_model_type]['basic_specs']['n_ODE']:
					LogUtil.logger.log(40,'{}:DER configuration with ID:{} has {} ODE equations which is invalid for {} DER model!'.format(self.name,self.parameter_ID,n_ODE,self.DER_model_type))
				
				if not n_ODE == len(templates.DER_design_template[self.DER_model_type]['initial_states']):
					LogUtil.logger.log(40,'{}:DER configuration with ID:{} needs {} states, but only {} states were found for {} DER model!'.
									 format(self.name,self.parameter_ID,n_ODE,len(templates.DER_design_template[self.DER_model_type]['initial_states']),self.DER_model_type))
						  
			else:
				LogUtil.logger.log(40,'{}:{} is an invalid DER model class'.format(self.name,self.DER_model_type)) 
		except:
			LogUtil.exception_handler()

	def read_config(self,configFile):
		"""Load config json file and return dictionary."""
		try:
			LogUtil.logger.log(10,'Reading configuration file:{}'.format(configFile))
			confDict = utility_functions.read_json(configFile)
			
			return confDict
		except:
			LogUtil.exception_handler()

class PVModule(object):
	"""
	Class for describing PV module.
	
	Attributes:
		Iph (float):Photocurrent from a single cell.
		Ipv (float): PV module current.
	
	"""
	
	#Select either NN or polyfit model for MPP
	USE_POLYNOMIAL_MPP = True
	_MPP_fit_points = 50
	_Tactual_min = 273.15 #0 degree celsius
	_S_min = 10.0
	_Tactual_max = _Tactual_min + 35.0 #35 degree celsius
	_S_max = 100.0
	
	Iscr = 8.03 #Cell short-circuit current at reference temperature and radiation
	Kv = 0.0017  #Short-circuit current temperature co-efficient 
	T0 = 273.15 + 25.0 #Cell reference temperature in Kelvin
	Irs = 1.2e-7		 #Cell reverse saturation current
	q = 1.602e-19		#Charge of an electron
	k = 1.38e-23		#Boltzmann's constant
	A = 1.92	  #p-n junction ideality factor
	
	"""
	module_parameters = {'1':{'Np':2,'Ns':500,'Vdcmpp0':250.0,'Vdcmpp_min': 225.0,'Vdcmpp_max': 300.0},
						 '10':{'Np':2,'Ns':1000,'Vdcmpp0':750.0,'Vdcmpp_min': 650.0,'Vdcmpp_max': 800.0},
						 '50':{'Np':11,'Ns':735,'Vdcmpp0':550.0,'Vdcmpp_min': 520.0,'Vdcmpp_max': 650.0},
						 '250':{'Np':45,'Ns':1000,'Vdcmpp0':750.0,'Vdcmpp_min': 750.0,'Vdcmpp_max': 1000.0}}
	
	module_parameters_list = module_parameters.keys()
	"""
	def __init__(self,Sinsol):
		"""Creates an instance of `PV_Module`.
		Args:
		   Sinsol (float): Solar insolation in percentage.				
		"""
		try:
			super(PVModule,self).__init__()
			self.initialize_module_parameters()
			self.Sinsol = Sinsol #PV module initial conditions
			self.Tactual = defaults.Tactual
			
			#Fit polynomial
			if self.MPPT_ENABLE and self.USE_POLYNOMIAL_MPP:
				self.fit_MPP_poly()
			self.Iph = self.Iph_calc()
		except:
			raise
	
	@property
	def Vdcmpp(self):
		"""Voltage at maximum power point for given insolation and temperature"""
		try:
			#sol, = fsolve(lambda x: 2.5 - np.sqrt(x), 8)
			self.Iph = self.Iph_calc() #This function uses solar insolation
			
			return fsolve(lambda Vdc0:-((self.Np*self.Irs*(scipy.exp((self.q*Vdc0)/(self.k*self.Tactual*self.A*self.Ns))))*(self.q/(self.k*self.Tactual*self.A*self.Ns))*Vdc0)-((self.Np*self.Irs*(scipy.exp((self.q*Vdc0)/(self.k*self.Tactual*self.A*self.Ns))-1)))\
						   +(self.Np*self.Iph),self.Vdcmpp0)[0] #This is a time consuming operation 
		except:
			LogUtil.exception_handler()

	def Iph_calc(self):
		"""Photocurrent from a single cell for given insolation and temperature."""
		try:
			return (self.Iscr+(self.Kv*(self.Tactual-self.T0)))*(self.Sinsol/100.0)
		except:
			LogUtil.exception_handler()

	def Ppv_calc(self,Vdc_actual):
		"""PV panel power output from  solar insolation.
				
		Args:
		  Vdc_actual (float): DC link voltage in volts
		  
		Returns:
			 float: Power output from PV module in p.u.
		"""
		try:
			self.Iph = self.Iph_calc()
			self.Ipv = (self.Np*self.Iph)-(self.Np*self.Irs*(math.exp((self.q*Vdc_actual)/(self.k*self.Tactual*self.A*self.Ns))-1))   #Faster  with Pure Python functions
		
			return max(0,(self.Ipv*Vdc_actual))/BaseValues.Sbase
			#return utility_functions.Ppv_calc(self.Iph,self.Np,self.Ns,Vdc_actual,self.Tactual,Grid.Sbase)
		except:
			LogUtil.exception_handler()

	def fit_MPP_poly(self):
		"""Method to fit MPP to a polynomial function."""
		try:
			self.Tactual =  298.15  #Use constant temperature
			self._MPP_fit_points = 10
			Srange = np.linspace(self._S_min,self._S_max,self._MPP_fit_points+1)
			Vdcmpp_list = []
			Ppvmpp_list =[]
			Sinsol_list= []
			LogUtil.logger.log(20,'Calculating {} values for MPP polynomial fit!'.format(len(Srange)))
			for S in Srange:
				self.Sinsol = S
				Vdcmpp = self.Vdcmpp
				Ppvmpp = self.Ppv_calc(self.Vdcmpp)*BaseValues.Sbase
				
				Sinsol_list.append(S)
				Vdcmpp_list.append(Vdcmpp)
				Ppvmpp_list.append(Ppvmpp)
			
			x = np.array(Sinsol_list)
			y = np.array(Vdcmpp_list)
			self.z = np.polyfit(x, y, 3)
			LogUtil.logger.log(20,'Found polynomial for MPP :{:.4f}x^3 + {:.4f}x^2 +{:.4f}x^1 + {:.4f}!'.format(self.z[0],self.z[1],self.z[2], self.z[3]))		
		except:
			LogUtil.exception_handler()

	def initialize_module_parameters(self):
		"""Initialize PV module parameters."""
		try:
			if self.parameter_ID in self.module_parameters:
				
				self.Np = self.module_parameters[self.parameter_ID]['Np']
				self.Ns = self.module_parameters[self.parameter_ID]['Ns']
				self.Vdcmpp0 = self.module_parameters[self.parameter_ID]['Vdcmpp0']
				self.Vdcmpp_min = self.module_parameters[self.parameter_ID]['Vdcmpp_min']
				self.Vdcmpp_max = self.module_parameters[self.parameter_ID]['Vdcmpp_max']
			else:
				LogUtil.logger.log(40,'PV module parameters not available for parameter ID {} '.format(self.parameter_ID))
		except:
			LogUtil.exception_handler()


