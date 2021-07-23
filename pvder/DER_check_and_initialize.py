"""Code for initializing and validating PV-DER model instances."""

from __future__ import division

import math
import cmath
import numpy as np
from scipy.optimize import fsolve, minimize

from pvder.grid_components import BaseValues
from pvder import utility_functions
from pvder import defaults,templates,specifications,properties
from pvder.logutil import LogUtil

PHASE_DIFFERENCE_120 = 120.0*(math.pi/180.0)


class PVDER_SetupUtilities(BaseValues):
	"""
	 Utility class for error checking during model initialization.
	"""
	module_parameters = {}
	inverter_ratings = {}
	circuit_parameters = {}
	controller_gains ={}
	steadystate_values= {}	
	
	solver_spec = specifications.steadystate_solver_spec
	steadystate_solver = defaults.STEADYSTATE_SOLVER
	
	def creation_message(self):
		"""Message after PV-DER instance was created."""		
		try:
			LogUtil.logger.info('{}:Instance created with DER parameter ID: {}; Specifications - Srated:{:.1f} kVA, Ppv:{:.1f} kW, Vrms:{:.1f} V, Steady state:{},LVRT Enable:{},HVRT Enable:{}'.format(self.name,self.parameter_ID,self.Sinverter_rated/1e3,(self.Ppv*BaseValues.Sbase)/1e3,self.Vrms_rated,self.steady_state_initialization,self.LVRT_ENABLE,self.HVRT_ENABLE))
		except:
			LogUtil.exception_handler()
		
	def update_DER_config(self,DER_config,DER_parent_config,DER_arguments,DER_id,DER_parent_id):
		"""Update PV-DER design."""
		try:
			self.check_DER_config(DER_config,DER_id)
	
			#parameters_from_arguments = {}
			#parameters_from_config = {}
			#parameters_from_default = {}
	
			for DER_component in self.DER_design_template:
				if DER_component in DER_arguments: 
					self.DER_config[DER_component] = DER_arguments[DER_component]
				else:
					for DER_parameter in self.DER_design_template[DER_component]:
						_ = self.update_DER_parameter(DER_config,DER_parent_config,DER_arguments,DER_id,DER_parent_id,DER_component,DER_parameter)
	
			for RT_component in list(templates.VRT_config_template.keys()) + list(templates.FRT_config_template.keys()):
				if RT_component in DER_arguments:
					LogUtil.logger.debug('{}:Updating {} from DER arguments with settings:{}.'.format(self.name,RT_component,DER_arguments[RT_component]))
					self.DER_config[RT_component] = DER_arguments[RT_component]
				elif RT_component in DER_config:
					self.DER_config[RT_component] = DER_config[RT_component]
					LogUtil.logger.debug('{}:Updating {} from DER config {} with settings:{}.'.format(self.name,RT_component,DER_id,DER_config[RT_component]))
				elif RT_component in DER_parent_config: #Check if RT setting exists in parent config file
					LogUtil.logger.debug('{}:Updating {} from parent DER config {} with settings {}.'.format(self.name,RT_component,DER_parent_id,DER_parent_config[RT_component])) 
					self.DER_config[RT_component] = DER_parent_config[RT_component]
				else:
					if RT_component in list(templates.VRT_config_template.keys()):
						default_RT = templates.VRT_config_template[RT_component]
					if RT_component in list(templates.FRT_config_template.keys()):
						default_RT = templates.FRT_config_template[RT_component]
					LogUtil.logger.debug('{}:Updating {} from template with  settings:{}.'.format(self.name,RT_component,default_RT))
					self.DER_config[RT_component] = default_RT
	
			self.basic_specs = {DER_id:self.DER_config['basic_specs']}
			self.module_parameters = {DER_id:self.DER_config['module_parameters']}
			self.inverter_ratings = {DER_id:self.DER_config['inverter_ratings']}
			self.circuit_parameters = {DER_id:self.DER_config['circuit_parameters']}
			self.circuit_parameters[DER_id].update({'Z1_actual':self.DER_config['circuit_parameters']['R1_actual'] + 1j*self.DER_config['circuit_parameters']['X1_actual']})
	
			self.controller_gains = {DER_id:self.DER_config['controller_gains']}
			self.steadystate_values = {DER_id:self.DER_config['steadystate_values']} 
		except:
			LogUtil.exception_handler()

	def update_DER_parameter(self,DER_config,DER_parent_config,DER_arguments,DER_id,DER_parent_id,DER_component,DER_parameter):
		"""Update DER config using both config file and arguments."""
		try:
			parameters ={}
			DER_parameter_type = properties.parameter_properties[DER_parameter]['type']
			if DER_parameter in DER_arguments: #Check if parameter exists in DER key word arguments
				if isinstance(DER_arguments[DER_parameter],DER_parameter_type):
					LogUtil.logger.debug('{}:Parameter {} in {} for ID:{} - updating from DER arguments with value {}.'.format(self.name,DER_parameter,DER_component,DER_id,DER_arguments[DER_parameter])) 
					self.DER_config[DER_component].update({DER_parameter:DER_arguments[DER_parameter]})
					parameters.update({DER_parameter:DER_arguments[DER_parameter]})
					source = 'DER arguments'
				else:
					raise ValueError('Found {} to have type {} - expected type:{}!'.format(DER_parameter,type(DER_arguments[DER_parameter]),DER_parameter_type))
		
			elif DER_parameter in DER_config[DER_component]: #Check if parameter exists in config file
				if isinstance(DER_config[DER_component][DER_parameter],DER_parameter_type):
					LogUtil.logger.debug('{}:Parameter {} in {} for ID:{} - updating from DER config {} with value {}.'.format(self.name,DER_parameter,DER_component,DER_id,DER_id,DER_config[DER_component][DER_parameter])) 
					self.DER_config[DER_component].update({DER_parameter:DER_config[DER_component][DER_parameter]})
					parameters.update({DER_parameter:DER_config[DER_component][DER_parameter]})
					source = 'DER config'
				else:
					raise ValueError('Found {} to have type {} - expected type:{}!'.format(DER_parameter,type(DER_config[DER_component][DER_parameter]),DER_parameter_type))
		
			elif DER_parameter in DER_parent_config[DER_component]: #Check if parameter exists in parent config file
				if isinstance(DER_parent_config[DER_component][DER_parameter],DER_parameter_type):
					LogUtil.logger.debug('{}:Parameter {} in {} for ID:{} - updating from parent DER config {} with value {}.'.format(self.name,DER_parameter,DER_component,DER_id,DER_parent_id,DER_parent_config[DER_component][DER_parameter])) 
					self.DER_config[DER_component].update({DER_parameter:DER_parent_config[DER_component][DER_parameter]})
					parameters.update({DER_parameter:DER_parent_config[DER_component][DER_parameter]})
					source = 'DER parent config'
				else:
					raise ValueError('Found {} to have type {} - expected type:{{}}!'.format(DER_parameter,type(DER_parent_config[DER_component][DER_parameter]),DER_parameter_type))
		
			elif DER_parameter in self.DER_design_template[DER_component]:#Check if parameter exists in default DER config
				LogUtil.logger.debug('{}:Parameter {} in {} for ID:{} - updating from DER template with value {}.'.format(self.name,DER_parameter,DER_component,DER_id,self.DER_design_template[DER_component][DER_parameter])) 
				self.DER_config[DER_component].update({DER_parameter:self.DER_design_template[DER_component][DER_parameter]})
				parameters.update({DER_parameter:self.DER_design_template[DER_component][DER_parameter]})
				source = 'DER default values'
			else:
				raise ValueError('{}: Parameter {} in {} not found for ID:{} - update config file and try again.'.format(self.name,DER_parameter,DER_component,DER_id))
		
			#LogUtil.logger.debug('{}:Following parameters were obtained from {}:{}'.format(self.name,source,parameters))
		
			return parameters
		except:
			LogUtil.exception_handler()


	def check_DER_config(self,DER_config,DER_id):
		"""Check DER config."""
		try:
			available_keys = sorted(list(DER_config.keys()))
			required_keys = sorted(self.DER_design_template.keys())
			missing_keys = set(required_keys)-set(available_keys)
				
			if not bool(missing_keys):
				LogUtil.logger.debug('{}:DER configuration with ID {} contained all required keys.'.format(self.name,DER_id))		
			else:
				if not DER_config['parent_config']:
					raise KeyError('{}:DER configuration with ID:{} did not contain the required config keys:{}!'.format(self.name, DER_id,missing_keys))			
		except:
			LogUtil.exception_handler()


	def check_config_file(self,config_file,config_dict):
		"""Check whether DER config file contains necessary fields."""
		try:
			LogUtil.logger.debug('Checking config file {}.'.format(config_file))				 
			for DER_id,DER_config in config_dict.items():
				self.check_DER_config(DER_config,DER_id)
		except:
			LogUtil.exception_handler()


	def modify_DER_parameters(self,parameter_ID):
		"""Modify the DER parameters to parameters corresponding to the given parameter ID.
		Args:
			 parameter_ID (str): User specified parameter ID (can be None). 
		"""
		try:
			self.parameter_ID = self.get_DER_id(derId=parameter_ID)
			self.initialize_module_parameters()
			self.initialize_DER_model()
		
			self.initialize_states(ia0 = 0+0j, xa0 = 0+0j, ua0 = 0+0j,\
									xDC0 = 0, xQ0 = 0, xPLL0 = 0.0,wte0 = 2*math.pi)
		
			self.initialize_derived_quantities()
		
			LogUtil.logger.info('{}:PV-DER parameters updated with parameters fromparameter dictionary {}!'.format(self.name,self.parameter_ID))
		except:
			LogUtil.exception_handler()

	def initialize_basic_specs(self):
		"""Initialize number of ODEs and phases"""
		try:
			self.n_ODE = templates.DER_design_template[self.DER_model_type]['basic_specs']['n_ODE'] #23#Number of ODE's
			self.n_phases = templates.DER_design_template[self.DER_model_type]['basic_specs']['n_phases'] #3 #Number of phases
		except:
			LogUtil.exception_handler()

	def initialize_basic_options(self):
		"""Initialize basic options"""
		try:
			self.t_stable = self.DER_config['basic_options']['t_stable'] #0.5#Time delay before activating logic for MPP, Volt-VAR control,LVRT/LFRT 
			self.m_steady_state = self.DER_config['basic_options']['m_steady_state'] #0.96 #Expected duty cycle at steady state	
			self.current_gradient_limiter = self.DER_config['basic_options']['current_gradient_limiter'] #0.96 #Expected duty cycle at steady state	
		except:
			LogUtil.exception_handler()

	def initialize_grid_measurements(self,DER_arguments):
		"""Initialize inverter states.
		Args:
			 gridVoltagePhaseA (complex): Value of gridVoltagePhaseA
			 gridVoltagePhaseB (complex): Value of gridVoltagePhaseB
			 gridVoltagePhaseC (complex): Value of gridVoltagePhaseC
			 gridFrequency (float): Value of gridFrequency
		
		"""		
		try:
			if not self.standAlone:
				if 'gridFrequency' in DER_arguments:
					self.gridFrequency = DER_arguments['gridFrequency']
				else:
					raise ValueError('Grid voltage source Frequency need to be supplied if model is not stand alone!')
			
				if templates.DER_design_template[self.DER_model_type]['basic_specs']['n_phases'] >=1: #Check if model has one phase
					if 'gridVoltagePhaseA' in DER_arguments:
						self.gridVoltagePhaseA = DER_arguments['gridVoltagePhaseA']/self.Vbase
					else:
						raise ValueError('Grid voltage source phase A need to be supplied if model is not stand alone!')
					
				if templates.DER_design_template[self.DER_model_type]['basic_specs']['n_phases'] >=2: #Check if model has 2 phases
					if templates.DER_design_template[self.DER_model_type]['basic_specs']['unbalanced']: #Check if model is unbalanced
						if 'gridVoltagePhaseB' in DER_arguments:
							self.gridVoltagePhaseB = DER_arguments['gridVoltagePhaseB']/self.Vbase
						else:
							raise ValueError('Grid voltage source phase B need to be supplied if model is not stand alone!')
					else:
						self.gridVoltagePhaseB = utility_functions.Ub_calc(self.gridVoltagePhaseA)
					
				if templates.DER_design_template[self.DER_model_type]['basic_specs']['n_phases'] >=3: #Check if model has 3 phases
					if templates.DER_design_template[self.DER_model_type]['basic_specs']['unbalanced']: #Check if model is unbalanced
						if 'gridVoltagePhaseC' in DER_arguments:
							self.gridVoltagePhaseC = DER_arguments['gridVoltagePhaseC']/self.Vbase
						else:
							raise ValueError('Grid voltage source phase C need to be supplied if model is not stand alone!')
					else:
						self.gridVoltagePhaseC = utility_functions.Uc_calc(self.gridVoltagePhaseA)
			
				if templates.DER_design_template[self.DER_model_type]['basic_specs']['n_phases'] >=4: #Check if model has 3 phases
					raise ValueError('Model has more than 3 phases!')
		except:
			LogUtil.exception_handler()


	def initialize_states(self,DER_arguments):
		"""Initialize inverter states.

		Args:
			ia0 (float): Initial current
			xa0 (float): Initial controller state
			ua0 (float): Initial controller state

		"""
		try:
			if 'ia' in DER_arguments:
				ia0 = DER_arguments['ia']
			else:
				ia0 = self.DER_config['initial_states']['iaR'] + 1j*self.DER_config['initial_states']['iaI']
			if 'xa' in DER_arguments:
				xa0 = DER_arguments['xa']
			else:
				xa0 = self.DER_config['initial_states']['xaR'] + 1j*self.DER_config['initial_states']['xaI']
			if 'ua' in DER_arguments:
				ua0 = DER_arguments['ua']
			else:
				ua0 = self.DER_config['initial_states']['uaR'] + 1j*self.DER_config['initial_states']['uaI']
		
			if 'xDC' in templates.DER_design_template[self.DER_model_type]['initial_states']:
				if 'xDC' in DER_arguments:
					xDC0 = DER_arguments['xDC']
				else:
					xDC0 = self.DER_config['initial_states']['xDC']
		
			if 'xP' in templates.DER_design_template[self.DER_model_type]['initial_states']:
				if 'xP' in DER_arguments:
					xP0 = DER_arguments['xP']
				else:
					xP0 = self.DER_config['initial_states']['xP']
				
			if 'xQ' in templates.DER_design_template[self.DER_model_type]['initial_states']:
				if 'xQ' in DER_arguments:
					xQ0 = DER_arguments['xQ']
				else:
					xQ0 = self.DER_config['initial_states']['xQ']
		
			xPLL0 = self.DER_config['initial_states']['xPLL']
			wte0 = self.DER_config['initial_states']['wte']		
		
			self.Vdc = self.Vdc_ref#DC link voltage		
			self.Ppv = self.Ppv_calc(self.Vdc_actual) #PV module power output	
		
			 #Initialize all states with steady state values at current operating point
			if self.steady_state_initialization:			
				self.steady_state_calc()
			else:
				#Phase a
				self.ia = ia0
				self.xa = xa0
				self.ua = ua0
			
				if 'xDC' in templates.DER_design_template[self.DER_model_type]['initial_states']:
					self.xDC = xDC0 #DC link voltage controller
				if 'xP' in templates.DER_design_template[self.DER_model_type]['initial_states']:
					self.xP = xP0 #Active power controller 
				if 'xQ' in templates.DER_design_template[self.DER_model_type]['initial_states']:
					self.xQ = xQ0 #Reactive power controller
			
				self.xPLL = xPLL0 #PLL
				self.wte = wte0
			
				if self.DER_model_type in templates.three_phase_models:
					ib0 = utility_functions.Ub_calc(ia0)
					xb0 = utility_functions.Ub_calc(xa0)
					ub0 = utility_functions.Ub_calc(ua0)
				
					ic0 = utility_functions.Uc_calc(ia0)
					xc0 = utility_functions.Uc_calc(xa0)
					uc0 = utility_functions.Uc_calc(ua0)
		
					#Phase b
					self.ib = ib0
					self.xb = xb0#Shift by -120 degrees
					self.ub = ub0

					#Phase c
					self.ic = ic0
					self.xc = xc0 #Shift by +120 degrees
					self.uc = uc0
		except:
			LogUtil.exception_handler()

	def initialize_derived_quantities(self):
		"""Initialize quantities other than states."""
		try:
			self.update_voltages()
			self.update_power()
			self.update_RMS()
			if self.DER_model_type in templates.constant_Vdc_models:
				self.ia_ref  = self.ia_ref_activepower_control()
			else:
				self.ia_ref = self.ia_ref_calc()
			self.t_iref = 0.0
			self.update_iref(t=0.000001) #Reference currents
			self.update_inverter_frequency(t=0.0)
		except:
			LogUtil.exception_handler()

	def initialize_DER_model(self):
		"""Initialize DER ratings.

		Args:
			 DER_arguments (dict): Key word arguments.

		Raises:
			 ValueError: If specified parameters correponding to `parameter_ID` is not available.
		"""
		try:
			LogUtil.logger.log(10,'Creating inverter instance for DER with parameter ID:{}!'.format(self.parameter_ID))
		
			self.initialize_module_parameters()
			self.initialize_inverter_ratings() #Initialize inverter ratings according to DER rating
			self.check_voltage() #Check if voltage is feasible			 
		
			self.initialize_circuit_parameters()
			self.check_circuit_parameters()#Check PV-DER circuit parameters
			
			self.initialize_controller_gains()#Initialize control loop gains according to DER rating		 
		except:
			LogUtil.exception_handler()


	def initialize_inverter_ratings(self):
		"""Initialize inverter voltage and power ratings."""
		try:
			if self.parameter_ID in self.inverter_ratings:
				self.initialize_Sinverter()
				self.initialize_Vdc()
				self.initialize_Vac()
				self.initialize_Iac()
			else:
				raise ValueError('Inverter voltage, current, power ratings not available for parameter ID {}!'.format(self.parameter_ID))
		except:
			LogUtil.exception_handler()


	def initialize_Sinverter(self):
		"""Initialize inverter power rating."""
		try:
			self.Sinverter_rated = self.inverter_ratings[self.parameter_ID]['Srated']
			self.Sinverter_nominal = (self.Sinverter_rated/BaseValues.Sbase) #Converting to p.u. value		 
		except:
			LogUtil.exception_handler()

	def initialize_Vdc(self):
		"""Initialize DC side voltages."""
		try:
			self.Vdcrated = self.inverter_ratings[self.parameter_ID]['Vdcrated'] #Rated DC link voltage in V
			self.Vdcnominal = self.Vdcrated/self.Vdcbase #Converting to p.u. value 
			self.Vdc_ref = self.set_Vdc_ref()
			self.Vdc_ref_new = self.set_Vdc_ref()		
		except:
			LogUtil.exception_handler()

	def initialize_Vac(self):
		"""Initialize AC side voltages."""
		try:
			self.Vrms_rated = self.inverter_ratings[self.parameter_ID]['Vrmsrated']
			self.Varated = self.Vrms_rated*math.sqrt(2) #L-G peak to peak 
			self.Vanominal = self.Varated/BaseValues.Vbase #Converting to p.u. value
			self.Vrms_ref =self.Vanominal/math.sqrt(2) #Reference voltage		
		except:
			LogUtil.exception_handler()

	def initialize_Iac(self):
		"""Initialize AC side currents."""
		try:
			self.Iarated = (self.Sinverter_rated/(self.n_phases*(self.Varated/math.sqrt(2))))*math.sqrt(2)
			self.Ioverload = self.inverter_ratings[self.parameter_ID]['Ioverload']#Inverter current overload rating (Max 10s)			
			self.iref_limit = (self.Iarated/self.Ibase)*self.Ioverload #Maximum current reference
			self.iR_ramp_up_max_gradient = (self.inverter_ratings[self.parameter_ID]['Iramp_max_gradient_real']*self.Iarated)/self.Ibase
			self.iI_ramp_up_max_gradient = (self.inverter_ratings[self.parameter_ID]['Iramp_max_gradient_imag']*self.Iarated)/self.Ibase
		except:
			LogUtil.exception_handler()

	def initialize_circuit_parameters(self):
		"""Initialize C, Lf, and Rf parameters."""
		try:
			if self.parameter_ID in self.inverter_ratings:
				if self.standAlone:
					self.a = self.grid_model.Vgridrated/self.Varated#Transformer turns ratio

				if 'Rf_actual' in templates.DER_design_template[self.DER_model_type]['circuit_parameters']:
					self.Rf_actual =self.circuit_parameters[self.parameter_ID]['Rf_actual'] #Filter resistance
					self.Rf = self.Rf_actual/BaseValues.Zbase		# VSC Filter resistance 

				if 'Lf_actual' in templates.DER_design_template[self.DER_model_type]['circuit_parameters']:
					self.Lf_actual = self.circuit_parameters[self.parameter_ID]['Lf_actual']#Filter inductance
					self.Lf = self.Lf_actual/BaseValues.Lbase		# VSC Filter inductance
					self.Xf = (self.Lf_actual*BaseValues.wbase)/BaseValues.Zbase # VSC Filter reactance			

				if 'C_actual' in templates.DER_design_template[self.DER_model_type]['circuit_parameters']:
					self.C_actual =self.circuit_parameters[self.parameter_ID]['C_actual'] #DC link capacitance
					self.C = self.C_actual/BaseValues.Cbase		#DC link capacitor capacitance
			
				self.Zf_actual =self.Rf_actual + 1j*self.Lf_actual*BaseValues.wbase
				self.Zf = self.Zf_actual/BaseValues.Zbase

				#Interconnection to PCC - HV side
				#Actual values
				self.Z1_actual = self.circuit_parameters[self.parameter_ID]['Z1_actual']
				self.R1_actual = self.Z1_actual.real
				self.L1_actual = self.Z1_actual.imag/(2*math.pi*60.0)

				#Per-unit values
				self.R1 = self.R1_actual/BaseValues.Zbase#Line/transformer resistance
				self.L1 = self.L1_actual/BaseValues.Lbase#Line/transformer inductance
				self.Z1 =self.Z1_actual/BaseValues.Zbase	#Line/transformer impedance

				#Exhibit different behavior when used in standalone mode
				if self.standAlone:
					self.n_total_ODE = self.n_ODE + self.grid_model.n_ODE
					self.Zload1_actual =self.events.load_events(t=0.0)		#Load at PCC 

				else:
					self.n_total_ODE = self.n_ODE
					self.Zload1_actual =10e6+1j*0.0	 #Using a large value of impedance to represent a no load condition

				self.Rload1_actual = self.Zload1_actual.real
				self.Lload1_actual = self.Z1_actual.imag/(2*math.pi*60.0)

				#Per-unit values
				self.Rload1 = self.Rload1_actual/BaseValues.Zbase#resistance in per unit value
				self.Lload1 = self.Lload1_actual/BaseValues.Lbase#inductance in per unit value
				self.Zload1 = self.Zload1_actual/BaseValues.Zbase

				self.transformer_name = 'transformer_'+str(self.ID)			
			
			else:
				raise ValueError('Inverter circuit_parameters not available for parameter ID: {}!'.format(self.parameter_ID))		
		except:
			LogUtil.exception_handler()

	def initialize_controller_gains(self):
		"""Initialize controller settings."""
		try:
			if self.parameter_ID in self.controller_gains:
				#Current controller gains
				self.Kp_GCC = self.controller_gains[self.parameter_ID]['Kp_GCC'] #Current controller Proportional constant
				self.Ki_GCC = self.controller_gains[self.parameter_ID]['Ki_GCC']#Current controller Integral constant
				self.wp =self.controller_gains[self.parameter_ID]['wp']	#First order filter gain for GCC

				#Power (active and reactive) controller gains
				if 'Kp_DC' in templates.DER_design_template[self.DER_model_type]['controller_gains']:
					self.Kp_DC = self.controller_gains[self.parameter_ID]['Kp_DC'] #Active power controller Proportional constant
				if 'Kp_DC' in templates.DER_design_template[self.DER_model_type]['controller_gains']:
					self.Ki_DC = self.controller_gains[self.parameter_ID]['Ki_DC'] #Active power controller Integral constant
			
				if 'Kp_P' in templates.DER_design_template[self.DER_model_type]['controller_gains']:
					self.Kp_P = self.controller_gains[self.parameter_ID]['Kp_P']#Active power controller Proportional constant
				if 'Ki_P' in templates.DER_design_template[self.DER_model_type]['controller_gains']:
					self.Ki_P = self.controller_gains[self.parameter_ID]['Ki_P'] #Active power controller Integral constant 
			
				if 'Kp_Q' in templates.DER_design_template[self.DER_model_type]['controller_gains']:
					self.Kp_Q = self.controller_gains[self.parameter_ID]['Kp_Q']#Reactive power controller Proportional constant
				if 'Ki_Q' in templates.DER_design_template[self.DER_model_type]['controller_gains']:
					self.Ki_Q = self.controller_gains[self.parameter_ID]['Ki_Q'] #Reactive power controller Integral constant 
			else:
				raise ValueError('Controller gains not available for parameter ID: {}!'.format(self.parameter_ID))	
		except:
			LogUtil.exception_handler()

	def initialize_jacobian(self):
		"""Create a Jacobian matrix with zero values."""
		try:
			self.J = np.zeros((self.n_total_ODE,self.n_total_ODE))
			self.varInd={}
			n=0
			if self.DER_model_type == 'SolarPVDERSinglePhase':
				state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
							'Vdc','xDC','xQ','xPLL','wte']
		
			elif self.DER_model_type == 'SolarPVDERThreePhase':
				state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
							'ibR','ibI','xbR','xbI','ubR','ubI',
							'icR','icI','xcR','xcI','ucR','ucI',
							'Vdc','xDC','xQ','xPLL','wte']
		
			elif self.DER_model_type == 'SolarPVDERSinglePhaseConstantVdc':
				state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
							'xP','xQ',
							'xPLL','wte']
		
			elif self.DER_model_type == 'SolarPVDERThreePhaseConstantVdc':
				state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
							'ibR','ibI','xbR','xbI','ubR','ubI',
							'icR','icI','xcR','xcI','ucR','ucI',
							'xP','xQ',
							'xPLL','wte']
			elif self.DER_model_type == 'SolarPVDERThreePhaseBalanced':
				state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
							'Vdc',
							'xDC','xQ',
							'xPLL','wte']
			elif self.DER_model_type == 'SolarPVDERThreePhaseNumba':
				state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
							'ibR','ibI','xbR','xbI','ubR','ubI',
							'icR','icI','xcR','xcI','ucR','ucI',
							'Vdc','xDC','xQ','xPLL','wte']		
			else:
				raise ValueError('{} is not a valid model type! - Valid model types:{}'.format(self.DER_model_type,templates.model_types))
			for entry in state_list:
				self.varInd[entry]=n
				n+=1
		except:
			LogUtil.exception_handler()

	def attach_grid_model(self,DER_arguments):
		"""Attach a grid model to the PV-DER instance.
		Args:
			 grid_model: An instance of `GridModel`.
		Returns:
			 bool: Description of return value
		"""
		try:
			#Connect grid instance only if working in stand alone mode
			if self.standAlone:
				if 'gridModel' in DER_arguments:
					self.grid_model = DER_arguments['gridModel']
					LogUtil.logger.debug('{}:Grid model {} attached since PV-DER instance is stand alone!'.format(self.name,DER_arguments['gridModel'].name))
				else:
					raise ValueError('`Grid` instance need to be provided in stand alone mode for creating `SolarPV_DER` instance`!')
			else: #Grid model is not connected
				LogUtil.logger.debug('{}:No grid model attached since PV-DER instance is not stand alone!'.format(self.name))
		except:
			LogUtil.exception_handler()

	def check_voltage(self):
		"""Method to check whether inverter voltage ratings are feasible.
		Raises:
			 ValueError: If any of the specificied voltage ratings is infeasible.
		"""
		try:
			if not self.standAlone and abs(abs(self.gridVoltagePhaseA) - self.Vanominal)/self.Vanominal > 0.1:
				raise ValueError('The rated PCC-LV voltage {:.2f} V has more than 10% deviation from the voltage input from external program {:.2f} V!'.format(self.Vanominal, abs(self.gridVoltagePhaseA)))
								 
			if self.m_steady_state*(self.Vdcrated/2) < self.Varated:
				raise ValueError('The nominal DC link voltage {:.1f} V is not sufficient for the inverter to generate the nominal voltage at PCC - LV side {:.1f} V (L-G peak). Increase nominal DC link voltage to {:.1f} V.'.format(self.Vdcrated,self.Varated,math.ceil((self.Varated/self.m_steady_state)*2)))
		
			if self.m_steady_state*(self.Vdcmpp_min/2) < self.Varated:
				raise ValueError('The minimum DC link voltage {:.1f} V is not sufficient for the inverter to generate the nominal voltage at PCC - LV side {:.1f} V (L-G peak). Increase minimum DC link voltage to {:.1f} V.'.format(self.Vdcmpp_min,self.Varated,math.ceil((self.Varated/self.m_steady_state)*2)))
		except:
			LogUtil.exception_handler()

	def check_circuit_parameters(self):
		"""Method to check whether inverter circuit parameter's are feasible."""
		try:
			del_I1max = 0.1*self.Iarated
			Lf_min = self.Vdcrated/(16*self.fswitching*del_I1max)
			del_I1max_actual = self.Vdcrated/(16*self.fswitching*self.Lf_actual)
			if del_I1max_actual > del_I1max: #Check if ripple current is less than 10 %
				LogUtil.logger.debug('{}:Filter inductance {:.4f} H is acceptable since AC side current ripple is {:.2f}% (< 10%)'.format(self.name,Lf_min,del_I1max_actual/self.Iarated))
			else:
				LogUtil.logger.debug('{}:Warning:Filter inductance {:.4f} H results in AC side current ripple of {:.2f}% (> 10%)'.format(self.name,Lf_min,del_I1max_actual/self.Iarated))

			if hasattr(self,'C_actual'):
				I_ripple = (0.25*self.Vdcrated)/(self.Lf_actual*self.fswitching)#Maximum ripple voltage (p-p) at DC link
				V_ripple = self.Vdcrated/(32*self.Lf_actual*self.C_actual*(self.fswitching**2))#Maximum ripple voltage (p-p) at DC link
				V_ripple_percentage = (V_ripple/self.Vdcrated)*100
				if V_ripple_percentage <= 1.0: #Check if voltage ripple on DC link is less than 1%
				 LogUtil.logger.debug('{}:DC link capacitance of {:.4f} F is acceptable since voltage ripple is only {:.2f}% (< 1%)'.format(self.name,self.C_actual,V_ripple_percentage))
				else:
				 V_ripple_ideal = self.Vdcrated*0.01#1% ripple is acceptable
				 C = self.Vdcrated/(32*self.Lf_actual*V_ripple_ideal*(self.fswitching**2))
				 LogUtil.logger.debug('{}:Warning:DC link capacitance of {:.4f} F results in DC link voltage ripple of {:.2}% (> 1%)!Please use at least {} F.'.format(self.name,self.C_actual,V_ripple_percentage,C))
				 #warnings.warn('Warning:DC link capacitance of {} F results in DC link voltage ripple of {:.3}% (> 1%)!Please use at least {} F.'.format(self.C_actual,_V_ripple_percentage,_C))
		except:
			LogUtil.exception_handler()

	def power_error_calc(self,x):
		"""Function for power."""
		try:
			maR = x[0]
			maI = x[1]
			iaR = x[2]
			iaI = x[3]
			ma = maR + 1j*maI
			ia = self.ia = iaR + 1j*iaI
			vta = ma*(self.Vdc/2)
			va = self.va_calc()
			St = (vta*ia.conjugate())/2
			S_PCC = (va*ia.conjugate())/2

			#Sloss_filter = ((abs(ia)/math.sqrt(2))**2)*self.Zf 
			if self.DER_model_type in templates.three_phase_models:
				if not self.allow_unbalanced_m:
					mb = utility_functions.Ub_calc(ma)
					mc = utility_functions.Uc_calc(ma)
					ib = self.ib = utility_functions.Ub_calc(ia)
					ic = self.ic = utility_functions.Uc_calc(ia)
				elif self.allow_unbalanced_m:
					mbR = x[4]
					mbI = x[5]
					ibR = x[6]
					ibI = x[7]

					mcR = x[8]
					mcI = x[9]
					icR = x[10]
					icI = x[11]

					mb = mbR + 1j*mbI
					mc = mcR + 1j*mcI

					ib = self.ib = ibR + 1j*ibI
					ic = self.ic = icR + 1j*icI

				vtb = mb*(self.Vdc/2)
				vtc = mc*(self.Vdc/2)

				vb = self.vb_calc()
				vc = self.vc_calc()

				St = St + (vtb*ib.conjugate() + vtc*ic.conjugate())/2
				S_PCC = S_PCC + (vb*ib.conjugate() + vc*ic.conjugate())/2
			
				#Sloss_filter = Sloss_filter + ((abs(ib)/math.sqrt(2))**2 + (abs(ic)/math.sqrt(2))**2)*self.Zf		 
		
			P_error = (St.real - self.Ppv)**2
			Q_error = (S_PCC.imag - self.Q_ref)**2
		
			#Qloss_filter_expected = self.n_phases*((self.Ppv/(self.n_phases*self.Vrms_ref))**2)*self.Xf
			#Ploss_filter_expected = self.n_phases*((self.Ppv/(self.n_phases*self.Vrms_ref))**2)*self.Rf
				
			#P_PCC_error = ((S_PCC.real + Sloss_filter.real) - self.Ppv)**2		 
			#S_error = (abs(St -(S_PCC+Sloss_filter)))**2
			#Q_error_filter_expected =(St.imag - Qloss_filter_expected)**2 
				
			if self.DER_model_type in templates.three_phase_models:
				del_1 = utility_functions.relative_phase_calc(ma,mb)
				del_2 = utility_functions.relative_phase_calc(ma,mc)
				del_3 = utility_functions.relative_phase_calc(mb,mc)
			
				if del_1 > math.pi:
					del_1 = abs(del_1 - 2*math.pi)

				if del_2 > math.pi:
					del_2 = abs(del_2 - 2*math.pi)

				if del_3> math.pi:
					del_3 = abs(del_3 - 2*math.pi)
			
				if self.allow_unbalanced_m:
					m_error = (abs((va+vb+vc)- (vta+vtb+vtc)))**2 
				else:
					m_error = (del_1 - PHASE_DIFFERENCE_120)**2 + (del_2 - PHASE_DIFFERENCE_120)**2 + (del_3 - PHASE_DIFFERENCE_120)**2 +\
							(abs(ma) - abs(mb))**2 + (abs(ma) - abs(mc))**2 + (abs(mb) - abs(mc))**2 
			
				m_error = m_error
							
				i_error = (abs(ia + ib+ic) +abs(vta-va - ia*self.Zf)+abs(vtb-vb - ib*self.Zf)+abs(vtc-vc - ic*self.Zf))**2
			else:
				m_error = 0.0
				i_error = (abs(vta-va - ia*self.Zf))**2

			return P_error + Q_error + m_error + i_error #+S_error + P_PCC_error
		except:
			LogUtil.exception_handler()

	def steady_state_calc(self):
		"""Find duty cycle and inverter current that minimize steady state error and return steady state values."""		
		try:
			LogUtil.logger.debug('Solving for steady state at current operating point.') 
			if self.standAlone:
				x0 = [self.steadystate_values[self.parameter_ID]['maR0'],self.steadystate_values[self.parameter_ID]['maI0'],
					self.steadystate_values[self.parameter_ID]['iaR0'],self.steadystate_values[self.parameter_ID]['iaI0']]		
		
			else:
				va = self.va_calc()
				x0 = [va.real,va.imag,
					self.steadystate_values[self.parameter_ID]['iaR0'],self.steadystate_values[self.parameter_ID]['iaI0']]	
		
			if self.allow_unbalanced_m:
				LogUtil.logger.info('Using unbalanced option - steady state duty cycles may be unbalanced.') 
				if self.standAlone:
					x0.extend([self.steadystate_values[self.parameter_ID]['mbR0'],self.steadystate_values[self.parameter_ID]['mbI0'],
							 self.steadystate_values[self.parameter_ID]['ibR0'],self.steadystate_values[self.parameter_ID]['ibI0'],
							 self.steadystate_values[self.parameter_ID]['mcR0'],self.steadystate_values[self.parameter_ID]['mcI0'],
							 self.steadystate_values[self.parameter_ID]['icR0'],self.steadystate_values[self.parameter_ID]['icI0']])	 
				else:
					vb = self.vb_calc()
					vc = self.vc_calc()
			
					x0.extend([vb.real,vb.imag,
							 self.steadystate_values[self.parameter_ID]['ibR0'],self.steadystate_values[self.parameter_ID]['ibI0'],
							 vc.real,vc.imag,
							 self.steadystate_values[self.parameter_ID]['icR0'],self.steadystate_values[self.parameter_ID]['icI0']])	
			
			x0 = np.array(x0)
		
			self.solver_spec[self.steadystate_solver].update({'disp':bool(self.verbosity == 'DEBUG')})				 
			result = minimize(self.power_error_calc, x0,
							method=self.steadystate_solver,options=self.solver_spec[self.steadystate_solver])
		
			if not result.success:
				raise ValueError('Steady state solution did not converge! Change operating point or disable steady state flag and try again.')
			print(f"Rf:{self.Rf},Optimization results:{result.x}")
			if 'xDC' in templates.DER_design_template[self.DER_model_type]['initial_states']:
			 self.xDC = self.ia.real
			if 'xP' in templates.DER_design_template[self.DER_model_type]['initial_states']:
			 self.xP = self.ia.real
			if 'xQ' in templates.DER_design_template[self.DER_model_type]['initial_states']:
			 self.xQ = self.ia.imag
			if 'xPLL' in templates.DER_design_template[self.DER_model_type]['initial_states']:
			 self.xPLL = 0.0
			if 'wte' in templates.DER_design_template[self.DER_model_type]['initial_states']:
			 self.wte = 2*math.pi
		 
			if self.DER_model_type in templates.single_phase_models or self.DER_model_type in templates.three_phase_models:
				ma0 = result.x[0] + 1j*result.x[1]
				self.ia = result.x[2] + 1j*result.x[3]
				self.xa = ma0
				self.ua = 0.0+0.0j			
				self.vta = self.vta_calc()
				self.va = self.va_calc()
		
			if self.DER_model_type in templates.three_phase_models:
			
				if not self.allow_unbalanced_m:
					mb0 = utility_functions.Ub_calc(ma0)
					mc0 = utility_functions.Uc_calc(ma0)
				
					self.ib = utility_functions.Ub_calc(self.ia)
					self.ic = utility_functions.Uc_calc(self.ia)			
			
				else:
					mb0 = result.x[4] + 1j*result.x[5]
					mc0 = result.x[8] + 1j*result.x[9]
					self.ib = result.x[6] + 1j*result.x[7]
					self.ic = result.x[10] + 1j*result.x[11]				
			
				self.xb = mb0
				self.xc = mc0
			
				self.ub = 0.0+0.0j
				self.uc = 0.0+0.0j			
			
				self.vtb = self.vtb_calc()
				self.vtc = self.vtc_calc()
			
				self.vb = self.vb_calc()
				self.vc = self.vc_calc()	 
			
			if self.DER_model_type == 'SolarPVDERThreePhaseNumba':
				from pvder import utility_functions_numba #Import numba lazily
				
				self.S = utility_functions_numba.S_calc(self.vta,self.vtb,self.vtc,self.ia,self.ib,self.ic)
				self.S_PCC = utility_functions_numba.S_calc(self.va,self.vb,self.vc,self.ia,self.ib,self.ic)
				self.Vtrms = utility_functions_numba.Urms_calc(self.vta,self.vtb,self.vtc)
				self.Vrms = utility_functions_numba.Urms_calc(self.va,self.vb,self.vc)
				self.Irms = utility_functions_numba.Urms_calc(self.ia,self.ib,self.ic)
			else:
				self.S =self.S_calc()
				self.S_PCC = self.S_PCC_calc()
				self.Vtrms = self.Vtrms_calc()
				self.Vrms = self.Vrms_calc()
				self.Irms = self.Irms_calc()
		
			LogUtil.logger.debug('{}:Steady state values for operating point defined by Ppv:{:.2f} W, Vdc:{:.2f} V, va:{:.2f} V found at:'.format(self.name,self.Ppv*self.Sbase,self.Vdc*self.Vdcbase,self.va*self.Vbase))
			
			if self.verbosity == 'DEBUG':
				self.show_PV_DER_states(quantity='power')
				self.show_PV_DER_states(quantity='duty cycle')
				self.show_PV_DER_states(quantity='voltage')
				self.show_PV_DER_states(quantity='current')
		except:
			LogUtil.exception_handler()

	def check_jacobian(self,t=0.0):
		"""Compare analytical and numerical Jacobian of the ODE model."""
		try:
			#Calculate numerical Jacobian using finite differences
			x = self.y0
			x0 = x.copy()

			eps = 1e-6
			Jn = np.zeros([len(x), len(x)], dtype = np.float)

			for i in range(len(x)):
				x1 = x.copy()
				x2 = x.copy()

				x1[i] += eps
				x2[i] -= eps

				f1 = self.ODE_model(x1,t)
				f2 = self.ODE_model(x2,t)

				Jn[ : , i] = (f1 - f2) / (2 * eps)
		
			Ja = self.jac_ODE_model(x0,t) #Calculate analytical Jacobian
		
			#Compare numerical Jacobian with analytical Jacobian
			err=np.abs(Ja-Jn)
			err[err<1e-2]=0
			ind=np.where(err>0)

			nErr=len(ind[0])
			err=[]
		
			if nErr > 0:
				print('Differences in analytical and numerical Jacobian found in {} entries!'.format(nErr))
				jac_CHECK = False
				for r,c in zip(ind[0],ind[1]):
					err.append([r,c,Jn[r,c],Ja[r,c]])
					#print(r,c,Jn[r,c],Ja[r,c])
					r_char = list(self.varInd.keys())[list(self.varInd.values()).index(r)]
					c_char = list(self.varInd.keys())[list(self.varInd.values()).index(c)]
					print('J[{}][{}]--Jn:{},Ja:{}'.format(r_char,c_char,Jn[r,c],Ja[r,c]))
			else:
				jac_CHECK = True
				print('No differences in analytical and numerical Jacobian!')
		
			return jac_CHECK,Jn,Ja
		except:
			LogUtil.exception_handler()


