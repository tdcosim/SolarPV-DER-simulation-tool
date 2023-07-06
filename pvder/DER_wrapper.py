# -*- coding: utf-8 -*-
"""
Created on Thu May 14 13:10:31 2020

@author: splathottam
"""

from pvder.DER_components import SolarPVDER
from pvder.DER_components_three_phase  import SolarPVDERThreePhase
from pvder.DER_components_three_phase_constant_Vdc import SolarPVDERThreePhaseConstantVdc
from pvder.DER_components_three_phase_no_Vrms_filter import SolarPVDERThreePhaseNoVrmsFilter
from pvder.DER_components_three_phase_balanced  import SolarPVDERThreePhaseBalanced
from pvder.DER_components_single_phase import SolarPVDERSinglePhase
from pvder.DER_components_single_phase_constant_Vdc import SolarPVDERSinglePhaseConstantVdc
from pvder import defaults,templates,specifications
from pvder import utility_functions
from pvder.logutil import LogUtil


class DERModel(SolarPVDER):
	"""
	Class providing a wrapper to all the DER models.
	"""

	#def __init__(self,modelType,events,configFile,**kwargs):
	def __init__(self,events,configFile,derID,createDERModel=True,**kwargs):
		"""Creates an instance of `SolarPV_DER_SinglePhase`.
		
		Args:
		  modelType (str): Name of the DER model type.
		  events (SimulationEvents): An instance of `SimulationEvents`.		  
		  gridModel (Grid): An instance of `Grid`(only need to be suppled for stand alone simulation).		  
		  gridVoltatePhaseA,gridVoltatePhaseA,gridVoltatePhaseA (float): Initial voltage phasor (V) at PCC - LV side from external program (only need to be suppled if model is not stand alone).
		  standAlone (bool): Specify if the DER instance is a stand alone simulation or part of a larger simulation.
		  steadyStateInitialization (bool): Specify whether states in the DER instance will be initialized to steady state values.
		  allowUnbalancedM (bool): Allow duty cycles to take on unbalanced values during initialization (default: False).
		  derConfig (dict): Configuration parameters that may be supplied from an external program.
		  identifier (str): An identifier that can be used to name the instance (default: None).
		  
		Raises:
		  ValueError: If parameters corresponding to `Sinverter_rated` are not available.
		  ValueError: If rated DC link voltage is not sufficient.
		
		"""
		try:
			"""
			if modelType == 'ThreePhaseUnbalanced':
				self.DER_model = SolarPVDERThreePhase(events,configFile,**kwargs)
			elif modelType == 'ThreePhaseUnbalancedNoVrmsFilter':
				self.DER_model = SolarPVDERThreePhaseNoVrmsFilter(events,configFile,**kwargs)
			elif modelType == 'ThreePhaseUnbalancedConstantVdc':
				self.DER_model = SolarPVDERThreePhaseConstantVdc(events,configFile,**kwargs)
			elif modelType == 'ThreePhaseBalanced':
				self.DER_model = SolarPVDERThreePhaseBalanced(events,configFile,**kwargs)
			elif modelType == 'SinglePhase':   
				self.DER_model = SolarPVDERSinglePhase(events,configFile,**kwargs)
			elif modelType == 'SinglePhaseConstantVdc':   
				self.DER_model = SolarPVDERSinglePhaseConstantVdc(events,configFile,**kwargs)
			elif modelType == 'ThreePhaseUnbalancedNumba':
				from pvder.DER_components_three_phase_numba  import SolarPVDERThreePhaseNumba
				self.DER_model = SolarPVDERThreePhaseNumba(events,configFile,**kwargs)
			else:
				raise ValueError('{} is not a valid model type! - Valid model types:{}'.format(modelType,templates.model_types))
			"""
			if createDERModel:
				DER_config = self.get_config(configFile,derID)
				modelType = DER_config["basic_specs"]["model_type"] #self.get_DER_model_type(DER_config,DER_parent_config)
				self.create_DER_model(events,modelType,configFile,derID,**kwargs)
						
		except:
			LogUtil.exception_handler()

	def get_config(self,configFile,derID):
		"""Return a config file.		
		
		Args:
			 configFile (str): The inverter power rating in kVA.
			 derId (str): User specified parameter ID (can be None). 

		Returns:
			 dict: Config file
		"""	   
		try:
			DER_config,config_dict = self.get_DER_config(configFile,derID)
			DER_parent_ID = self.get_DER_parent_id(DER_config)
			DER_parent_config = self.get_DER_parent_config(configFile,DER_parent_ID)
			#DER_config, DER_parent_config = self.get_completed_DER_config(DER_config,DER_parent_config)
			DER_config = self.update_DER_config_specs(DER_config,DER_parent_config)
			return DER_config
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
			assert 'derId' in kwargs, 'DER parameter ID was not provided' # nor DER power rating was provided!' # or 'powerRating' in kwargs,
			#if 'derId' in kwargs:
			DER_id = kwargs['derId']   
			#elif 'powerRating' in kwargs:
			#	DER_id = str(int(kwargs['powerRating']/1e3)) 
			return DER_id
		except:
			LogUtil.exception_handler()

	def create_DER_model(self,events,modelType,configFile,derID,**kwargs):
		"""Create DER model object."""
		try:
			if modelType == 'SolarPVDERThreePhase': #'ThreePhaseUnbalanced':
				self.DER_model = SolarPVDERThreePhase(events,configFile,derID,**kwargs)
			elif modelType == 'SolarPVDERThreePhaseNoVrmsFilter': #'ThreePhaseUnbalancedNoVrmsFilter':
				self.DER_model = SolarPVDERThreePhaseNoVrmsFilter(events,configFile,derID,config_dict,**kwargs)
			elif modelType == 'SolarPVDERThreePhaseConstantVdc': #'ThreePhaseUnbalancedConstantVdc':
				self.DER_model = SolarPVDERThreePhaseConstantVdc(events,configFile,derID,config_dict,**kwargs)
			elif modelType == 'SolarPVDERThreePhaseBalanced': #'ThreePhaseBalanced':
				self.DER_model = SolarPVDERThreePhaseBalanced(events,configFile,derID,**kwargs)
			elif modelType == 'SolarPVDERSinglePhase': #'SinglePhase':   
				self.DER_model = SolarPVDERSinglePhase(events,configFile,derID,**kwargs)
			elif modelType == 'SolarPVDERSinglePhaseConstantVdc': #SinglePhaseConstantVdc':   
					self.DER_model = SolarPVDERSinglePhaseConstantVdc(events,configFile,derID,**kwargs)
			elif modelType == 'SolarPVDERThreePhaseNumba': #'ThreePhaseUnbalancedNumba':
				from pvder.DER_components_three_phase_numba  import SolarPVDERThreePhaseNumba
				self.DER_model = SolarPVDERThreePhaseNumba(events,configFile,derID,**kwargs)
			else:
				raise ValueError('{} is not a valid model type! - Valid model types:{}'.format(modelType,list(templates.DER_design_template.keys()))) #templates.model_types
		except:
			LogUtil.exception_handler()

	