# -*- coding: utf-8 -*-
"""
Created on Thu May 14 13:10:31 2020

@author: splathottam
"""

from pvder.DER_components_three_phase  import SolarPVDERThreePhase
from pvder.DER_components_three_phase_constant_Vdc import SolarPVDERThreePhaseConstantVdc
from pvder.DER_components_three_phase_balanced  import SolarPVDERThreePhaseBalanced
from pvder.DER_components_single_phase import SolarPVDERSinglePhase
from pvder.DER_components_single_phase_constant_Vdc import SolarPVDERSinglePhaseConstantVdc
from pvder import defaults,templates,specifications
from pvder.logutil import LogUtil


class DERModel(object):
	"""
	Class providing a wrapper to all the DER models.
	"""

	def __init__(self,modelType,events,configFile,**kwargs): 
		"""Creates an instance of `SolarPV_DER_SinglePhase`.
		
		Args:
		  modelType (str): Name of the DER model type.
		  events (SimulationEvents): An instance of `SimulationEvents`.		  
		  gridModel (Grid): An instance of `Grid`(only need to be suppled for stand alone simulation).
		  powerRating (float): A scalar specifying the rated power (VA) of the DER.
		  VrmsRating (float): A scalar specifying the rated RMS L-G voltage (V) of the DER.
		  ia0,xa0,ua0 (complex): Initial value of inverter states in p.u. of the DER instance.
		  xDC0,xQ0,xPLL0,wte0 (float): Initial value of inverter states in the DER instance.
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
			if modelType == 'ThreePhaseUnbalanced':			
				self.DER_model = SolarPVDERThreePhase(events,configFile,**kwargs)		
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
		except:
			LogUtil.exception_handler()


