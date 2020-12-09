from __future__ import division
import sys
import os
import argparse
import unittest

import math

import matplotlib.pyplot as plt

from pvder.DER_components_single_phase import SolarPVDERSinglePhase
from pvder.grid_components import Grid
from pvder.dynamic_simulation import DynamicSimulation
from pvder.simulation_events import SimulationEvents
from pvder.simulation_utilities import SimulationResults
from unittest_utilities import show_DER_status, plot_DER_trajectories
config_file = r'..\config_der.json'


def suite():
	"""Define a test suite."""
	all_tests = ['test_init','test_parameter_dict','test_jacobian']
	
	avoid_tests = ['test_jacobian']
   
	tests = list(set(all_tests) - set(avoid_tests))
	print('Following unittest scenarios will be run:{}'.format(tests))
	suite = unittest.TestSuite()
	
	for test in tests:
		suite.addTest(TestPVDER(test))
											   
	return suite

class TestPVDER(unittest.TestCase):
	
	Vnominal = 1.0
	
	scaler = 0.835
	Va = (.50+0j)*Grid.Vbase*scaler
	Vb = (-.25-.43301270j)*Grid.Vbase*scaler
	Vc = (-.25+.43301270j)*Grid.Vbase*scaler
	
	power_rating = 10.0e3
	Vrms = abs(Va)/math.sqrt(2)
	
	wgrid = 2*math.pi*60.0
	
	flag_arguments = {'standAlone': False,
					  'steadyStateInitialization':True,
					  'verbosity':'DEBUG'}
	ratings_arguments ={'powerRating':power_rating,
						'VrmsRating':Vrms}
	voltage_arguments = {'gridVoltagePhaseA':Va,
						 'gridFrequency':wgrid}		
	
	def test_init(self):
		"""Test PV-DER three phase mode."""		  
						
		events = SimulationEvents()
		kwargs={}
		kwargs.update(self.flag_arguments)
		kwargs.update(self.ratings_arguments)
		kwargs.update(self.voltage_arguments)
		PVDER = SolarPVDERSinglePhase(events = events,configFile=config_file,**kwargs)
	
		self.assertIsInstance(PVDER, SolarPVDERSinglePhase)
		self.assertTrue(PVDER.steady_state_initialization)	  
		
	def test_parameter_dict(self):
		"""Test initalization and update of paraemter dictionary.""" 
		
		source_ID = '10'
		new_ID ='test_DER'
		new_module_parameters = {'Np':5}
		new_circuit_parameters = {'C_actual':100.0e-6}
				
		events = SimulationEvents()
		kwargs={}
		kwargs.update(self.flag_arguments)
		kwargs.update(self.ratings_arguments)
		kwargs.update(self.voltage_arguments)
		PVDER = SolarPVDERSinglePhase(events = events,configFile=config_file,**kwargs)
	
		PVDER.initialize_parameter_dict(parameter_ID = new_ID,source_parameter_ID = source_ID)
		
		self.assertEqual(PVDER.module_parameters[source_ID]['Np'],PVDER.module_parameters[new_ID]['Np'])
		self.assertEqual(PVDER.inverter_ratings[source_ID]['Srated'],PVDER.inverter_ratings[new_ID]['Srated'])
		self.assertEqual(PVDER.circuit_parameters[source_ID]['Rf_actual'],PVDER.circuit_parameters[new_ID]['Rf_actual'])
		
		PVDER.update_parameter_dict(parameter_ID = new_ID,parameter_type = 'module_parameters',parameter_dict = new_module_parameters)
		PVDER.update_parameter_dict(parameter_ID = new_ID,parameter_type = 'circuit_parameters',parameter_dict = new_circuit_parameters)		
						 
		self.assertEqual(PVDER.module_parameters[new_ID]['Np'],new_module_parameters['Np']) 
		self.assertEqual(PVDER.circuit_parameters[new_ID]['C_actual'],new_circuit_parameters['C_actual'])		
	
	def test_jacobian(self):
		"""Test PV-DER Jacobian."""		  
						
		events = SimulationEvents()
		kwargs={}
		kwargs.update(self.flag_arguments)
		kwargs.update(self.ratings_arguments)
		kwargs.update(self.voltage_arguments)
		PVDER = SolarPVDERSinglePhase(events = events,configFile=config_file,**kwargs)
		jac_CHECK,Jn,Ja = PVDER.check_jacobian()
		self.assertTrue(jac_CHECK,'Analytical and numerical Jacobian should be same.')
		self.assertEqual(Jn.shape,(PVDER.n_ODE,PVDER.n_ODE))	
	
if __name__ == '__main__':
	runner = unittest.TextTestRunner()
	runner.run(suite())
