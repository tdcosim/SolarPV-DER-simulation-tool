from __future__ import division
import sys
import os
import argparse
import logging
import unittest

import math
import cmath

import matplotlib.pyplot as plt

from pvder.DER_components_three_phase import SolarPVDERThreePhase
from pvder.grid_components import Grid
from pvder.dynamic_simulation import DynamicSimulation
from pvder.simulation_events import SimulationEvents
from pvder.simulation_utilities import SimulationResults
from pvder import utility_functions
from unittest_utilities import show_DER_status, plot_DER_trajectories

config_file = r'..\config_der.json'


def suite():
	"""Define a test suite."""
	all_tests = ['test_init','test_parameter_dict','test_jacobian','test_steady_state_calc']
	avoid_tests = []
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
	
	power_rating = 50.0e3
	Vrms = abs(Va)/math.sqrt(2)
	
	wgrid = 2*math.pi*60.0
	
	flag_arguments = {'standAlone': False,
					  'steadyStateInitialization':True,
					  'verbosity':'DEBUG'}
	ratings_arguments ={'powerRating':power_rating,
						'VrmsRating':Vrms}
	voltage_arguments = {'gridVoltagePhaseA': Va,
						 'gridVoltagePhaseB':Vb,
						 'gridVoltagePhaseC':Vc,
						 'gridFrequency':wgrid}			  
	kwargs={}
	kwargs.update(flag_arguments)
	kwargs.update(ratings_arguments)
	kwargs.update(voltage_arguments)

	def test_init(self):
		"""Test PV-DER three phase mode."""		  
						
		events = SimulationEvents()
				
		PVDER = SolarPVDERThreePhase(events = events,configFile=config_file,**self.kwargs)
	
		self.assertIsInstance(PVDER, SolarPVDERThreePhase)
		self.assertTrue(PVDER.steady_state_initialization)	  
		
	def test_parameter_dict(self):
		"""Test initalization and update of paraemter dictionary.""" 
		
		source_ID = '50'
		new_ID ='test_DER'
		new_module_parameters = {'Np':5}
		new_circuit_parameters = {'C_actual':500.0e-6}
				
		events = SimulationEvents()				
				
		PVDER = SolarPVDERThreePhase(events = events,configFile=config_file,**self.kwargs)
	
		PVDER.initialize_parameter_dict(parameter_ID = new_ID,source_parameter_ID = source_ID)
		
		self.assertEqual(PVDER.module_parameters[source_ID]['Np'],PVDER.module_parameters[new_ID]['Np'])
		self.assertEqual(PVDER.inverter_ratings[source_ID]['Srated'],PVDER.inverter_ratings[new_ID]['Srated'])
		self.assertEqual(PVDER.circuit_parameters[source_ID]['Rf_actual'],PVDER.circuit_parameters[new_ID]['Rf_actual'])
		
		PVDER.update_parameter_dict(parameter_ID = new_ID,parameter_type = 'module_parameters',parameter_dict = new_module_parameters)
		PVDER.update_parameter_dict(parameter_ID = new_ID,parameter_type = 'circuit_parameters',parameter_dict = new_circuit_parameters)		
						 
		self.assertEqual(PVDER.module_parameters[new_ID]['Np'],new_module_parameters['Np']) 
		self.assertEqual(PVDER.circuit_parameters[new_ID]['C_actual'],new_circuit_parameters['C_actual'])			
	
	def test_jacobian(self):
		"""Test PV-DER three phase mode."""		  
						
		events = SimulationEvents()				
				
		PVDER = SolarPVDERThreePhase(events = events,configFile=config_file,**self.kwargs)
	
		jac_CHECK,Jn,Ja = PVDER.check_jacobian()
		self.assertTrue(jac_CHECK,'Analytical and numerical Jacobian should be same.')	  
		self.assertEqual(Jn.shape,(PVDER.n_ODE,PVDER.n_ODE))	
	
	def test_steady_state_calc(self):
		"""Test PV-DER three phase mode."""		  
						
		events = SimulationEvents()
		
		voltage_list = [(self.Va,self.Vb,self.Vc),
					(cmath.rect(206.852, math.radians(-36.9906)),cmath.rect(206.128, math.radians(-157.745)),cmath.rect(208.387, math.radians(82.7291))),(169.18+118.52j,utility_functions.Ub_calc(169.18+118.52j),utility_functions.Uc_calc(169.18+118.52j))]	
					
		for voltages in voltage_list:
			Va = voltages[0]
			Vb = voltages[1]
			Vc = voltages[2]
			
			print('Testing voltages:{}'.format(voltages))
					
			PVDER = SolarPVDERThreePhase(events = events,configFile=config_file,**self.kwargs)	
		
			self.assertAlmostEqual(PVDER.Ppv,PVDER.S.real,delta=0.001,msg='Inverter power output must be equal to PV module power output at steady-state!')
			self.assertAlmostEqual(PVDER.S_PCC.imag,PVDER.Q_ref,delta=0.001,msg='Inverter reactive power output must be equal to Q reference!')

			self.assertAlmostEqual(PVDER.Vdc,PVDER.Vdc_ref,delta=0.001,msg='DC link voltage should be equal to reference!')

			self.assertAlmostEqual(PVDER.ma+PVDER.mb+PVDER.mc,0.0+1j*0.0,delta=0.001,msg='Duty cycles should sum to zero!')
			self.assertLess(abs(PVDER.ma),1.0,msg='Magnitude of duty cycle should be less than 1!')

if __name__ == '__main__':
	runner = unittest.TextTestRunner()
	runner.run(suite())
