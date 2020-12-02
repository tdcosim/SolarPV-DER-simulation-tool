from __future__ import division
import sys
import os
import argparse
import logging
import unittest

import math

import matplotlib.pyplot as plt

from pvder.DER_components_three_phase import SolarPVDERThreePhase
from pvder.grid_components import Grid
from pvder.dynamic_simulation import DynamicSimulation
from pvder.simulation_events import SimulationEvents
from pvder.simulation_utilities import SimulationResults

from unittest_utilities import show_DER_status, plot_DER_trajectories
config_file = r'..\config_der.json'

def suite():
	"""Define a test suite."""
	
	all_tests = ['test_init','test_run_simulation']
	
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
	events = SimulationEvents()
	
	flag_arguments = {'standAlone': False,
					  'steadyStateInitialization':True,
					  'verbosity':'DEBUG'}
	ratings_arguments ={'powerRating':power_rating,
						'VrmsRating':Vrms}
	voltage_arguments = {'gridVoltagePhaseA': Va,
						 'gridVoltagePhaseB':Vb,
						 'gridVoltagePhaseC':Vc,
						 'gridFrequency':wgrid}						 
	
	def test_init(self):
		"""Test Dynamic Simulation initialization."""
		events = SimulationEvents()
		kwargs={}
		kwargs.update(self.flag_arguments)
		kwargs.update(self.ratings_arguments)
		kwargs.update(self.voltage_arguments)
		PVDER = SolarPVDERThreePhase(events = events,configFile=config_file,**kwargs)

		sim = DynamicSimulation(PV_model=PVDER,events = events,jacFlag = True,verbosity = 'DEBUG',solverType='odeint')
		self.assertIsInstance(sim, DynamicSimulation)
		self.assertTrue(sim.jacFlag)
		
	def test_run_simulation(self):
		"""Test run simulation method.""" 
		
		events = SimulationEvents()
		kwargs={}
		kwargs.update(self.flag_arguments)
		kwargs.update(self.ratings_arguments)
		kwargs.update(self.voltage_arguments)
		PVDER = SolarPVDERThreePhase(events = events,configFile=config_file,**kwargs)

		sim = DynamicSimulation(PV_model=PVDER,events = events,
								jacFlag = True,verbosity = 'DEBUG',solverType='odeint')
		sim.tStop = 10.0
		sim.tInc = 1/120.
		sim.run_simulation()

		self.assertEqual(sim.t[-1],sim.tStop)
		self.assertTrue(sim.SOLVER_CONVERGENCE)  


if __name__ == '__main__':
	#unittest.main()
	logging.debug('test')
	runner = unittest.TextTestRunner()
	runner.run(suite())
