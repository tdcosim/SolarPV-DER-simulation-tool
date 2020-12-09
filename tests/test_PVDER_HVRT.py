from __future__ import division
import sys
import os
import math
import argparse

import matplotlib.pyplot as plt
import unittest

from pvder.DER_components_single_phase import SolarPVDERSinglePhase
from pvder.DER_components_three_phase import SolarPVDERThreePhase
from pvder.grid_components import Grid
from pvder.dynamic_simulation import DynamicSimulation
from pvder.simulation_events import SimulationEvents
from pvder.simulation_utilities import SimulationResults

from unittest_utilities import show_DER_status, plot_DER_trajectories

#working_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
#module_folder = os.path.dirname(os.path.abspath(__file__))

#sys.path.append(module_folder)
#os.environ['PATH']+=';'+module_folder
config_file = r'..\config_der.json'

def suite():
	"""Define a test suite."""
	avoid_tests = ['HVRT2','HVRT3']
	tests = test_options.options 
	tests = list(set(tests) - set(avoid_tests))
	print('Following unittest scenarios will be run:{}'.format(tests))
	suite = unittest.TestSuite()
	
	for test in tests:
		suite.addTest(TestPVDER('test_PVDER_'+ test))
											   
	return suite

class TestPVDER(unittest.TestCase):
	
	Vnominal = 1.0
	default_setup_scenario = 'default'
	default_test_scenario = 'default'
	
	setup_scenarios ={'default':{'powerRating':50e3,'SinglePhase':False,'steadyStateInitialization':True,'n_DER':1},
					 'case1':{'SinglePhase':True},
					 'case2':{'SteadyState':False}}
	
	test_scenarios = {'default':{'HVRT_ENABLE':True,'tEnd':5.0,'Vnominal':1.0,'Vspike':1.1},
					  'HVRT1':{'HVRT_ENABLE':True,'tEnd':10.0,'tspike_start':2.0,'tspike_duration':0.5,'Vspike':1.1},
					  'HVRT2':{'HVRT_ENABLE':True,'tEnd':10.0,'tspike_start':2.0,'tspike_duration':4.0,'Vspike':1.1},
					  'HVRT3':{'HVRT_ENABLE':True,'tEnd':10.0,'tspike_start':2.0,'tspike_duration':4.0,'Vspike':1.1}}
		
	
	
	def test_PVDER_default(self):
		"""Test PV-DER without HVRT."""
			
		self.setup_simulation(scenario='default')
		
		n_time_steps = self.run_simulation(scenario='default')		
		
		self.loop_and_check(n_time_steps) 
	
	def test_PVDER_HVRT1(self):
		"""Test PV-DER with HVRT."""
			
		self.setup_simulation(scenario='default')
		
		n_time_steps = self.run_simulation(scenario='HVRT1')		
		
		self.loop_and_check(n_time_steps) 
	
	def test_PVDER_HVRT2(self):
		"""Test PV-DER with HVRT."""			
		
		self.setup_simulation(scenario='default')
		
		n_time_steps = self.run_simulation(scenario='HVRT2')		
		
		self.loop_and_check(n_time_steps)	 
	
	def test_PVDER_HVRT3(self):
		"""Test PV-DER with HVRT."""
			
		self.setup_simulation(scenario='default')
		
		n_time_steps = self.run_simulation(scenario='HVRT3')		
		
		self.loop_and_check(n_time_steps)		
	
	def setup_simulation(self,scenario='default'):
		"""Setup a simulation."""		
		
		self.events_list = []
		self.grid_list = []
		self.DER_model_list =[]
		self.sim_list = []
		self.results_list = []
		
		self.n_instances = self.return_settings(scenario=scenario,parameter='n_DER',settings_type='setup')
		
		flag_arguments = {'standAlone': True,
						  'steadyStateInitialization':self.return_settings(scenario=scenario,parameter='steadyStateInitialization',settings_type='setup'),
						  'verbosity':'DEBUG'}
		ratings_arguments ={'powerRating': self.return_settings(scenario=scenario,parameter='powerRating',settings_type='setup')}
		   
		
		SinglePhase = self.return_settings(scenario=scenario,parameter='SinglePhase',settings_type='setup')
		
		for i in range(self.n_instances):
			
			self.events_list.append(SimulationEvents())
			self.grid_list.append(Grid(events=self.events_list[-1]))

			kwargs={"gridModel":self.grid_list[-1],"identifier":scenario}
			kwargs.update(flag_arguments)
			kwargs.update(ratings_arguments)

			if SinglePhase:
				self.DER_model_list.append(SolarPVDERSinglePhase(events=self.events_list[-1],\
				configFile=config_file,**kwargs))
			else:
				self.DER_model_list.append(SolarPVDERThreePhase(events=self.events_list[-1],\
				configFile=config_file,**kwargs))

			self.sim_list.append(DynamicSimulation(gridModel=self.grid_list[-1],PV_model=self.DER_model_list[-1],events = self.events_list[-1],LOOP_MODE=False,COLLECT_SOLUTION=True))
			self.sim_list[-1].jacFlag = False	  #Provide analytical Jacobian to ODE solver
			self.sim_list[-1].DEBUG_SOLVER = False #Check whether solver is failing to converge at any step
			
			self.results_list.append(SimulationResults(simulation = self.sim_list[-1],PER_UNIT=True))
	
	def run_simulation(self,scenario):
		"""Test PV-DER with HVRT."""
		
		dt=1/120
		
		for i in range(self.n_instances):
			self.specify_scenario(self.events_list[i],self.DER_model_list[i],self.sim_list[i],scenario=scenario)
			n_time_steps = int(self.sim_list[i].tStop/dt) 
			self.sim_list[i].run_simulation()		
				
		return n_time_steps  
	
	def specify_scenario(self,events,DER_model,sim,scenario='HVRT1'):
		"""Specify scenario for unit test."""
		
		Vnominal = self.return_settings(scenario=scenario,parameter='Vnominal',settings_type='test')
		Vspike = self.return_settings(scenario=scenario,parameter='Vspike',settings_type='test')
		
		tspike_start = self.return_settings(scenario=scenario,parameter='tspike_start',settings_type='test')
		tspike_duration = self.return_settings(scenario=scenario,parameter='tspike_duration',settings_type='test')
		
		DER_model.HVRT_ENABLE = self.return_settings(scenario=scenario,parameter='HVRT_ENABLE',settings_type='test')	   
		
		sim.tStop  = self.return_settings(scenario=scenario,parameter='tEnd',settings_type='test')	   
		sim.name = scenario+'-'+sim.name
		if tspike_start is not None:
			events.add_grid_event(tspike_start,Vspike)
		
		if tspike_duration is not None:
			events.add_grid_event(tspike_start+tspike_duration,Vnominal)	
				
	def return_settings(self,scenario,parameter,settings_type):
		"""Check and repalce with correct parameter from scenario."""
		
		if settings_type == 'test':
			if parameter in self.test_scenarios[scenario]:
				_variable = self.test_scenarios[scenario][parameter] 
			else:
				if parameter in self.test_scenarios[self.default_test_scenario]:
					_variable = self.test_scenarios[self.default_test_scenario][parameter]
				else:
					_variable = None
		
		elif settings_type == 'setup':
			if parameter in self.setup_scenarios[scenario]:
				_variable = self.setup_scenarios[scenario][parameter] 
			else:
				if parameter in self.setup_scenarios[self.default_setup_scenario]:
					_variable = self.setup_scenarios[self.default_setup_scenario][parameter]
				else:
					_variable = None			
	   
		return _variable
	
	def loop_and_check(self,n_time_steps):
		"""Loop trhough DER instances and check."""		
		
		for i in range(self.n_instances):
			pvder_object =  self.DER_model_list[i]
			sim_object =  self.sim_list[i]
			results_object =  self.results_list[i]
			
			for convergence_failure in sim_object.convergence_failure_list:
				print('Failure event:{}'.format(convergence_failure))
			
			show_DER_status(pvder_object)
			self.check_DER_state(pvder_object)
			plot_DER_trajectories(results_object)
			self.check_HVRT_status(pvder_object)
			 
			#self.assertTrue(len(sim_object.t_t) == len(sim_object.Vdc_t) == n_time_steps+1, msg='{}:The number of states collected should be {} but it is actually {}!'.format(sim_object.name,n_time_steps+1,len(sim_object.t_t)))
	
	def check_DER_state(self,pvder_object):
		"""Check whether DER states are nominal."""
		
		#Check if DC link voltage within inverter limits.
		self.assertTrue(pvder_object.Vdc*pvder_object.Vdcbase >= pvder_object.Vdcmpp_min or pvder_object.Vdc*pvder_object.Vdcbase <= pvder_object.Vdcmpp_max, msg='{}:DC link voltage {:.2f} V exceeded limit!'.format(pvder_object.name,pvder_object.Vdc*pvder_object.Vdcbase))
		
		#Check current reference and power output within inverter limits.
		self.assertTrue(abs(pvder_object.ia_ref)<= pvder_object.iref_limit, msg='{}:Inverter current exceeded limit by {:.2f} A!'.format(pvder_object.name,(abs(pvder_object.ia_ref) - pvder_object.iref_limit)*pvder_object.Ibase))
		
		self.assertTrue(abs(pvder_object.S_PCC)<=pvder_object.Sinverter_nominal, msg='{}:Inverter power output exceeded limit by {:.2f}  VA!'.format(pvder_object.name,(abs(pvder_object.S_PCC) -pvder_object.Sinverter_nominal)*pvder_object.Sbase))
	
	def check_HVRT_status(self,pvder_object):
		"""Check whether ride through is working."""
		
				 
		if pvder_object.DER_TRIP: #Check if DER trip flag is True
			#Check if HVRT momentary cessation is True is connected
			self.assertTrue(pvder_object.HVRT_TRIP, msg='{}: HVRT trip should be true!'.format(pvder_object.name))
			#Check if DER is connected
			self.assertFalse(pvder_object.DER_CONNECTED, msg='{}: DER connected despite trip!'.format(pvder_object.name))
			#Check if DER stopped supplying power
			self.assertAlmostEqual(abs(pvder_object.S_PCC), 0.0, places=4, msg='{}:Inverter power output is {:.2f} VA despite trip status!'.format(pvder_object.name,pvder_object.S_PCC*pvder_object.Sbase))
			#Check if DC link voltage limits are breached
			self.assertTrue(pvder_object.Vdc*pvder_object.Vdcbase >= pvder_object.Vdcmpp_min or pvder_object.Vdc*pvder_object.Vdcbase <= pvder_object.Vdcmpp_max, msg='{}:DC link voltage exceeded limits!'.format(pvder_object.name))
				
		elif pvder_object.DER_MOMENTARY_CESSATION:
			#Check if HVRT momentary cessation is True is connected
			self.assertTrue(pvder_object.HVRT_MOMENTARY_CESSATION, msg='{}: HVRT momentary cessation should be true!'.format(pvder_object.name))
			#Check if DER is connected
			self.assertFalse(pvder_object.DER_CONNECTED, msg='{}: DER connected despite momentary cessation!'.format(pvder_object.name))
		
parser = argparse.ArgumentParser(description='Unit tests for HVRT operation in OpenDSS - PVDER simulation.')

test_options = TestPVDER.test_scenarios.keys()
#test_options.sort()
test_options = sorted(test_options)

parser.add_argument('-i', '--item', action='store', dest='options',
					type=str, nargs='*', default=test_options,choices=test_options,
					help="Examples: -i HVRT1 HVRT2, -i HVRT3")
test_options = parser.parse_args()
		
if __name__ == '__main__':
	#unittest.main()
	
	runner = unittest.TextTestRunner()
	runner.run(suite())
