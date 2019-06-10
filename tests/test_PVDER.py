from __future__ import division
import sys
import os
import math
import argparse

import matplotlib.pyplot as plt
import unittest

from pvder.DER_components_single_phase import SolarPV_DER_SinglePhase
from pvder.DER_components_three_phase import SolarPV_DER_ThreePhase
from pvder.grid_components import Grid
from pvder.dynamic_simulation import DynamicSimulation
from pvder.simulation_events import SimulationEvents
from pvder.simulation_utilities import SimulationResults

#working_folder = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
#module_folder = os.path.dirname(os.path.abspath(__file__))

#sys.path.append(module_folder)
#os.environ['PATH']+=';'+module_folder


def suite():
    """Define a test suite."""
    
    tests = test_options.options #['LVRT1'] #,,'LVRT2','LVRT3'
    print('Following unittest scenarios will be run:{}'.format(tests))
    suite = unittest.TestSuite()
    
    for test in tests:
        suite.addTest(TestPVDER('test_PVDER_'+ test))
                                               
    return suite

class TestPVDER(unittest.TestCase):
    
    Vnominal = 1.0
    default_setup_scenario = 'default'
    default_test_scenario = 'default'
    
    setup_scenarios ={'default':{'power_rating':50e3,'SinglePhase':False,'SteadyState':True,'n_DER':1},
                     'case1':{'SinglePhase':True},
                     'case2':{'SteadyState':False}}
    
    test_scenarios = {'default':{'LVRT_ENABLE':True,'LVRT_INSTANTANEOUS_TRIP':False,'LVRT_MOMENTARY_CESSATION':False,'tEnd':5.0,'Vnominal':1.0,'Vfault':0.8},
                      'LVRT1':{'LVRT_ENABLE':True,'LVRT_INSTANTANEOUS_TRIP':False,'LVRT_MOMENTARY_CESSATION':False,'tEnd':10.0,'tfault_start':4.0,'tfault_duration':0.5,'Vfault':0.7},
                      'LVRT2':{'LVRT_ENABLE':True,'LVRT_INSTANTANEOUS_TRIP':False,'LVRT_MOMENTARY_CESSATION':True,'tEnd':10.0,'tfault_start':4.0,'tfault_duration':2.0,'Vfault':0.7},
                      'LVRT3':{'LVRT_ENABLE':True,'LVRT_INSTANTANEOUS_TRIP':True,'LVRT_MOMENTARY_CESSATION':False,'tEnd':10.0,'tfault_start':4.0,'tfault_duration':2.0,'Vfault':0.7}}
        
    
    def test_PVDER_default(self):
        """Test PV-DER without LVRT."""
            
        self.setup_simulation(scenario='default')
        
        n_time_steps = self.run_simulation(scenario='default')        
        
        self.loop_and_check(n_time_steps) 
    
    def test_PVDER_LVRT1(self):
        """Test PV-DER with LVRT."""
            
        self.setup_simulation(scenario='default')
        
        n_time_steps = self.run_simulation(scenario='LVRT1')        
        
        self.loop_and_check(n_time_steps) 
    
    def test_PVDER_LVRT2(self):
        """Test PV-DER with LVRT."""            
        
        self.setup_simulation(scenario='default')
        
        n_time_steps = self.run_simulation(scenario='LVRT2')        
        
        self.loop_and_check(n_time_steps)     
    
    def test_PVDER_LVRT3(self):
        """Test PV-DER with LVRT."""
            
        self.setup_simulation(scenario='default')
        
        n_time_steps = self.run_simulation(scenario='LVRT3')        
        
        self.loop_and_check(n_time_steps)        
    
    def setup_simulation(self,scenario='default'):
        """Setup a simulation."""        
        
        self.events_list = []
        self.grid_list = []
        self.DER_model_list =[]
        self.sim_list = []
        self.results_list = []
        
        self.n_instances = self.return_settings(scenario=scenario,parameter='n_DER',settings_type='setup')
        
        power_rating = self.return_settings(scenario=scenario,parameter='power_rating',settings_type='setup')
        SinglePhase = self.return_settings(scenario=scenario,parameter='SinglePhase',settings_type='setup')
        SteadyState = self.return_settings(scenario=scenario,parameter='SteadyState',settings_type='setup')
        
        for i in range(self.n_instances):
            
            self.events_list.append(SimulationEvents())
            self.grid_list.append(Grid(events=self.events_list[-1]))
        
            if SinglePhase:
                self.DER_model_list.append(SolarPV_DER_SinglePhase(grid_model=self.grid_list[-1],events=self.events_list[-1],
                                                                   standAlone=True,
                                                                   Sinverter_rated = power_rating,
                                                                   STEADY_STATE_INITIALIZATION=SteadyState))
            else:
                self.DER_model_list.append(SolarPV_DER_ThreePhase(grid_model=self.grid_list[-1],events=self.events_list[-1],
                                                                  standAlone=True,
                                                                  Sinverter_rated = power_rating,
                                                                  STEADY_STATE_INITIALIZATION=SteadyState))

            self.sim_list.append(DynamicSimulation(grid_model=self.grid_list[-1],PV_model=self.DER_model_list[-1],events = self.events_list[-1],LOOP_MODE=False,COLLECT_SOLUTION=True))
            self.sim_list[-1].jacFlag = True      #Provide analytical Jacobian to ODE solver
            self.sim_list[-1].DEBUG_SOLVER = False #Check whether solver is failing to converge at any step
            
            self.results_list.append(SimulationResults(simulation = self.sim_list[-1],PER_UNIT=True))
    
    def run_simulation(self,scenario):
        """Test PV-DER with LVRT."""
        
        dt=1/120
        
        for i in range(self.n_instances):
            self.specify_scenario(self.events_list[i],self.DER_model_list[i],self.sim_list[i],scenario=scenario)
            n_time_steps = int(self.sim_list[i].tStop/dt) 
            self.sim_list[i].run_simulation()        
                
        return n_time_steps  
    
    def specify_scenario(self,events,DER_model,sim,scenario='LVRT1'):
        """Specify scenario for unit test."""
        
        _Vnominal = self.return_settings(scenario=scenario,parameter='Vnominal',settings_type='test')
        _Vfault = self.return_settings(scenario=scenario,parameter='Vfault',settings_type='test')
        
        _tfault_start = self.return_settings(scenario=scenario,parameter='tfault_start',settings_type='test')
        _tfault_duration = self.return_settings(scenario=scenario,parameter='tfault_duration',settings_type='test')
        
        DER_model.LVRT_ENABLE = self.return_settings(scenario=scenario,parameter='LVRT_ENABLE',settings_type='test')
        DER_model.LVRT_INSTANTANEOUS_TRIP = self.return_settings(scenario=scenario,parameter='LVRT_INSTANTANEOUS_TRIP',settings_type='test')
        DER_model.LVRT_MOMENTARY_CESSATION = self.return_settings(scenario=scenario,parameter='LVRT_MOMENTARY_CESSATION',settings_type='test')
        
        sim.tStop  = self.return_settings(scenario=scenario,parameter='tEnd',settings_type='test')       
        
        if _tfault_start is not None:
            events.add_grid_event(_tfault_start,_Vfault)
        
        if _tfault_duration is not None:
            events.add_grid_event(_tfault_start+_tfault_duration,_Vnominal)    
                
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
            self.show_DER_status(pvder_object)
            self.check_DER_state(pvder_object)
            self.plot_DER_trajectories(results_object)
            self.check_LVRT_status(pvder_object)
             
            #self.assertTrue(len(sim_object.t_t) == len(sim_object.Vdc_t) == n_time_steps+1, msg='{}:The number of states collected should be {} but it is actually {}!'.format(sim_object.name,n_time_steps+1,len(sim_object.t_t)))
    
    def check_DER_state(self,pvder_object):
        """Check whether DER states are nominal."""
        
        #Check if DC link voltage within inverter limits.
        self.assertTrue(pvder_object.Vdc*pvder_object.Vdcbase >= pvder_object.Vdcmpp_min or pvder_object.Vdc*pvder.PV_model.Vdcbase <= pvder_object.Vdcmpp_max, msg='{}:DC link voltage {:.2f} V exceeded limit!'.format(pvder_object.name,pvder_object.Vdc*pvder_object.Vdcbase))
        
        #Check current reference and power output within inverter limits.
        self.assertTrue(abs(pvder_object.ia_ref)<= pvder_object.iref_limit, msg='{}:Inverter current exceeded limit by {:.2f} A!'.format(pvder_object.name,(abs(pvder_object.ia_ref) - pvder_object.iref_limit)*pvder_object.Ibase))
        
        self.assertTrue(abs(pvder_object.S_PCC)<=pvder_object.Sinverter_nominal, msg='{}:Inverter power output exceeded limit by {:.2f}  VA!'.format(pvder_object.name,(abs(pvder_object.S_PCC) -pvder_object.Sinverter_nominal)*pvder_object.Sbase))
    
    def check_LVRT_status(self,pvder_object):
        """Check whether ride through is working."""
                 
        if pvder_object.Vrms <= pvder_object.V_LV2:
                #Check if LVRT trip flag is True
                self.assertTrue(pvder_object.LVRT_TRIP, msg='{}: Inverter trip flag  not set despite low voltage!'.format(pvder_object.name))
                #Check if Inverter stopped supplying power
                self.assertAlmostEqual(abs(pvder_object.S_PCC), 0.0, places=4, msg='{}:Inverter power output is {:.2f} VA despite trip status!'.format(pvder_object.name,pvder_object.S_PCC*pvder_object.Sbase))
                #Check if DC link voltage limits are breached
                self.assertTrue(pvder_object.Vdc*pvder_object.Vdcbase >= pvder_object.Vdcmpp_min or pvder_object.Vdc*pvder.PV_model.Vdcbase <= pvder_object.Vdcmpp_max, msg='{}:DC link voltage exceeded limits!'.format(pvder_object.name))
                
        elif pvder_object.Vrms > pvder_object.V_LV2 and pvder_object.LVRT_MOMENTARY_CESSATION:
                #Check if LVRT trip flag is False
                self.assertFalse(pvder_object.LVRT_TRIP, msg='{}: Inverter trip flag set despite nominal voltage!'.format(pvder_object.name))
        
    def show_DER_status(self,pvder_object):
        """Show DER states."""     
        
        pvder_object.show_PV_DER_states(quantity='power')
        pvder_object.show_PV_DER_states(quantity='current')
        pvder_object.show_PV_DER_states(quantity='voltage')                
        pvder_object.show_PV_DER_states(quantity='duty cycle')
        pvder_object.show_RT_settings(settings_type='LVRT')    
    
    def plot_DER_trajectories(self,results_object):
        """PLot DER trajectory."""
        
        results_object.PER_UNIT = False
        results_object.PLOT_TITLE = True
        results_object.font_size = 18
        
        results_object.plot_DER_simulation(plot_type='active_power_Ppv_Pac_PCC')#
        results_object.plot_DER_simulation(plot_type='reactive_power')
        results_object.plot_DER_simulation(plot_type='current')
        results_object.plot_DER_simulation(plot_type='voltage_Vdc')
        results_object.plot_DER_simulation(plot_type='voltage_LV')
        results_object.plot_DER_simulation(plot_type='duty_cycle')              
   
        
parser = argparse.ArgumentParser(description='Unit tests for LVRT operation in OpenDSS - PVDER simulation.')

test_options = TestPVDER.test_scenarios.keys()
#test_options.sort()
test_options = sorted(test_options)
print(test_options)
parser.add_argument('-i', '--item', action='store', dest='options',
                    type=str, nargs='*', default=test_options,choices=test_options,
                    help="Examples: -i LVRT1 LVRT2, -i LVRT3")
test_options = parser.parse_args()
        
if __name__ == '__main__':
    #unittest.main()
    
    runner = unittest.TextTestRunner()
    runner.run(suite())
