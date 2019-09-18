from __future__ import division
import sys
import os
import argparse
import logging
import unittest

import math

import matplotlib.pyplot as plt

from pvder.DER_components_single_phase import SolarPV_DER_SinglePhase
from pvder.grid_components import Grid
from pvder.dynamic_simulation import DynamicSimulation
from pvder.simulation_events import SimulationEvents
from pvder.simulation_utilities import SimulationResults

from unittest_utilities import show_DER_status, plot_DER_trajectories

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
    
    def test_init(self):
        """Test PV-DER three phase mode."""          
                        
        events = SimulationEvents()
                
        PVDER = SolarPV_DER_SinglePhase(events = events,
                                       Sinverter_rated = self.power_rating,Vrms_rated = self.Vrms, #175
                                       gridVoltagePhaseA = self.Va,
                                       gridVoltagePhaseB = self.Vb,
                                       gridVoltagePhaseC = self.Vc,
                                       gridFrequency = self.wgrid,
                                       standAlone = False,STEADY_STATE_INITIALIZATION=True,verbosity = 'DEBUG')
    
        self.assertIsInstance(PVDER, SolarPV_DER_SinglePhase)
        self.assertTrue(PVDER.STEADY_STATE_INITIALIZATION)      
        
    def test_parameter_dict(self):
        """Test initalization and update of paraemter dictionary.""" 
        
        source_ID = '10'
        new_ID ='test_DER'
        new_module_parameters = {'Np':5}
        new_circuit_parameters = {'C_actual':100.0e-6}
                
        events = SimulationEvents()
                
        PVDER = SolarPV_DER_SinglePhase(events = events,
                                       Sinverter_rated = self.power_rating,Vrms_rated = self.Vrms, #175
                                       gridVoltagePhaseA = self.Va,
                                       gridVoltagePhaseB = self.Vb,
                                       gridVoltagePhaseC = self.Vc,
                                       gridFrequency = self.wgrid,
                                       standAlone = False,STEADY_STATE_INITIALIZATION=True,verbosity = 'INFO')   
    
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
                
        PVDER = SolarPV_DER_SinglePhase(events = events,
                                       Sinverter_rated = self.power_rating,Vrms_rated = self.Vrms, #175
                                       gridVoltagePhaseA = self.Va,
                                       gridVoltagePhaseB = self.Vb,
                                       gridVoltagePhaseC = self.Vc,
                                       gridFrequency = self.wgrid,
                                       standAlone = False,STEADY_STATE_INITIALIZATION=True,verbosity = 'DEBUG')
    
        
        jac_CHECK,Jn,Ja = PVDER.check_jacobian()
        self.assertTrue(jac_CHECK,'Analytical and numerical Jacobian should be same.')      
        self.assertEqual(Jn.shape,(PVDER.n_ODE,PVDER.n_ODE))    
    
if __name__ == '__main__':
    #unittest.main()
    logging.debug('test')
    runner = unittest.TextTestRunner()
    runner.run(suite())