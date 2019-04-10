import unittest
import pdb
import numpy as np
import sys, os

from pvder.DER_components_single_phase import SolarPV_DER_SinglePhase
from pvder.DER_components_three_phase import SolarPV_DER_ThreePhase
from pvder.grid_components import Grid
from pvder.simulation_events import SimulationEvents
from pvder.dynamic_simulation import GridSimulation
from pvder.utilities import SimulationUtilities,SimulationResults
from pvder import utility_functions

class TestPVDER(unittest.TestCase):
    
    def test_PVDER(self):
        """Test individual PV-DER model with change in solar insolation events."""
        events = SimulationEvents()
        grid = Grid(events=events,unbalance_ratio_b=1.0,unbalance_ratio_c=1.0)
        
        PV_DER_single_phase = SolarPV_DER_SinglePhase(grid_model=grid,events=events,Sinverter_rated = 10.0e3,standAlone = True,STEADY_STATE_INITIALIZATION=False)
        sim_single_phase = GridSimulation(grid_model=grid,PV_model=PV_DER_single_phase,events = events)
        results_single_phase = SimulationResults(simulation = sim_single_phase,PER_UNIT=True)
        
        PV_DER_three_phase = SolarPV_DER_ThreePhase(grid_model=grid,events=events,Sinverter_rated = 50.0e3,standAlone = True,STEADY_STATE_INITIALIZATION=False)
        sim_three_phase = GridSimulation(grid_model=grid,PV_model=PV_DER_three_phase,events = events)
        results_three_phase = SimulationResults(simulation = sim_three_phase,PER_UNIT=True)
        
        PV_DER_list = [PV_DER_single_phase,PV_DER_three_phase]
        sim_list = [sim_single_phase,sim_three_phase]
        results_list = [results_single_phase,results_three_phase]
        
        #Loop through all PV-DER instances and initialize
        events.add_solar_event(1.0,50,300)  #Add new solar even at 1.0 s
        events.add_solar_event(1.5,100,300)  #Add new solar even at 1.0 s        
        
        for PV_DER in PV_DER_list:
            
            PV_DER.MPPT_ENABLE=False
            PV_DER.RAMP_ENABLE = False
            PV_DER.VOLT_VAR_ENABLE = False
            PV_DER.LVRT_ENABLE = False
            PV_DER.LFRT_ENABLE = False
            PV_DER.DO_EXTRA_CALCULATIONS = True
        
        for sim in sim_list:
            sim.jacFlag = True
            sim.DEBUG_SIMULATION = False
            sim.DEBUG_VOLTAGES = True
            sim.DEBUG_CURRENTS = True
            sim.DEBUG_POWER = False
            sim.DEBUG_CONTROLLERS  = True
            sim.DEBUG_PLL = False
            sim.PER_UNIT = True
            sim.DEBUG_SOLVER  = True
            sim.tStop = 2.0
            sim.tInc = 0.001
            sim.run_simulation()
        
        for PV_DER in PV_DER_list:
            #Print values of current and voltage
            PV_DER.show_PV_DER_states(quantity='voltage')
            PV_DER.show_PV_DER_states(quantity='current')
            PV_DER.show_PV_DER_states(quantity='power')
            
            #Test if PV-DER power outputs validation is successful.
            PV_DER.validate_model(PRINT_ERROR = False)
            
            #Test if power output of PV module is always greater than zero.
            self.assertTrue(PV_DER.Ppv_calc(PV_DER.Vdcrated)>=0.0 and PV_DER.Ppv_calc(PV_DER.Vdcrated)<=PV_DER.Sinverter_rated,
                            msg='{}: PV module power output {:.3f} W is not within expected limits ({} to {})!'.format(PV_DER.name,PV_DER.Ppv_calc(PV_DER.Vdcrated),0.0,PV_DER.Sinverter_rated))
            
            #Test if PV module actual MPP is in the correct range.
            self.assertTrue(PV_DER.Vdcmpp>=PV_DER.Vdcmpp_min and PV_DER.Vdcmpp<=PV_DER.Vdcmpp_max,
                            msg='{}: MPP voltage {:.3f} V is not within expected limits ({} to {})!'.format(PV_DER.name,PV_DER.Vdcmpp,PV_DER.Vdcmpp_min,PV_DER.Vdcmpp_max))
        
            self.assertAlmostEqual(PV_DER.Pt_phasor, PV_DER.Pt_RMS, places=3, msg='{}:{:.3f} Watt error in active power output!'.format(PV_DER.name,PV_DER.Pt_phasor- PV_DER.Pt_RMS))
            
            self.assertAlmostEqual(PV_DER.Qt_phasor, PV_DER.Qt_RMS, places=3, msg='{}:{:.3f} VAR error in reactive power output!'.format(PV_DER.name,PV_DER.Qt_phasor-PV_DER.Qt_RMS))
        
            self.assertAlmostEqual(PV_DER.Pf_phasor, PV_DER.Pf_RMS, places=3, msg='{}:{:.3f} Watt error in filter power loss!'.format(PV_DER.name,PV_DER.Pf_phasor- PV_DER.Pf_RMS))
        
            self.assertAlmostEqual(PV_DER.Qf_phasor, PV_DER.Qf_RMS, places=3, msg='{}:{:.3f} VAR error in filter reactive consumption!'.format(PV_DER.name,PV_DER.Qf_phasor-PV_DER.Qf_RMS))
        
            self.assertTrue(abs(PV_DER.ma)<=1.0, msg='Duty cycle {} is  greater than 1!'.format(PV_DER.ma))
            if type(PV_DER).__name__ == 'SolarPV_DER_ThreePhase':
                self.assertTrue(abs(PV_DER.mb)<=1.0 and abs(PV_DER.mc)<=1.0, msg='Duty cycle {},{} is  greater than 1!'.format(PV_DER.mb,PV_DER.mc))
        
            #Check if DC link voltage within inverter limits.
            self.assertTrue(PV_DER.Vdc*PV_DER.Vdcbase >= PV_DER.Vdcmpp_min or PV_DER.Vdc*PV_DER.Vdcbase <= PV_DER.Vdcmpp_max, msg='DC link voltage {} exceeded limit!'.format(PV_DER.Vdc*PV_DER.Vdcbase))
            print(abs(PV_DER.S_PCC),PV_DER.Sinverter_nominal)
            #Check current reference and power output within inverter limits.
            self.assertTrue(abs(PV_DER.ia_ref)<= PV_DER.iref_limit, msg='Inverter current exceeded limit by {:.2f} A!'.format((abs(PV_DER.ia_ref) - PV_DER.iref_limit)*PV_DER.Ibase))
            self.assertTrue(abs(PV_DER.S_PCC)<= PV_DER.Sinverter_nominal, msg='Inverter power output exceeded limit by {:.2f}  VA!'.format((abs(PV_DER.S_PCC) - PV_DER.Sinverter_nominal)*PV_DER.Sbase))
        
if __name__ == '__main__':
    unittest.main()