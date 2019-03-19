import unittest
import pdb
import numpy as np
import sys, os

from openDSSAPI import Worker

from pvder.DER_components import SolarPV_DER
from pvder.grid_components import Grid
from pvder.simulation_events import SimulationEvents
from pvder.dynamic_simulation import GridSimulation
import pvder.utility_functions as utility_functions


class TestPVDER(unittest.TestCase):
    
    def test_PVDER(self):
        """Test individual PV-DER model with change in solar insolation events."""
        events = SimulationEvents()
        grid = Grid(events=events1,unbalance_ratio_b=1.0,unbalance_ratio_c=1.0)
        PV_DER = SolarPV_DER(grid_model=grid,events=events,Sinverter_rated = 50.0e3,standAlone = True,STEADY_STATE_INITIALIZATION=False)
        sim = GridSimulation(grid_model=grid,PV_model=PV_DER,simulation_events = events1)
        results = SimulationResults(simulation = sim,PER_UNIT=True)
        
        #Loop through all PV-DER instances and initialize
        events.add_solar_event(1.0,50,300)  #Add new solar even at 1.0 s
        events.add_solar_event(1.5,100,300)  #Add new solar even at 1.0 s
        PV_DER.MPPT_ENABLE=False
        PV_DER.RAMP_ENABLE = False
        PV_DER.VOLT_VAR_ENABLE = False
        PV_DER.LVRT_ENABLE = False
        PV_DER.LFRT_ENABLE = False
        PV_DER.DO_EXTRA_CALCULATIONS = True
        sim.jacFlag = True
        sim.DEBUG_SIMULATION = False
        sim.DEBUG_VOLTAGES = True
        sim.DEBUG_CURRENTS = True
        sim.DEBUG_POWER = False
        sim.DEBUG_CONTROLLERS  = True
        sim.DEBUG_PLL = False
        sim.PER_UNIT = True
        sim.DEBUG_SOLVER  = True
        sim.tStop = 20.0
        sim.tInc = 0.001
        sim.run_simulation()
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
        
        #Check if duty cycle less that 1.0.
        self.assertTrue(abs(PV_DER.ma)<=1.0 and abs(PV_DER.mb)<=1.0 and abs(PV_DER.mc)<=1.0, msg='Duty cycle {},{},{} is  greater than 1!'.format(PV_DER.ma,PV_DER.mb,PV_DER.mc))
        #Check if DC link voltage within inverter limits.
        self.assertTrue(PV_DER.Vdc*PV_DER.Vdcbase >= PV_DER.Vdcmpp_min or PV_DER.Vdc*PV_DER.Vdcbase <= PV_DER.Vdcmpp_max, msg='DC link voltage {} exceeded limit!'.format(PV_DER.Vdc))
        #Check current reference and power output within inverter limits.
        self.assertTrue(abs(PV_DER.ia_ref)<= PV_DER.iref_limit, msg='Inverter current exceeded limit!')
        self.assertTrue(abs(PV_DER.S_PCC)<= PV_DER.Sinverter_nominal, msg='Inverter power output exceeded limit!')
        
if __name__ == '__main__':
    unittest.main()