"""Code for setting up, and running, and collecting data from PV-DER simulations."""

from __future__ import division
import numpy as np
import math
import cmath
import time

import pdb
import six
import logging

####from graphviz import Digraph
from pvder.utility_classes import Logging
from pvder.grid_components import Grid
from pvder.simulation_utilities import SimulationUtilities
from pvder import utility_functions
from pvder import config

class ModelUtilities():
    """Class for model wide utilities."""
    
    def __init__(self,PV_model,simulation,grid_model=None):
        self.PV_model = PV_model
        
        self.simulation = simulation
        if self.PV_model.standAlone and grid_model is not None:
            self.grid_model = grid_model
        elif self.PV_model.standAlone and grid_model is None:
            raise ValueError('`Grid` instance need to provided in stand alone mode for creating `GridSimulation` instance!')
    
    def draw_model(self,display_value_type='per_unit'):
        """Draw and render the model using graphs."""
        
        assert display_value_type in ['per_unit','actual'],'Display can only be in per unit or actual values'
        
        self.display_value_type =display_value_type
        if self.display_value_type  == 'per_unit': 
            _C = self.PV_model.C
            _Lf = self.PV_model.Lf
            _Rf = self.PV_model.Rf
            _Z1 = self.PV_model.Z1
            _Z2 = self.grid_model.Z2
        else: 
            _C = self.PV_model.C_actual
            _Lf = self.PV_model.Lf_actual
            _Rf = self.PV_model.Rf_actual
            _Z1 = self.PV_model.Z1_actual
            _Z2 = self.grid_model.Z2_actual
        
        
        dot = Digraph(comment='PV_DER and grid_model.')
        dot.node('Base','Vbase={},Sbase={},Zbase={:.4f},Lbase={:.4f},Cbase={:.4f}'.format(self.grid_model.Vbase,self.grid_model.Sbase,self.grid_model.Zbase,self.grid_model.Lbase,self.grid_model.Cbase),shape='rectangle')
        dot.node('Value_type','{}'.format(self.display_value_type))
        dot.node('DER','{}\nC={:.5f},Lf={:.5f},Rf= {:.4f}'.format(self.PV_model.name,_C,_Lf,_Rf),shape='box')
        dot.node('Transformer','{}\nZl = {:.5f}'.format(self.PV_model.transformer_name,_Z1),shape='rectangle')
        dot.node('Transmission_Line','{}\nZt = {:.5f}'.format(self.grid_model.transmission_name,_Z2),shape='rectangle')
        dot.node('Grid',self.grid_model.name)
       
        dot.edge('DER', 'Transformer')
        dot.edge('Transformer', 'Transmission_Line')
        dot.edge('Transmission_Line', 'Grid')
        dot.render('model_graphs/model.gv', view=True)

class DynamicSimulation(Grid,SimulationUtilities,Logging):
    """ Utility class for running simulations."""
    
    count = 0
    tStart = 0.0
    tInc = config.DEFAULT_DELTA_T
    #jacFlag = False
    DEBUG_SOLVER = False
    DEBUG_SIMULATION = False
    DEBUG_CONTROLLERS = False
    DEBUG_VOLTAGES = False
    DEBUG_CURRENTS = False
    DEBUG_POWER = False
    DEBUG_PLL = False
        
    def __init__(self,PV_model,events,grid_model = None,tStop = 0.5,
                 LOOP_MODE = False,COLLECT_SOLUTION = True,jacFlag = False,
                 verbosity ='INFO',solver_type ='odeint',identifier = None):
        """Creates an instance of `GridSimulation`.
        
        Args:
          PV_model: An instance of `SolarPV_DER`.
          events: An instance of `SimulationEvents`.
          grid_model: An instance of `GridModel` (only need to be suppled in stand alone simulation).
          tStop: A scalar specifying the end time for simulation.
          tInc: A scalar specifying the time step for simulation.
          LOOP_MODE: A boolean specifying whether simulation is run in loop.
        """
        #if LOOP_MODE:
        #    assert not PV_model.standAlone, 'Loop mode can only be true if PV-DER model is stand alone.'
        
        #Increment count to keep track of number of simulation instances
        DynamicSimulation.count = DynamicSimulation.count + 1
        self.name_instance(identifier) #Generate a name for the instance
        
        self.initialize_logger(logging_level=verbosity) #Set logging level - {DEBUG,INFO,WARNING,ERROR}  
       
        self.tStop = tStop       
        self.t = self.t_calc()
        self.PV_model = PV_model
        self.simulation_events = events
        self.simulation_events.del_t_event = self.tInc
        
        self.initialize_solver(solver_type=solver_type)
        
        self.SOLVER_CONVERGENCE = False
        self.convergence_failure_list =[]
        
        self.LOOP_MODE = LOOP_MODE
        self.COLLECT_SOLUTION = COLLECT_SOLUTION
        self.jacFlag = jacFlag
        
        if self.PV_model.standAlone and grid_model is not None:
            self.grid_model = grid_model
        elif self.PV_model.standAlone and grid_model is None:
            raise ValueError('`Grid` instance need to provided in stand alone mode for creating `GridSimulation` instance!')
        
        #Remove existing simulation events
        #self.simulation_events.remove_solar_event(3.0)
        #self.simulation_events.remove_load_event(4.0)
        #self.simulation_events.remove_grid_event(5.0)
        
        self.solution_time = None #Always reset solution time to None
        
        if self.LOOP_MODE:
            self.reset_stored_trajectories()
            
        self.initialize_y0_t()
    
    #@property
    #def y0(self):
    #    """ Combine all initial conditions."""

    #    y0 = self.PV_model.y0
        
    #    return y0
   
    @property
    def y0(self):
        """ Combine all initial conditions from solution."""
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            y0 = [self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],\
                    self.ibR_t[-1], self.ibI_t[-1], self.xbR_t[-1], self.xbI_t[-1], self.ubR_t[-1],self.ubI_t[-1],\
                    self.icR_t[-1], self.icI_t[-1], self.xcR_t[-1], self.xcI_t[-1], self.ucR_t[-1],self.ucI_t[-1],\
                    self.Vdc_t[-1],self.xDC_t[-1],self.xQ_t[-1],self.xPLL_t[-1],self.wte_t[-1]]    
        
        elif type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
            y0 =[self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],\
                   self.Vdc_t[-1],self.xDC_t[-1],self.xQ_t[-1],self.xPLL_t[-1],self.wte_t[-1]]
        
        return y0
    
    #@property
    def t_calc(self):
        """Vector of time steps for simulation"""
        
        #if (self.tStop - self.tStart) <= self.tInc:
        #    self.tStop =  self.tStart  + 1e-6 #+ self.tInc
        return np.arange(self.tStart, self.tStop + self.tInc, self.tInc)
    
    def initialize_y0_t(self):
        """Initialize y0_t."""
        
        self.iaR_t = np.array([self.PV_model.y0[0]])
        self.iaI_t = np.array([self.PV_model.y0[1]])
        self.xaR_t = np.array([self.PV_model.y0[2]])
        self.xaI_t = np.array([self.PV_model.y0[3]])
        self.uaR_t = np.array([self.PV_model.y0[4]])
        self.uaI_t = np.array([self.PV_model.y0[5]])
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
            #DC link voltage variables
            self.Vdc_t = np.array([self.PV_model.y0[6]])
            self.xDC_t = np.array([self.PV_model.y0[7]])
            self.xQ_t = np.array([self.PV_model.y0[8]])
            #PLL variables
            self.xPLL_t = np.array([self.PV_model.y0[9]])
            #Frequency integration to get angle
            self.wte_t = np.array([self.PV_model.y0[10]])
        
        elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.ibR_t = np.array([self.PV_model.y0[6]])
            self.ibI_t = np.array([self.PV_model.y0[7]])
            self.xbR_t = np.array([self.PV_model.y0[8]])
            self.xbI_t = np.array([self.PV_model.y0[9]])
            self.ubR_t = np.array([self.PV_model.y0[10]])
            self.ubI_t = np.array([self.PV_model.y0[11]])

            self.icR_t = np.array([self.PV_model.y0[12]])
            self.icI_t = np.array([self.PV_model.y0[13]])
            self.xcR_t = np.array([self.PV_model.y0[14]])
            self.xcI_t = np.array([self.PV_model.y0[15]])
            self.ucR_t = np.array([self.PV_model.y0[16]])
            self.ucI_t = np.array([self.PV_model.y0[17]])

            self.Vdc_t = np.array([self.PV_model.y0[18]])
            self.xDC_t = np.array([self.PV_model.y0[19]])
            self.xQ_t = np.array([self.PV_model.y0[20]])
            self.xPLL_t = np.array([self.PV_model.y0[21]])
            self.wte_t = np.array([self.PV_model.y0[22]])

    def reset_stored_trajectories(self):
        """Reset for plotting."""
        
        self._t_t = np.array(0.0)
        self.Vdc_t = self._Vdc_t = np.array(self.PV_model.Vdc)
        
        self.ia_t = self._ia_t = np.array(self.PV_model.ia)
        self.ma_t = self._ma_t = np.array(self.PV_model.ma)
        self.vta_t = self._vta_t = np.array(self.PV_model.vta)
        self.va_t = self._va_t = np.array(self.PV_model.va)
        
        self.phita_t = self._phita_t = np.array(np.angle(self.PV_model.vta))
        self.phia_t = self._phia_t = np.array(np.angle(self.PV_model.va))
        
        self.wgrid_t = self._wgrid_t = np.array(self.PV_model.wgrid_measured)
        self.we_t = self._we_t = np.array(self.PV_model.we)
        self.wte_t = self._wte_t = np.array(self.PV_model.wte)
        self.vd_t = self._vd_t = np.array(self.PV_model.vd)
        self.vq_t = self._vq_t = np.array(self.PV_model.vq)
                
        self.ma_absolute_t = self._ma_absolute_t = np.array(abs(self.PV_model.ma))
        self.Varms_t  = self._Varms_t = np.array(abs(self.PV_model.va)/math.sqrt(2))
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.mb_absolute_t = self._mb_absolute_t = np.array(abs(self.PV_model.mb))
            self.mc_absolute_t = self._mc_absolute_t = np.array(abs(self.PV_model.mc))
            
            self.Vbrms_t  = self._Vbrms_t = np.array(abs(self.PV_model.vb)/math.sqrt(2))
            self.Vcrms_t  = self._Vcrms_t = np.array(abs(self.PV_model.vc)/math.sqrt(2))
            
            self.ib_t = self._ib_t = np.array(self.PV_model.ib)
            self.mb_t = self._mb_t = np.array(self.PV_model.mb)
            self.vtb_t = self._vtb_t = np.array(self.PV_model.vtb)
            self.vb_t = self._vb_t = np.array(self.PV_model.vb)
            
            self.phitb_t = self._phitb_t = np.array(np.angle(self.PV_model.vtb))
            self.phib_t = self._phib_t = np.array(np.angle(self.PV_model.vb))
            
            self.ic_t = self._ic_t = np.array(self.PV_model.ic)
            self.mc_t = self._mc_t = np.array(self.PV_model.mc)
            self.vtc_t = self._vtc_t = np.array(self.PV_model.vtc)
            self.vc_t = self._vc_t = np.array(self.PV_model.vc)
            
            self.phitc_t = self._phitc_t = np.array(np.angle(self.PV_model.vtc))
            self.phic_t = self._phic_t = np.array(np.angle(self.PV_model.vc))
            
        self.Irms_t = self._Irms_t = np.array(self.PV_model.Irms)
        self.Ppv_t = self._Ppv_t = np.array(self.PV_model.Ppv)
        self.S_PCC_t = self._S_PCC_t = np.array(self.PV_model.S_PCC)
        self.S_t = self._S_t = np.array(self.PV_model.S)
        self.Vtrms_t = self._Vtrms_t = np.array(self.PV_model.Vtrms)
        self.Vrms_t = self._Vrms_t = np.array(self.PV_model.Vrms)
    
    def ODE_model(self,y,t):
        """ Combine all derivatives."""
        
        y1 = y[0:self.PV_model.n_ODE]
        
        if self.PV_model.standAlone:
            self.grid_model.steady_state_model(t)
            
        y = self.PV_model.ODE_model(y1,t)
        
        if self.DEBUG_SIMULATION:
            self.debug_simulation(t)

        return y
        
    def jac_ODE_model(self,y,t):
        """ Combine all derivatives."""
        
        y1 = y[0:self.PV_model.n_ODE]
        
        if self.PV_model.standAlone:
            self.grid_model.steady_state_model(t)
            
        y = self.PV_model.jac_ODE_model(y1,t)

        return y
        
    def debug_simulation(self,t):
        """ Print to terminal for debugging."""
      
        utility_functions.print_to_terminal('t:{:.4f}'.format(t))
       
        if self.DEBUG_VOLTAGES:
            utility_functions.print_to_terminal('Vdc_ref:{:.3f},Vdc:{:.3f},Vat:{:.3f},Va:{:.3f},Vag:{:.3f},Vagrid:{:.3f},Vagrid_setpoint:{:.3f}'. format(self.PV_model.Vdc_ref,self.PV_model.Vdc,self.PV_model.Vtrms,self.PV_model.Vrms,self.grid_model.Vgrms,abs(self.grid_model.Vagrid_no_conversion)/math.sqrt(2),abs(self.grid_model.Vagrid)/math.sqrt(2)))
            
        if self.DEBUG_CURRENTS:
            utility_functions.print_to_terminal('ia_ref:{:.3f},ia:{:.3f},iload1:{:.3f}'. format(self.PV_model.ia_ref,self.PV_model.ia,self.PV_model.iaload1))
                
        if self.DEBUG_POWER:
            utility_functions.print_to_terminal('Sinsol:{:.3f},Q_ref:{:.3f},Ppv:{:.3f},S:{:.3f},S_PCC:{:.3f},S_load1:{:.3f},S_G:{:.3f}'.format(self.PV_model.Sinsol,self.PV_model.Q_ref,self.PV_model.Ppv,self.PV_model.S,self.PV_model.S_PCC,self.PV_model.S_load1,self.PV_model.S_G))
        
        if self.DEBUG_CONTROLLERS:
            utility_functions.print_to_terminal('xDC:{:.3f},xQ:{:.3f},ua:{:.3f},xa:{:.3f},ma:{:.3f}'. format(self.PV_model.xdc,self.PV_model.xQ,self.PV_model.ua,self.PV_model.xa,self.PV_model.ma))
            
        if self.DEBUG_PLL:
            utility_functions.print_to_terminal("we:{:.3f}, wte:{:.3f} rad, vd: {:.3f} V, vq {:.3f} V".format(self.PV_model.we,self.PV_model.wte,self.PV_model.vd,self.PV_model.vq))
    
    def time_series_PCC_HV_side_voltage(self):
        """Calculate time series PCC voltage."""
        
        assert len(self.ia_t) == len(self.vag_t) != None, "States must be available from simulation."
        
        self.vaHV_t = self.vag_t + (self.ia_t/self.PV_model.a)*self.grid_model.Z2 - (self.va_t/self.PV_model.a)*(self.grid_model.Z2/self.Zload1_t)
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
            self.vbHV_t = self.vbg_t 
            self.vcHV_t = self.vcg_t 
            
        elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.vbHV_t = self.vbg_t + (self.ib_t/self.PV_model.a)*self.grid_model.Z2 - (self.vb_t/self.PV_model.a)*(self.grid_model.Z2/self.Zload1_t)
            self.vcHV_t = self.vcg_t + (self.ic_t/self.PV_model.a)*self.grid_model.Z2 - (self.vc_t/self.PV_model.a)*(self.grid_model.Z2/self.Zload1_t)

    def time_series_duty_cycle(self):
        """Calculate time series PCC voltage."""
        
        self.ma_t =  utility_functions.m_time_series(self.ua_t,self.xa_t,self.PV_model.Kp_GCC)
        self.maR_t = self.ma_t.real 
        self.maI_t = self.ma_t.imag
        
        self.ma_absolute_t = utility_functions.Uabsolute_time_series(self.ma_t)
                
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.mb_t =  utility_functions.m_time_series(self.ub_t,self.xb_t,self.PV_model.Kp_GCC)
            self.mc_t =  utility_functions.m_time_series(self.uc_t,self.xc_t,self.PV_model.Kp_GCC)

            self.mbR_t = self.mb_t.real 
            self.mbI_t = self.mb_t.imag
            self.mcR_t = self.mc_t.real 
            self.mcI_t = self.mc_t.imag
            
            self.mb_absolute_t = utility_functions.Uabsolute_time_series(self.mb_t)
            self.mc_absolute_t = utility_functions.Uabsolute_time_series(self.mc_t)

    def time_series_PCC_LV_side_voltage(self):
        """Calculate time series PCC voltage."""
         
        if self.PV_model.standAlone:
            self.va_t = ((self.vag_t+(self.ia_t/self.PV_model.a)*self.grid_model.Z2)/(self.PV_model.a) +self.ia_t*self.PV_model.Z1)*((self.Zload1_t*self.PV_model.a*self.PV_model.a)/((self.PV_model.a*self.PV_model.a*(self.PV_model.Z1+self.Zload1_t))+self.grid_model.Z2))
            
        else:
            self.va_t = np.repeat(self.PV_model.gridVoltagePhaseA,len(self.t))
            
        self.vaR_t = self.va_t.real
        self.vaI_t = self.va_t.imag
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            if self.PV_model.standAlone:
                self.vb_t = ((self.vbg_t+(self.ib_t/self.PV_model.a)*self.grid_model.Z2)/(self.PV_model.a) +self.ib_t*self.PV_model.Z1)*((self.Zload1_t*self.PV_model.a*self.PV_model.a)/((self.PV_model.a*self.PV_model.a*(self.PV_model.Z1+self.Zload1_t))+self.grid_model.Z2))
                self.vc_t = ((self.vcg_t+(self.ic_t/self.PV_model.a)*self.grid_model.Z2)/(self.PV_model.a) +self.ic_t*self.PV_model.Z1)*((self.Zload1_t*self.PV_model.a*self.PV_model.a)/((self.PV_model.a*self.PV_model.a*(self.PV_model.Z1+self.Zload1_t))+self.grid_model.Z2))
              
            else:
                self.vb_t = np.repeat(self.PV_model.gridVoltagePhaseB,len(self.t))
                self.vc_t = np.repeat(self.PV_model.gridVoltagePhaseC,len(self.t))
                
            self.vbR_t = self.vb_t.real
            self.vbI_t = self.vb_t.imag
            self.vcR_t = self.vc_t.real
            self.vcI_t = self.vc_t.imag
        
    def time_series_inv_terminal_voltage(self):
        """Calculate time series inverter terminal voltage."""
        
        self.vta_t =  utility_functions.Vinv_terminal_time_series(self.ma_t,self.Vdc_t)
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.vtb_t =  utility_functions.Vinv_terminal_time_series(self.mb_t,self.Vdc_t)
            self.vtc_t =  utility_functions.Vinv_terminal_time_series(self.mc_t,self.Vdc_t)
    
    def time_series_RMS(self):
        """Calculate time series RMS quantities."""
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
            self.Vtrms_t = utility_functions.Urms_time_series(self.vta_t,self.vta_t,self.vta_t)
            self.Vrms_t = utility_functions.Urms_time_series(self.va_t,self.va_t,self.va_t)
            self.Irms_t = utility_functions.Urms_time_series(self.ia_t,self.ia_t,self.ia_t)
            
            self.Varms_t = self.Vrms_t
    
        elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.Vtrms_t = utility_functions.Urms_time_series(self.vta_t,self.vtb_t,self.vtc_t)
            self.Vrms_t = utility_functions.Urms_time_series(self.va_t,self.vb_t,self.vc_t)
            self.Irms_t = utility_functions.Urms_time_series(self.ia_t,self.ib_t,self.ic_t)       
            
            self.Varms_t = utility_functions.Uphrms_time_series(self.va_t)
            self.Vbrms_t = utility_functions.Uphrms_time_series(self.vb_t)
            self.Vcrms_t = utility_functions.Uphrms_time_series(self.vc_t)
       
        if self.PV_model.standAlone:
            self.Vgrms_t = utility_functions.Urms_time_series(self.vag_t,self.vbg_t,self.vcg_t)
            self.Vhvrms_t= utility_functions.Urms_time_series(self.vaHV_t,self.vbHV_t,self.vcHV_t)
    
    def time_series_grid(self):
        """Time series grid voltage and frequency for standalone model."""
        
        self.wgrid_t = []
        
        if self.PV_model.standAlone:
            self.vag_t = []
            self.vbg_t = []
            self.vcg_t = []
            
            for i,t in enumerate(self.t):   #Loop through grid events and append voltage and frequency at each time step
                Vagrid_new,self.grid_model.wgrid = self.simulation_events.grid_events(t)
                self.grid_model.vag = Vagrid_new*(self.grid_model.Vgridrated/self.Vbase)
                self.grid_model.vbg = utility_functions.Ub_calc(self.grid_model.vag*self.grid_model.unbalance_ratio_b)
                self.grid_model.vcg = utility_functions.Uc_calc(self.grid_model.vag*self.grid_model.unbalance_ratio_c)  

                self.vag_t.append(self.grid_model.vag)
                self.vbg_t.append(self.grid_model.vbg)
                self.vcg_t.append(self.grid_model.vcg)
                
                self.wgrid_t.append(self.grid_model.wgrid)

            self.vag_t = np.asarray(self.vag_t)
            self.vbg_t = np.asarray(self.vbg_t)
            self.vcg_t = np.asarray(self.vcg_t)
            self.wgrid_t = np.asarray(self.wgrid_t)
                
            self.vagR_t = self.vag_t.real
            self.vagI_t = self.vag_t.imag
                
            self.vbgR_t = self.vbg_t.real
            self.vbgI_t = self.vbg_t.imag
                
            self.vcgR_t = self.vcg_t.real
            self.vcgI_t = self.vcg_t.imag
        
        else:
            self.wgrid_t = np.repeat(self.PV_model.wgrid_measured,len(self.t)) #only frequency comes from outside in standalone mode
        
        if not self.LOOP_MODE:
            self.simulation_events.reset_event_counters() #reset event counters
   
    def time_series_PLL(self):
        """Calculate time series PLL and d-q quantities."""
        
        if self.PV_model.standAlone:
            self.we_t = []
            self.vd_t = []
            self.vq_t = []
        
            for i,t in enumerate(self.t):     #Loop through time steps and calculate d-q values at each time step

                self.PV_model.wte = self.wte_t[i]
                self.PV_model.xPLL = self.xPLL_t[i]

                if type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
                    self.PV_model.valpha = utility_functions.phasor_to_time_1phase(self.va_t[i],self.wgrid_t[i],t)
                    self.PV_model.vbeta =utility_functions.phasor_to_time_1phase(self.va_t[i]*pow(math.e,-1j*(math.pi/2)),self.wgrid_t[i],t)
                    self.PV_model.vd,self.PV_model.vq = utility_functions.alpha_beta_to_d_q(self.PV_model.valpha,self.PV_model.vbeta,self.PV_model.wte)

                elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
                    self.PV_model.vat = utility_functions.phasor_to_time_1phase(self.va_t[i],self.wgrid_t[i],t)
                    self.PV_model.vbt = utility_functions.phasor_to_time_1phase(self.vb_t[i],self.wgrid_t[i],t)
                    self.PV_model.vct = utility_functions.phasor_to_time_1phase(self.vc_t[i],self.wgrid_t[i],t)

                    self.PV_model.vd,self.PV_model.vq,self.PV_model.v0 = utility_functions.abc_to_dq0(self.PV_model.vat,self.PV_model.vbt,self.PV_model.vct,self.PV_model.wte)

                self.PV_model.we = self.PV_model.we_calc() #Calculate inverter frequency from PLL equation
                
                self.we_t.append(self.PV_model.we)
                self.vd_t.append(self.PV_model.vd)
                self.vq_t.append(self.PV_model.vq)
                
            self.we_t = np.asarray(self.we_t)
            self.vd_t = np.asarray(self.vd_t)
            self.vq_t = np.asarray(self.vq_t)
        else:
            self.we_t = np.repeat(self.PV_model.we,len(self.t))
            self.vd_t = np.repeat(self.PV_model.vd,len(self.t))
            self.vq_t = np.repeat(self.PV_model.vq,len(self.t))
        
    def time_series_Zload1(self,tOverride=None):
        """Calculate time series load impedance."""
        
        if tOverride:
            tAll=tOverride
        else:
            tAll=self.t

        self.Zload1_t = []
        for i,t in enumerate(tAll):       #Loop through load events
            _Zload1_actual =  self.simulation_events.load_events(t)  
            _Zload1 = _Zload1_actual/self.PV_model.Zbase
            
            self.Zload1_t.append(_Zload1)
        self.Zload1_t = np.asarray(self.Zload1_t)
        if not self.LOOP_MODE:
           self.simulation_events.reset_event_counters() #reset event counters
    
    def time_series_Ppv(self):
        """Calculate time series Solar PV power output."""
        
        self.Ppv_t = []
        self.Sinsol_t = []
        for i,t in enumerate(self.t):       #Loop through solar events and calculate Ppv at each time step
           self.PV_model.Sinsol,self.PV_model.Tactual =  self.simulation_events.solar_events(t)  #Parse through solar events
           self.PV_model.Vdc = self.Vdc_t[i]
           self.PV_model.Iph = self.PV_model.Iph_calc()
           self.Sinsol_t.append(self.PV_model.Sinsol)
           self.Ppv_t.append(self.PV_model.Ppv_calc(self.PV_model.Vdc_actual))
        self.Sinsol_t = np.asarray(self.Sinsol_t)
        self.Ppv_t = np.asarray(self.Ppv_t)
        if not self.LOOP_MODE:
            self.simulation_events.reset_event_counters() #reset event counters
    
    def time_series_S(self):
        """Calculate time series apparent power."""
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
            self.S_t = (1/2)*(self.vta_t*self.ia_t.conjugate())
            self.S_PCC_t = (1/2)*(self.va_t*self.ia_t.conjugate())            
            
            if self.PV_model.standAlone:
                self.S_G_t = (1/2)*(-(self.ia_t - (self.va_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vag_t
                self.S_load1_t = (1/2)*(self.va_t*(-(self.va_t/self.Zload1_t)).conjugate())
        
        elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.S_t = (1/2)*(self.vta_t*self.ia_t.conjugate() + self.vtb_t*self.ib_t.conjugate() + self.vtc_t*self.ic_t.conjugate())
            self.S_PCC_t = (1/2)*(self.va_t*self.ia_t.conjugate()+self.vb_t*self.ib_t.conjugate()+self.vc_t*self.ic_t.conjugate())
             
            if self.PV_model.standAlone:
                self.S_G_t = (1/2)*(-(self.ia_t - (self.va_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vag_t + (-(self.ib_t -(self.vb_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vbg_t + (-(self.ic_t -(self.vc_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vcg_t
                self.S_load1_t = (1/2)*(self.va_t*(-(self.va_t/self.Zload1_t)).conjugate() + self.vb_t*(-(self.vb_t/self.Zload1_t)).conjugate() + self.vc_t*(-(self.vc_t/self.Zload1_t)).conjugate())
    
    def time_series_phase_angle(self):
        """Calculate time series phase angle between grid and inverter voltage."""
        
        assert len(self.va_t) == len(self.vta_t) != None, "States must be available from simulation."
        
        self.phita_t = np.angle(self.vta_t)
        self.phia_t = np.angle(self.va_t)
        
        if self.PV_model.standAlone:
            self.phi_aHV_t = np.angle(self.vaHV_t) 
            self.phi_ag_t = np.angle(self.vag_t)        
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.phitb_t = np.angle(self.vtb_t)
            self.phitc_t = np.angle(self.vtc_t)
            self.phib_t = np.angle(self.vb_t)
            self.phic_t = np.angle(self.vc_t)
            if self.PV_model.standAlone:
                self.phi_bHV_t = np.angle(self.vbHV_t)
                self.phi_cHV_t = np.angle(self.vcHV_t)
                self.phi_bg_t = np.angle(self.vbg_t)
                self.phi_cg_t = np.angle(self.vcg_t)
        
    def time_series_power_transfer(self):
        """Calculate time series power transfer between grid and PCC."""
        
        self.P_transfer = (self.Vhvrms_t*self.Vgrms_t*np.sin(self.phi_aHV_t-self.phi_ag_t))/np.abs(self.grid_model.Z2)
    
    def collect_solution(self,solution,t=None):
        """Collect solution."""
        
        if self.LOOP_MODE:
            self.t = t
            self.collect_full_trajectory(solution)
            self.collect_last_states()
        else:
            self.collect_full_trajectory(solution)
        
        self.logger.debug('{}:Stored solution for {} time points starting at {:.3f} s and ending at {:.3f} s!'.format(self.name,len(self.t),self.t[0],self.t[-1]))
        
    def collect_last_states(self):
        """Collect states at last time step."""
        
        self._t_t = np.append(self._t_t,self.t[1:])
        self._Vdc_t = np.append(self._Vdc_t,self.Vdc_t[1:])
        
        self._ia_t = np.append(self._ia_t,self.ia_t[1:])
        self._ma_t = np.append(self._ma_t,self.ma_t[1:])
        self._vta_t = np.append(self._vta_t,self.vta_t[1:])
        self._va_t = np.append(self._va_t,self.va_t[1:])
        
        self._phita_t = np.append(self._phita_t,self.phita_t[1:])
        self._phia_t = np.append(self._phia_t,self.phia_t[1:])
        
        self._wgrid_t = np.append(self._wgrid_t,self.wgrid_t[1:])
        self._we_t = np.append(self._we_t,self.we_t[1:])
        self._wte_t = np.append(self._wte_t,self.wte_t[1:])
        self._vd_t = np.append(self._vd_t,self.vd_t[1:])
        self._vq_t = np.append(self._vq_t,self.vq_t[1:])        
        
        self._ma_absolute_t = np.append(self._ma_absolute_t,self.ma_absolute_t[1:])
        self._Varms_t = np.append(self._Varms_t,self.Varms_t[1:])
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self._mb_absolute_t = np.append(self._mb_absolute_t,self.mb_absolute_t[1:])
            self._mc_absolute_t = np.append(self._mc_absolute_t,self.mc_absolute_t[1:])
            
            self._Vbrms_t = np.append(self._Vbrms_t,self.Vbrms_t[1:])
            self._Vcrms_t = np.append(self._Vcrms_t,self.Vcrms_t[1:])
            
            self._ib_t = np.append(self._ib_t,self.ib_t[1:])
            self._mb_t = np.append(self._mb_t,self.mb_t[1:])
            self._vtb_t = np.append(self._vtb_t,self.vtb_t[1:])
            self._vb_t = np.append(self._vb_t,self.vb_t[1:])
            
            self._phitb_t = np.append(self._phitb_t,self.phitb_t[1:])
            self._phib_t = np.append(self._phib_t,self.phib_t[1:])
        
            self._ic_t = np.append(self._ic_t,self.ic_t[1:])
            self._mc_t = np.append(self._mc_t,self.mc_t[1:])
            self._vtc_t = np.append(self._vtc_t,self.vtc_t[1:])
            self._vc_t = np.append(self._vc_t,self.vc_t[1:])
            
            self._phitc_t = np.append(self._phitc_t,self.phitc_t[1:])
            self._phic_t = np.append(self._phic_t,self.phic_t[1:])
        
        self._Irms_t = np.append(self._Irms_t,self.Irms_t[1:])        

        self._Vtrms_t = np.append(self._Vtrms_t,self.Vtrms_t[1:])
        self._Vrms_t = np.append(self._Vrms_t,self.Vrms_t[1:])        
        
        self._Ppv_t = np.append(self._Ppv_t,self.Ppv_t[1:])
        self._S_t = np.append(self._S_t,self.S_t[1:])
        self._S_PCC_t = np.append(self._S_PCC_t,self.S_PCC_t[1:])

    def collect_full_trajectory(self,solution):
        """Collect full solution from solver."""
        
        #Collect states
        self.t_t  = self.t #Time stamps
        
        self.collect_states(solution)
                
        #Time series current
        self.ia_t = self.iaR_t + 1j*self.iaI_t
        #Time series u
        self.ua_t = self.uaR_t + 1j*self.uaI_t
        #Time series x
        self.xa_t = self.xaR_t+1j*self.xaI_t
        assert len(self.ua_t) == len(self.xa_t) != None, "States must be available from simulation."
                
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.ib_t = self.ibR_t + 1j*self.ibI_t
            self.ic_t = self.icR_t + 1j*self.icI_t
        
            self.ub_t = self.ubR_t + 1j*self.ubI_t
            self.uc_t = self.ucR_t + 1j*self.ucI_t
            
            self.xb_t = self.xbR_t+1j*self.xbI_t
            self.xc_t = self.xcR_t+1j*self.xcI_t
            assert len(self.ua_t) == len(self.ub_t) == len(self.uc_t) == len(self.xa_t) == len(self.xb_t) == len(self.xc_t) != None, "States must be available from simulation."
            
        #Time series PCC voltage
        self.time_series_duty_cycle()
        
        #Exhibit different behavior in stand alone mode
        if self.PV_model.standAlone:
            self.time_series_Zload1()        
        self.time_series_grid()
        self.time_series_inv_terminal_voltage()
        self.time_series_PCC_LV_side_voltage()
        
        if self.PV_model.standAlone:
            self.time_series_PCC_HV_side_voltage()
        
        self.time_series_S()
        self.time_series_RMS()
        self.time_series_Ppv()
        
        self.time_series_phase_angle()
        #self.time_series_power_transfer()
        self.time_series_PLL()
        
        self.logger.debug("All states collected!")
    
    def collect_states(self,solution):
        """Collect states from ode solution."""
        
         #Phase a states
        self.iaR_t = solution[:,0]
        self.iaI_t = solution[:,1]
        self.xaR_t = solution[:,2]
        self.xaI_t = solution[:,3]
        self.uaR_t = solution[:,4]
        self.uaI_t = solution[:,5]              
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
            #DC link voltage variables
            self.Vdc_t = solution[:,6]
            self.xDC_t = solution[:,7]
            self.xQ_t = solution[:,8]
            #PLL variables
            self.xPLL_t = solution[:,9]
            #Frequency integration to get angle
            self.wte_t = solution[:,10]

        elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            
            #Phase b states
            self.ibR_t = solution[:,6]
            self.ibI_t = solution[:,7]
            self.xbR_t = solution[:,8]
            self.xbI_t = solution[:,9]
            self.ubR_t = solution[:,10]
            self.ubI_t = solution[:,11]

            #Phase c states
            self.icR_t = solution[:,12]
            self.icI_t = solution[:,13]
            self.xcR_t = solution[:,14]
            self.xcI_t = solution[:,15]
            self.ucR_t = solution[:,16]
            self.ucI_t = solution[:,17]

            #DC link voltage variables
            self.Vdc_t = solution[:,18]
            self.xDC_t = solution[:,19]
            self.xQ_t = solution[:,20]

            #PLL variables
            self.xPLL_t = solution[:,21]
            #Frequency integration to get angle
            self.wte_t = solution[:,22]        
    
    def invert_arrays(self):
        """Inverter arrays before usage by plots."""
        
        self.t_t = self._t_t
        self.Vdc_t = self._Vdc_t
        
        self.ia_t = self._ia_t
        self.ma_t = self._ma_t
        self.vta_t = self._vta_t
        self.va_t = self._va_t
        
        self.phita_t = self._phita_t
        self.phia_t = self._phia_t
        
        self.wgrid_t = self._wgrid_t
        self.we_t = self._we_t
        self.wte_t = self._wte_t
        self.vd_t = self._vd_t
        self.vq_t = self._vq_t
        
        self.ma_absolute_t = self._ma_absolute_t
        self.Varms_t = self._Varms_t
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.mb_absolute_t = self._mb_absolute_t
            self.mc_absolute_t = self._mc_absolute_t
            
            self.Vbrms_t = self._Vbrms_t
            self.Vcrms_t = self._Vcrms_t
            
            self.ib_t = self._ib_t
            self.mb_t = self._mb_t
            self.vtb_t = self._vtb_t
            self.vb_t = self._vb_t
            
            self.phitb_t = self._phitb_t
            self.phib_t = self._phib_t
            
            self.ic_t = self._ic_t
            self.mc_t = self._mc_t
            self.vtc_t = self._vtc_t
            self.vc_t = self._vc_t
            
            self.phitc_t = self._phitc_t
            self.phic_t = self._phic_t
                
        self.Vtrms_t = self._Vtrms_t
        self.Vrms_t = self._Vrms_t
                
        self.Irms_t = self._Irms_t
        
        self.Ppv_t  = self._Ppv_t 
        self.S_t =  self._S_t 
        self.S_PCC_t = self._S_PCC_t        
        
    def run_simulation(self,gridVoltagePhaseA=None,gridVoltagePhaseB=None,gridVoltagePhaseC=None,gridFrequency=None,y0=None,t=None):
        """Call the ODE solver and collect states."""
        
        self.solution_time = None #Always reset simulation time to None
        
        if self.LOOP_MODE:
            if isinstance(gridVoltagePhaseA,complex) and isinstance(y0,list) and isinstance(t,list):
                self.update_grid_measurements(gridVoltagePhaseA=gridVoltagePhaseA, gridVoltagePhaseB=gridVoltagePhaseB, gridVoltagePhaseC=gridVoltagePhaseC,gridFrequency=gridFrequency)
        
            else:
                raise ValueError('Grid voltage (complex scalar), initial states (y0), and time steps (t) should be provided in loop mode.')
            
            if t[0] == 0.0:
                self.t = t
                six.print_("{}:Simulation started in loop mode with a step size of {:.4f} s!".format(self.name,self.t[-1]-self.t[0]))
                if self.jacFlag:
                    self.logger.debug("{}:Analytical Jacobian will be provided to ODE solver.".format(self.name))
            solution,_,_  = self.call_ODE_solver(self.ODE_model,self.jac_ODE_model,y0,t)
            
        else:
            self.t = self.t_calc()
            self.initialize_y0_t()
            
            timer_start = time.time()
            six.print_("{}:Simulation started at {} s and will end at {} s".format(self.name,self.tStart,self.tStop))
            self.logger.debug('{}:{} solver will be used.'.format(self.name,self.solver_type))
            if self.jacFlag:
                self.logger.debug("{}:Analytical Jacobian will be provided to ODE solver.".format(self.name))
                
            if self.solver_type != 'odeint':
                self.ode_solver.set_initial_value(self.y0,self.tStart) 
                self.logger.debug("{}:Resetting {} internal time step to {} and states to:\n{}.".format(self.name,self.solver_type,self.tStart,self.y0))
        
            solution,_,_  = self.call_ODE_solver(self.ODE_model,self.jac_ODE_model,self.y0,self.t)
            
            self.solution_time = time.time() - timer_start
            
            self.show_simulation_time()
           
            self.simulation_events.reset_event_counters() #Reset grid and solar event counters
            self.PV_model.reset_reference_counters() #Reset counters for reference change events
        
        if self.COLLECT_SOLUTION:
            self.collect_solution(solution,t=t)  #Collect full solution
        else:
            self.collect_states(solution)  #Atleast states must be collected    
    
    def update_grid_measurements(self,gridVoltagePhaseA, gridVoltagePhaseB, gridVoltagePhaseC,gridFrequency):
        """Update grid voltage and frequency in non-standalone model.
        Args:
             gridVoltagePhaseA (complex): Value of gridVoltagePhaseA
             gridVoltagePhaseB (complex): Value of gridVoltagePhaseB
             gridVoltagePhaseC (complex): Value of gridVoltagePhaseC             
        
        """
        
        if isinstance(gridFrequency,float):
            self.PV_model.gridFrequency = gridFrequency
        else:
            self.PV_model.gridFrequency = self.PV_model.gridFrequency
        
        self.PV_model.gridVoltagePhaseA = gridVoltagePhaseA
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.PV_model.gridVoltagePhaseB = gridVoltagePhaseB
            self.PV_model.gridVoltagePhaseC = gridVoltagePhaseC
        #elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhaseBalanced':
        #    self.PV_model.gridVoltagePhaseB = utility_functions.Ub_calc(gridVoltagePhaseA)
        #    self.PV_model.gridVoltagePhaseC = utility_functions.Uc_calc(gridVoltagePhaseA)
    
    def playback_3phgrid(self,Va_t,Vb_t,Vc_t,fgrid_t,n_time_steps=None):
        """Playback grid voltages and frequency.
         Args:
             Va_t (list): Time series of phase A grid voltages in complex phasor form (V).
             Vb_t (list): Time series of phase B grid voltages in complex phasor form (V).
             Vc_t (list): Time series of phase C grid voltages in complex phasor form (V).
             fgrid_t (list): Time series of frequency (Hz).             
        """
        
        assert self.LOOP_MODE and not self.PV_model.standAlone, 'Simulation must be in loop mode.'
        
        if not isinstance(n_time_steps,int):
            n_time_steps = len(Va_t) #Number of time steps
        
        print('Found {} elements in phase A grid voltage time series.'.format(n_time_steps))
        
        if isinstance(fgrid_t,float):
            fgrid_t = [fgrid_t]*n_time_steps
            print('Creating a list of constant frequency values.')
        
        assert len(Va_t) ==len(Vb_t) ==len(Vc_t) ==len(fgrid_t), 'All lists should have equal number of elements.'
        
        t0 = 0.0
        dt = self.tInc
        
        for i in range(n_time_steps):
            Va = Va_t[i]/self.Vbase
            Vb = Vb_t[i]/self.Vbase
            Vc = Vc_t[i]/self.Vbase
            wgrid = 2*math.pi*fgrid_t[i]
            t_sim= [t0,t0+dt]
            
            self.run_simulation(gridVoltagePhaseA=Va, gridVoltagePhaseB=Vb, gridVoltagePhaseC=Vc,gridFrequency=wgrid, y0=self.y0, t=t_sim)
            
            print('t:{:.3f},wgrid:{:.2f},winv:{:.2f},Va:{:.2f},Vb:{:.2f},Vc:{:.2f}'.format(t0,wgrid,self.PV_model.winv,Va,Vb,Vc))
            print('Va_actual:{:.2f}'.format(self.PV_model.va))
            print('_Va_t_actual:{:.2f},Va_t_actual:{:.2f}'.format(self._va_t[-1],self.va_t[-1]))
            print('wgrid_actual:{:.2f},wgrid_measured:{:.2f},_wgrid_t:{:.2f}'.format(wgrid,self.PV_model.gridFrequency,self._wgrid_t[-1]))
            print('Size:',len(self._va_t))
            
            t0 = t_sim[-1]

    def show_simulation_time(self):
        """Show simulation time."""
        
        print('{}:Simulation was completed in {}'.format(self.name,time.strftime("%H:%M:%S", time.gmtime(self.solution_time))))
        
    def get_trajectories(self):
        """Return trajectories as a dictionary."""
        
        if self.LOOP_MODE:
            self.invert_arrays()
        
        trajectory_dictionary = {'ia_t':self.ia_t,'ma_t':self.ma_t,
                                 'Vdc_t':self.Vdc_t,'vta_t':self.vta_t,'va_t':self.va_t,
                                 'Irms_t':self.Irms_t,'Vtrms_t':self.Vtrms_t,'Vrms_t':self.Vrms_t,
                                 'Ppv_t':self.Ppv_t,'S_t':self.S_t,'S_PCC_t':self.S_PCC_t
                                }
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            
            trajectory_dictionary.update({'ib_t':self.ib_t,'ic_t':self.ic_t,'mb_t':self.mb_t,'mc_t':self.mc_t,
                                          'vtb_t':self.vtb_t,'vtc_t':self.vtc_t,
                                          'vb_t':self.vb_t,'vc_t':self.vc_t})     
         
        for trajectory in trajectory_dictionary.values():
            
            assert len(self.t_t) == len(trajectory), 'The length of each trajectory should be equal to the time_steps.'
        
        return trajectory_dictionary