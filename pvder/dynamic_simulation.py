from __future__ import division
import numpy as np
import math
import cmath
import time

import pdb
import six

####from graphviz import Digraph
from pvder.grid_components import Grid
from pvder.utilities import SimulationUtilities
from pvder import utility_functions

class ModelUtilities():
    """Class for model wide utilities."""
    
    def __init__(self,grid_model,PV_model,simulation):
        self.PV_model = PV_model
        self.grid_model = grid_model
        self.simulation = simulation
    
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

class GridSimulation(Grid,SimulationUtilities):
    """ Utility class for running simulations."""
    
    sim_count = 0
    tStart = 0.0
    jacFlag = False
    DEBUG_SOLVER = False
    DEBUG_SIMULATION = False
    DEBUG_CONTROLLERS = False
    DEBUG_VOLTAGES = False
    DEBUG_CURRENTS = False
    DEBUG_POWER = False
    DEBUG_PLL = True
    
    def __init__(self,PV_model,events,grid_model=None,tStop = 0.5,tInc = 0.001):
        """Creates an instance of `GridSimulation`.
        
        Args:
          PV_model: An instance of `SolarPV_DER`.
          events: An instance of `SimulationEvents`.
          grid_model: An instance of `GridModel` (only need to be suppled in stand alone simulation).
          tStop: A scalar specifying the end time for simulation.
          tInc: A scalar specifying the time step for simulation.
        """
        
        #Increment count to keep track of number of simulation instances
        GridSimulation.sim_count = GridSimulation.sim_count+1
        self.sim_ID = GridSimulation.sim_count
        #Object name
        self.name = 'sim_'+str(self.sim_ID)
        
        self.tStop = tStop
        self.tInc = tInc
        self.PV_model = PV_model
        self.simulation_events = events
        if self.PV_model.standAlone and grid_model is not None:
            self.grid_model = grid_model
        elif self.PV_model.standAlone and grid_model is None:
            raise ValueError('`Grid` instance need to provided in stand alone mode for creating `GridSimulation` instance!')
        
        #Remove existing simulation events
        self.simulation_events.remove_solar_event(3.0)
        self.simulation_events.remove_load_event(4.0)
        self.simulation_events.remove_grid_event(5.0)
    
    @property
    def y0(self):
        """ Combine all initial conditions."""
        """
        if not self.PV_model.standAlone:
           y0 = self.PV_model.y0
        else:
           y0 = self.PV_model.y0+ self.grid_model.y0
        """
        y0 = self.PV_model.y0
        return y0
    
    @property
    def t(self):
        """Vector of time steps for simulation"""
        return np.arange(self.tStart, self.tStop, self.tInc)

    def ODE_model(self,y,t):
        """ Combine all derivatives."""
        
        y1 = y[0:self.PV_model.n_ODE]
        
        if self.PV_model.standAlone:
            self.grid_model.grid_model(t)
            
        y = self.PV_model.ODE_model(y1,t)
        
        """
        if not self.PV_model.standAlone:
           y1 = y[0:self.PV_model.n_ODE]
           y = self.PV_model.ODE_model(y1,t) 
        
        else:
           #y1 = y[0:self.PV_model.n_ODE]
           #y2 = y[self.PV_model.n_ODE:]
           #y = self.PV_model.ODE_model(y1,t) + self.grid_model.ODE_model(y2,t)
           
        """ 
        if self.DEBUG_SIMULATION:
            self.debug_simulation(t)

        return y

        
    def jac_ODE_model(self,y,t):
        """ Combine all derivatives."""
        
        y1 = y[0:self.PV_model.n_ODE]
        
        if self.PV_model.standAlone:
            self.grid_model.grid_model(t)
            
        y = self.PV_model.jac_ODE_model(y1,t)
        """
        if not self.PV_model.standAlone:
            y1 = y[0:self.PV_model.n_ODE]
            y = self.PV_model.jac_ODE_model(y1,t)
        else:
            y1 = y[0:self.PV_model.n_ODE]
            y2 = y[self.PV_model.n_ODE:]
            y = self.PV_model.jac_ODE_model(y1,t) + self.grid_model.jac_ODE_model(y2,t)
        """
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
            utility_functions.print_to_terminal('xdc:{:.3f},xQ:{:.3f},ua:{:.3f},xa:{:.3f},ma:{:.3f}'. format(self.PV_model.xdc,self.PV_model.xQ,self.PV_model.ua,self.PV_model.xa,self.PV_model.ma))
            
        if self.DEBUG_PLL:
            utility_functions.print_to_terminal("we:{:.3f}, wte:{:.3f} rad, vdg: {:.3f} V, vqg {:.3f} V, vd: {:.3f} V, vq {:.3f} V".format(self.PV_model.we,self.PV_model.wte,self.PV_model.vdg,self.PV_model.vqg,self.PV_model.vd,self.PV_model.vq))
    
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
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.mb_t =  utility_functions.m_time_series(self.ub_t,self.xb_t,self.PV_model.Kp_GCC)
            self.mc_t =  utility_functions.m_time_series(self.uc_t,self.xc_t,self.PV_model.Kp_GCC)

            self.mbR_t = self.mb_t.real 
            self.mbI_t = self.mb_t.imag
            self.mcR_t = self.mc_t.real 
            self.mcI_t = self.mc_t.imag

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
    
        elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.Vtrms_t = utility_functions.Urms_time_series(self.vta_t,self.vtb_t,self.vtc_t)
            self.Vrms_t = utility_functions.Urms_time_series(self.va_t,self.vb_t,self.vc_t)
            self.Irms_t = utility_functions.Urms_time_series(self.ia_t,self.ib_t,self.ic_t)       
       
        self.Vgrms_t = utility_functions.Urms_time_series(self.vag_t,self.vbg_t,self.vcg_t)
        if self.PV_model.standAlone:
            self.Vhvrms_t= utility_functions.Urms_time_series(self.vaHV_t,self.vbHV_t,self.vcHV_t)
        else:
            self.Vhvrms_t= self.Vgrms_t        
        
    def time_series_vgrid(self):
        """Time series grid frequency."""
        
        self.vag_t = []
        self.vbg_t = []
        self.vcg_t = []
        for i,t in enumerate(self.t):   #Loop through grid events and calculate wgrid at each time step

                Vagrid_new,_ = self.simulation_events.grid_events(t)
                
                #Conversion of grid voltage setpoint
                self.grid_model.vag = Vagrid_new*(self.grid_model.Vgridrated/self.Vbase)
                self.grid_model.vbg = utility_functions.Ub_calc(self.grid_model.vag*self.grid_model.unbalance_ratio_b)
                self.grid_model.vcg = utility_functions.Uc_calc(self.grid_model.vag*self.grid_model.unbalance_ratio_c)  
        
                self.vag_t.append(self.grid_model.vag)
                self.vbg_t.append(self.grid_model.vbg)
                self.vcg_t.append(self.grid_model.vcg)

        self.vag_t = np.asarray(self.vag_t)
        self.vbg_t = np.asarray(self.vbg_t)
        self.vcg_t = np.asarray(self.vcg_t)
        
        self.simulation_events.reset_event_counters() #reset event counters
    
    def time_series_wgrid(self):
        """Time series grid frequency."""
        
        self.wgrid_t = []
        for i,t in enumerate(self.t):   #Loop through grid events and calculate wgrid at each time step
            if self.PV_model.standAlone:
                _,self.grid_model.wgrid =  self.simulation_events.grid_events(t)
                self.wgrid_t.append(self.grid_model.wgrid)
            else:
                self.wgrid_t.append(self.PV_model.wgrid_measured)  
        self.wgrid_t = np.asarray(self.wgrid_t)
        self.simulation_events.reset_event_counters() #reset event counters
    
    def time_series_d_q(self):
        """Calculate time series d-q quantities."""
        
        self.vdg_t = []
        self.vqg_t = []
        
        self.vd_t = []
        self.vq_t = []
        
        for i,t in enumerate(self.t):     #Loop through time steps and calculate d-q values at each time step
            
            self.PV_model.vagt = utility_functions.phasor_to_time_domain(self.vag_t[i],self.wgrid_t[i],t)
            self.PV_model.vbgt = utility_functions.phasor_to_time_domain(self.vbg_t[i],self.wgrid_t[i],t)
            self.PV_model.vcgt = utility_functions.phasor_to_time_domain(self.vcg_t[i],self.wgrid_t[i],t)
            
            self.PV_model.vat = utility_functions.phasor_to_time_domain(self.va_t[i],self.wgrid_t[i],t)
            self.PV_model.vbt = utility_functions.phasor_to_time_domain(self.vb_t[i],self.wgrid_t[i],t)
            self.PV_model.vct = utility_functions.phasor_to_time_domain(self.vc_t[i],self.wgrid_t[i],t)
            
            self.PV_model.wte = self.wte_t[i]
            
            self.PV_model.vdg,self.PV_model.vqg,self.PV_model.v0g = utility_functions.abc_to_dq0(self.PV_model.vagt,self.PV_model.vbgt,self.PV_model.vcgt,self.PV_model.wte)
            self.vdg_t.append(self.PV_model.vdg)
            self.vqg_t.append(self.PV_model.vqg)
            
            self.PV_model.vd,self.PV_model.vq,self.PV_model.v0 = utility_functions.abc_to_dq0(self.PV_model.vat,self.PV_model.vbt,self.PV_model.vct,self.PV_model.wte)
            self.vd_t.append(self.PV_model.vd)
            self.vq_t.append(self.PV_model.vq)
            
        
        self.vdg_t = np.asarray(self.vdg_t)
        self.vqg_t = np.asarray(self.vqg_t)
        
        self.vd_t = np.asarray(self.vd_t)
        self.vq_t = np.asarray(self.vq_t)

    def time_series_PLL(self):
        """Calculate time series PLL output frequency."""
        
        self.we_t = []
        for i,t in enumerate(self.t):       #Loop through time steps and calculate we at each time step
            self.PV_model.vdg = self.vdg_t[i]
            self.PV_model.vd = self.vd_t[i]
            self.PV_model.xPLL = self.xPLL_t[i]
            self.we_t.append(self.PV_model.we_calc())
        self.we_t = np.asarray(self.we_t)
    
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
        self.simulation_events.reset_event_counters() #reset event counters
    
    def time_series_S(self):
        """Calculate time series apparent power."""
        
        if type(self.PV_model).__name__ == 'SolarPV_DER_SinglePhase':
            self.S_t = (1/2)*(self.vta_t*self.ia_t.conjugate())
            self.S_PCC_t = (1/2)*(self.va_t*self.ia_t.conjugate())
            self.S_load1_t = (1/2)*(self.va_t*(-(self.va_t/self.Zload1_t)).conjugate())
            
            if self.PV_model.standAlone:
                self.S_G_t = (1/2)*(-(self.ia_t - (self.va_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vag_t
                
        elif type(self.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
            self.S_t = (1/2)*(self.vta_t*self.ia_t.conjugate() + self.vtb_t*self.ib_t.conjugate() + self.vtc_t*self.ic_t.conjugate())
            self.S_PCC_t = (1/2)*(self.va_t*self.ia_t.conjugate()+self.vb_t*self.ib_t.conjugate()+self.vc_t*self.ic_t.conjugate())
            self.S_load1_t = (1/2)*(self.va_t*(-(self.va_t/self.Zload1_t)).conjugate() + self.vb_t*(-(self.vb_t/self.Zload1_t)).conjugate() + self.vc_t*(-(self.vc_t/self.Zload1_t)).conjugate())
            
            if self.PV_model.standAlone:
                self.S_G_t = (1/2)*(-(self.ia_t - (self.va_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vag_t + (-(self.ib_t -(self.vb_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vbg_t + (-(self.ic_t -(self.vc_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vcg_t
    
    def time_series_phase_angle(self):
        """Calculate time series phase angle between grid and inverter voltage."""
        assert len(self.va_t) == len(self.va_t) == len(self.va_t) == len(self.vag_t) == len(self.vbg_t) == len(self.vcg_t) != None, "States must be available from simulation."
        
        self.phi_at_t = np.angle(self.vta_t)
        self.phi_bt_t = np.angle(self.vtb_t)
        self.phi_ct_t = np.angle(self.vtc_t)
        
        self.phi_a_t = np.angle(self.va_t)
        self.phi_b_t = np.angle(self.vb_t)
        self.phi_c_t = np.angle(self.vc_t)
        
        self.phi_aHV_t = np.angle(self.vaHV_t)
        self.phi_bHV_t = np.angle(self.vbHV_t)
        self.phi_cHV_t = np.angle(self.vcHV_t)
        
        self.phi_ag_t = np.angle(self.vag_t)
        self.phi_bg_t = np.angle(self.vbg_t)
        self.phi_cg_t = np.angle(self.vcg_t)
    
    def time_series_power_transfer(self):
        """Calculate time series power transfer between grid and PCC."""
        self.P_transfer = (self.Vhvrms_t*self.Vgrms_t*np.sin(self.phi_aHV_t-self.phi_ag_t))/np.abs(self.grid_model.Z2)
        
    def run_simulation(self):
        """Call the ODE solver and collect states."""
        
        solution  = self.solve_ODE()
        
        self.simulation_events.reset_event_counters() #Reset grid and solar event counters
        #Collect states
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
        
        #Grid voltage
        if self.PV_model.standAlone:
            self.time_series_vgrid()
            """
            self.vagR_t = solution[:,self.PV_model.n_ODE+self.grid_model.n_ODE]
            self.vagI_t = solution[:,self.PV_model.n_ODE+self.grid_model.n_ODE+1]
            self.vbgR_t = solution[:,self.PV_model.n_ODE+self.grid_model.n_ODE+2]
            self.vbgI_t = solution[:,self.PV_model.n_ODE+self.grid_model.n_ODE+3]
            self.vcgR_t = solution[:,self.PV_model.n_ODE+self.grid_model.n_ODE+4]
            self.vcgI_t = solution[:,self.PV_model.n_ODE+self.grid_model.n_ODE+5]
            """
        else:
            self.vagR_t = np.repeat(self.PV_model.gridVoltagePhaseA.real,len(self.t))
            self.vagI_t = np.repeat(self.PV_model.gridVoltagePhaseA.imag,len(self.t))
            self.vbgR_t = np.repeat(self.PV_model.gridVoltagePhaseB.real,len(self.t))
            self.vbgI_t = np.repeat(self.PV_model.gridVoltagePhaseB.imag,len(self.t))
            self.vcgR_t = np.repeat(self.PV_model.gridVoltagePhaseC.real,len(self.t))
            self.vcgI_t = np.repeat(self.PV_model.gridVoltagePhaseC.imag,len(self.t))
            
            #Time series grid voltage
            self.vag_t = self.vagR_t + 1j*self.vagI_t
            self.vbg_t = self.vbgR_t + 1j*self.vbgI_t
            self.vcg_t = self.vcgR_t + 1j*self.vcgI_t
        
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
        self.time_series_Zload1()
        self.time_series_duty_cycle()
        
        self.time_series_inv_terminal_voltage()
        self.time_series_PCC_LV_side_voltage()
        if self.PV_model.standAlone:
            self.time_series_PCC_HV_side_voltage()
        
        self.time_series_S()
        self.time_series_RMS()
        self.time_series_Ppv()
        
        #self.time_series_phase_angle()
        #self.time_series_power_transfer()
        self.time_series_wgrid()
        #self.time_series_d_q()
        #self.time_series_PLL()
        six.print_("All states collected!")


