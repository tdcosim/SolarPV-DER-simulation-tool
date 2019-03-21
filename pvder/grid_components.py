from __future__ import division
import numpy as np
import math
import cmath
import six
from utilities import SimulationUtilities
import utility_functions as utility_functions

class Grid(SimulationUtilities):
    """ Class for grid"""
    
    grid_count = 0
    #Number of ODE's
    n_grid_ODE = 6
    
    #Grid voltage time constant
    Vgridrated =  20415.0 # L-G peak to peak equivalent to 25000 V L-L RMS
    Tgrid = 0.000001
    Vbase = 500.0  #L-G peak"
    Sbase = 50e3 #VA base
    wbase = 2*math.pi*60.0
    
    Ibase = Sbase/Vbase
    Zbase = (Vbase**2)/Sbase
    
    Lbase = Zbase/wbase
    Cbase = 1/(Zbase*wbase)
    
    #Simulation time steps
    tStart = 0.0
    tStop = 0.5
    tInc = 0.001
    t = np.arange(tStart, tStop, tInc)
    
    def __init__(self,events,unbalance_ratio_b=1.0,unbalance_ratio_c=1.0,Z2_actual = 1.61 + 1j*5.54):
        """Creates an instance of `GridSimulation`.
        
        Args:
          events: An instance of `SimulationEvents`.
          unbalance_ratio_b,unbalance_ratio_c: Scalar specifying difference in Phase B and Phase C voltage magnitude compared to phase A.
          Z2_actual: Complex scalar specifying the impedance of the feeder connecting the DER with the voltage source.
        """
        
        #Increment count
        Grid.grid_count = Grid.grid_count+1
        #Events object
        self.events = events
        
        #Object name
        self.name = 'grid'+str(Grid.grid_count)
        
        #Voltage unbalance
        self.unbalance_ratio_b = unbalance_ratio_b
        self.unbalance_ratio_c = unbalance_ratio_c
        
        #Grid impedance
        self.Z2_actual = Z2_actual
        self.R2_actual = self.Z2_actual.real
        self.L2_actual = self.Z2_actual.imag/(2*math.pi*60.0)
        
        #Converting to per unit
        self.R2 = self.R2_actual/self.Zbase  #Transmission line resistance
        self.L2 = self.L2_actual/self.Lbase  #Transmission line resistance

        self.Z2 = self.Z2_actual/self.Zbase   #Transmission line impedance

        self.transmission_name = 'transmission_'+str(Grid.grid_count)
        #Grid voltage/frequency events
        self.Vagrid_no_conversion,self.wgrid = events.grid_events(t=0.0) #Grid voltage and frequency set-point
        self.Vagrid_no_conversion = self.Vagrid_no_conversion*(self.Vgridrated/self.Vbase)
        #Grid voltage setpoint
        self.Vagrid = self.grid_voltage_setpoint_conversion(self.Vagrid_no_conversion,self.wgrid)
        #self.Vagrid = self.Vagrid_no_conversion
        self.Vbgrid = utility_functions.Ub_calc(self.Vagrid*self.unbalance_ratio_b)
        self.Vcgrid = utility_functions.Uc_calc(self.Vagrid*self.unbalance_ratio_c)
        #Actual Grid voltage
        self.vag = self.Vagrid_no_conversion
        self.vbg = utility_functions.Ub_calc(self.vag*self.unbalance_ratio_b)
        self.vcg = utility_functions.Uc_calc(self.vag*self.unbalance_ratio_c)
        self.Vgrms = self.Vgrms_calc()
    
    @property
    def y0(self):        
        """Grid states"""
        return [self.vag.real,self.vag.imag,\
                self.vbg.real,self.vbg.imag,\
                self.vcg.real,self.vcg.imag]

    def Vgrms_calc(self):
        """Grid side terminal voltage -  RMS"""
        return utility_functions.Urms_calc(self.vag,self.vbg,self.vcg)
        
    def update_grid_states(self,vag,vbg,vcg):
        """Update grid states"""

        self.vag = vag
        self.vbg = vbg
        self.vcg = vcg

    def ODE_model(self,y,t):
        """ODE's for grid voltage."""
        
        vagR, vagI, vbgR, vbgI, vcgR, vcgI = y
        
        self.update_grid_states(vagR+1j*vagI,vbgR+1j*vbgI,vcgR+1j*vcgI)
        #Grid conditions
        Vagrid_new,wgrid_new = self.events.grid_events(t)
        Vagrid_new = Vagrid_new*(self.Vgridrated/self.Vbase)

        if abs(self.Vagrid_no_conversion- Vagrid_new) or abs(self.wgrid- wgrid_new) > 0.0:
           utility_functions.print_to_terminal("Grid voltage changed from {:.3f} V to {:.3f} V at {:.3f}".format(self.Vagrid_no_conversion,Vagrid_new,t))
           utility_functions.print_to_terminal("Grid frequency changed from {:.3f} Hz to {:.3f} Hz at {:.3f}".format(self.wgrid/(2.0*math.pi),wgrid_new/(2.0*math.pi),t))
           self.Vagrid_no_conversion = Vagrid_new#*self.Vbase
           self.wgrid = wgrid_new
           #Conversion of grid voltage setpoint
           self.Vagrid = self.grid_voltage_setpoint_conversion(Vagrid_new,wgrid_new)      
           self.Vbgrid = utility_functions.Ub_calc(self.Vagrid*self.unbalance_ratio_b)
           self.Vcgrid = utility_functions.Uc_calc(self.Vagrid*self.unbalance_ratio_c)
        
        self.Vgrms = self.Vgrms_calc()
        #Dynamics for grid voltage
        
        dvagR = (1/self.Tgrid)*(self.Vagrid.real - vagR)   #Tgrid = time constant
        dvagI = (1/self.Tgrid)*(self.Vagrid.imag - vagI)   #Tgrid = time constant
        dvbgR = (1/self.Tgrid)*(self.Vbgrid.real - vbgR)   #Tgrid = time constant
        dvbgI = (1/self.Tgrid)*(self.Vbgrid.imag - vbgI)   #Tgrid = time constant
        dvcgR = (1/self.Tgrid)*(self.Vcgrid.real - vcgR)   #Tgrid = time constant
        dvcgI = (1/self.Tgrid)*(self.Vcgrid.imag - vcgI)   #Tgrid = time constant
        
        #utility_functions.print_to_terminal("t:{:.4f},Vagrid {:.3f}, vag {:.3f} V".format(t,self.Vagrid,self.vag))
        result=[dvagR,
                dvagI,
                dvbgR,
                dvbgI,
                dvcgR,
                dvcgI]

        return result

    def jac_ODE_model(self,y,t):
        """ODE's for grid voltage."""
        
        J=np.zeros((29,29))
        
        vagR, vagI, vbgR, vbgI, vcgR, vcgI = y
        
        self.update_grid_states(vagR+1j*vagI,vbgR+1j*vbgI,vcgR+1j*vcgI)
        #Grid conditions
        Vagrid_new,wgrid_new = self.events.grid_events(t)
        Vagrid_new = Vagrid_new*(self.Vgridrated/self.Vbase)

        if abs(self.Vagrid_no_conversion- Vagrid_new) or abs(self.wgrid- wgrid_new) > 0.0:
           utility_functions.print_to_terminal("Grid voltage changed from {:.3f} V to {:.3f} V at {:.3f}".format(self.Vagrid_no_conversion,Vagrid_new,t))
           utility_functions.print_to_terminal("Grid frequency changed from {:.3f} Hz to {:.3f} Hz at {:.3f}".format(self.wgrid/(2.0*math.pi),wgrid_new/(2.0*math.pi),t))
           self.Vagrid_no_conversion = Vagrid_new#*self.Vbase
           self.wgrid = wgrid_new
           #Conversion of grid voltage setpoint
           self.Vagrid = self.grid_voltage_setpoint_conversion(Vagrid_new,wgrid_new)      
           self.Vbgrid = utility_functions.Ub_calc(self.Vagrid*self.unbalance_ratio_b)
           self.Vcgrid = utility_functions.Uc_calc(self.Vagrid*self.unbalance_ratio_c)
        
        self.Vgrms = self.Vgrms_calc()
        #Dynamics for grid voltage
        
        J[23+0,23+0]=-(1/self.Tgrid)
        J[23+1,23+1]=-(1/self.Tgrid)
        J[23+2,23+2]=-(1/self.Tgrid)
        J[23+3,23+3]=-(1/self.Tgrid)
        J[23+4,23+4]=-(1/self.Tgrid)
        J[23+5,23+5]=-(1/self.Tgrid)
 
        return J
        
    def grid_voltage_setpoint_conversion(self,Vgrid_new,wgrid_new):
       
        #return (1+1j*(wgrid_new/self.wbase)*self.Tgrid)*Vgrid_new
        return Vgrid_new