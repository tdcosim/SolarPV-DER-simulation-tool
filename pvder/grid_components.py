"""Grid model and shared attributes."""

from __future__ import division
import numpy as np
import math
import cmath
import six
from pvder import utility_functions

class BaseValues():
	"""Class to store base values."""
	
	Vbase = 500.0  #L-G peak"
	Sbase = 50e3 #VA base
	wbase = 2*math.pi*60.0
	Vdcbase = Vbase #DC side base value is same as AC side base value
	Ibase = Sbase/Vbase
	Zbase = (Vbase**2)/Sbase
	
	Lbase = Zbase/wbase
	Cbase = 1/(Zbase*wbase)
	
class Grid(BaseValues):
	""" Class for grid"""
	
	grid_count = 0 #Count for grid objects
	
	n_ODE = 0 #Number of ODE's
	
	Vgridrated =  20415.0 # L-G peak to peak equivalent to 25000 V L-L RMS
	_t_voltage_previous = 0.0
	_t_frequency_previous = 0.0
	
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
		self.name = 'grid_'+str(Grid.grid_count)
		
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
		self.Vagrid,self.wgrid = events.grid_events(t=0.0) #Grid voltage and frequency set-point
		self.Vagrid = self.Vagrid*(self.Vgridrated/self.Vbase)
		self.Vbgrid = utility_functions.Ub_calc(self.Vagrid*self.unbalance_ratio_b)
		self.Vcgrid = utility_functions.Uc_calc(self.Vagrid*self.unbalance_ratio_c)
		#Actual Grid voltage
		self.vag = self.Vagrid
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
  
	def steady_state_model(self,t):
		"""Grid voltage change."""
		
		Vagrid_new,wgrid_new = self.events.grid_events(t)
		Vagrid_new = Vagrid_new*(self.Vgridrated/self.Vbase)
		
		if abs(self.Vagrid- Vagrid_new) > 0.0 and t >= self._t_voltage_previous:
			utility_functions.print_to_terminal("{}:Grid voltage changed from {:.3f} V to {:.3f} V at {:.3f} s".format(self.name,self.Vagrid,Vagrid_new,t))
			
			self.Vagrid = Vagrid_new
			self.Vbgrid = utility_functions.Ub_calc(self.Vagrid*self.unbalance_ratio_b)
			self.Vcgrid = utility_functions.Uc_calc(self.Vagrid*self.unbalance_ratio_c)
			
			self._t_voltage_previous = t
			
		if abs(self.wgrid- wgrid_new) > 0.0 and t >= self._t_frequency_previous:
			utility_functions.print_to_terminal("{}:Grid frequency changed from {:.3f} Hz to {:.3f} Hz at {:.3f} s".format(self.name,self.wgrid/(2.0*math.pi),wgrid_new/(2.0*math.pi),t))
					  
			self.wgrid = wgrid_new
			self._t_frequency_previous = t			
		
		self.vag = self.Vagrid
		self.vbg = self.Vbgrid
		self.vcg = self.Vcgrid					
   