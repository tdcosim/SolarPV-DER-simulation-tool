"""Code for setting up, and running, and collecting data from PV-DER simulations."""

from __future__ import division
import numpy as np
import math
import cmath
import time

import pdb
import six

from pvder.utility_classes import Utilities
from pvder.grid_components import Grid
from pvder.simulation_utilities import SimulationUtilities
#from pvder.simulation_utilities_experimental import SimulationUtilitiesExperimental
from pvder import utility_functions
from pvder import defaults,templates
from pvder.logutil import LogUtil


class DynamicSimulation(Grid,SimulationUtilities,Utilities):
	""" Utility class for running simulations."""
	
	count = 0
	tStart = 0.0
	tInc = defaults.DEFAULT_DELTA_T
	
	DEBUG_SOLVER = False
	DEBUG_SIMULATION = False
	DEBUG_CONTROLLERS = False
	DEBUG_VOLTAGES = False
	DEBUG_CURRENTS = False
	DEBUG_POWER = False
	DEBUG_PLL = False
	
	jac_list = ['SolarPVDERThreePhase','SolarPVDERSinglePhase','SolarPVDERThreePhaseBalanced']
		
	def __init__(self,PV_model,events,gridModel = None,tStop = 0.5,
				 LOOP_MODE = False,COLLECT_SOLUTION = True,jacFlag = False,
				 verbosity ='INFO',solverType ='odeint',identifier = None):
		"""Creates an instance of `GridSimulation`.
		Args:
		  PV_model: An instance of `SolarPV_DER`.
		  events: An instance of `SimulationEvents`.
		  grid_model: An instance of `GridModel` (only need to be suppled in stand alone simulation).
		  tStop: A scalar specifying the end time for simulation.
		  tInc: A scalar specifying the time step for simulation.
		  LOOP_MODE: A boolean specifying whether simulation is run in loop.
		"""
		try:
			DynamicSimulation.count = DynamicSimulation.count + 1 #Increment count to keep track of number of simulation instances
			self.name_instance(identifier) #Generate a name for the instance
			self.tStop = tStop
			self.t = self.t_calc()
			self.PV_model = PV_model
			self.DER_model_type = type(self.PV_model).__name__
			self.simulation_events = events
			self.simulation_events.del_t_event = self.tInc
			self.initialize_solver(solver_type=solverType)
			self.SOLVER_CONVERGENCE = False
			self.convergence_failure_list =[]
			self.LOOP_MODE = LOOP_MODE
			self.COLLECT_SOLUTION = COLLECT_SOLUTION
			self.jacFlag = jacFlag
			self.check_jac_availability()
		
			if self.PV_model.standAlone and gridModel is not None:
				self.grid_model = gridModel
			elif self.PV_model.standAlone and gridModel is None:
				raise ValueError('`Grid` instance need to provided in stand alone mode for creating `GridSimulation` instance!')
		
			#Remove existing simulation events
			#self.simulation_events.remove_solar_event(3.0)
			#self.simulation_events.remove_load_event(4.0)
			#self.simulation_events.remove_grid_event(5.0)
			self.solution_time = None #Always reset solution time to None
			if self.LOOP_MODE:
				self.reset_stored_trajectories()
			self.initialize_y0_t()
		except:
			LogUtil.exception_handler()


	@property
	def y0(self):
		""" Combine all initial conditions from solution."""
		try:
			if type(self.PV_model).__name__ == 'SolarPVDERThreePhase':
				y0 = [self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],\
					  self.ibR_t[-1], self.ibI_t[-1], self.xbR_t[-1], self.xbI_t[-1], self.ubR_t[-1],self.ubI_t[-1],\
					  self.icR_t[-1], self.icI_t[-1], self.xcR_t[-1], self.xcI_t[-1], self.ucR_t[-1],self.ucI_t[-1],\
					  self.Vdc_t[-1],
					  self.xDC_t[-1],self.xQ_t[-1],
					  self.xPLL_t[-1],self.wte_t[-1]]	
		
			elif type(self.PV_model).__name__ == 'SolarPVDERSinglePhase':
				y0 =[self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],\
					 self.Vdc_t[-1],
					 self.xDC_t[-1],self.xQ_t[-1],
					 self.xPLL_t[-1],self.wte_t[-1]]
		
			elif type(self.PV_model).__name__ == 'SolarPVDERSinglePhaseConstantVdc':
				y0 =[self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],\
					 self.xP_t[-1],self.xQ_t[-1],
					 self.xPLL_t[-1],self.wte_t[-1]]
		
			elif type(self.PV_model).__name__ == 'SolarPVDERThreePhaseConstantVdc':
				y0 =[self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],
					 self.ibR_t[-1], self.ibI_t[-1], self.xbR_t[-1], self.xbI_t[-1], self.ubR_t[-1],self.ubI_t[-1],
					 self.icR_t[-1], self.icI_t[-1], self.xcR_t[-1], self.xcI_t[-1], self.ucR_t[-1],self.ucI_t[-1],
					 self.xP_t[-1],self.xQ_t[-1],
					 self.xPLL_t[-1],self.wte_t[-1]]		
		
			elif self.DER_model_type == 'SolarPVDERThreePhaseBalanced':
				y0 =[self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],
					 self.Vdc_t[-1],
					 self.xDC_t[-1],self.xQ_t[-1],
					 self.xPLL_t[-1],self.wte_t[-1]]
		
			elif self.DER_model_type  == 'SolarPVDERThreePhaseNumba':
				y0 = [self.iaR_t[-1], self.iaI_t[-1], self.xaR_t[-1], self.xaI_t[-1], self.uaR_t[-1],self.uaI_t[-1],\
					  self.ibR_t[-1], self.ibI_t[-1], self.xbR_t[-1], self.xbI_t[-1], self.ubR_t[-1],self.ubI_t[-1],\
					  self.icR_t[-1], self.icI_t[-1], self.xcR_t[-1], self.xcI_t[-1], self.ucR_t[-1],self.ucI_t[-1],\
					  self.Vdc_t[-1],
					  self.xDC_t[-1],self.xQ_t[-1],
					  self.xPLL_t[-1],self.wte_t[-1]]
		
			return y0
		except:
			LogUtil.exception_handler()


	#@property
	def t_calc(self):
		"""Vector of time steps for simulation"""
		try:
			#if (self.tStop - self.tStart) <= self.tInc:
			#	self.tStop =  self.tStart  + 1e-6 #+ self.tInc
			return np.arange(self.tStart, self.tStop + self.tInc, self.tInc)
		except:
			LogUtil.exception_handler()


	def check_jac_availability(self):
		"""Check if Jacobian matrix is available."""		
		try:
			if self.jacFlag:
				if not self.DER_model_type in self.jac_list:
					raise ValueError('{}:Jacobian matrix is not available for DER model:{}'.format(self.name,self.DER_model_type))
		except:
			LogUtil.exception_handler()


	def initialize_y0_t(self):
		"""Initialize y0_t."""
		try:
			self.iaR_t = np.array([self.PV_model.y0[0]])
			self.iaI_t = np.array([self.PV_model.y0[1]])
			self.xaR_t = np.array([self.PV_model.y0[2]])
			self.xaI_t = np.array([self.PV_model.y0[3]])
			self.uaR_t = np.array([self.PV_model.y0[4]])
			self.uaI_t = np.array([self.PV_model.y0[5]])
		
			if type(self.PV_model).__name__ == 'SolarPVDERSinglePhase':			
				self.Vdc_t = np.array([self.PV_model.y0[6]]) #DC link voltage variable
			
				self.xDC_t = np.array([self.PV_model.y0[7]]) #DC link voltage control variable
				self.xQ_t = np.array([self.PV_model.y0[8]]) #Reactive power control variable
			
				self.xPLL_t = np.array([self.PV_model.y0[9]]) #PLL variables			
				self.wte_t = np.array([self.PV_model.y0[10]]) #Frequency integration to get angle
		
			elif (self.DER_model_type == 'SolarPVDERThreePhase') or (self.DER_model_type == 'SolarPVDERThreePhaseNumba'):
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
			
			elif type(self.PV_model).__name__ == 'SolarPVDERThreePhaseConstantVdc':
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
			
				self.Vdc_t = np.array([self.PV_model.Vdc]) #Voltage is constant
			
				self.xP_t = np.array([self.PV_model.y0[18]])  #Active power control variable
				self.xQ_t = np.array([self.PV_model.y0[19]])  #Reactive power control variable
			
				self.xPLL_t = np.array([self.PV_model.y0[20]]) #PLL variables			
				self.wte_t = np.array([self.PV_model.y0[21]]) #Frequency integration to get angle
		
			elif type(self.PV_model).__name__ == 'SolarPVDERThreePhaseBalanced':			
			
				ia_t = self.iaR_t+self.iaI_t*1j
				xa_t = self.xaR_t+self.xaI_t*1j
				ua_t = self.uaR_t+self.uaI_t*1j
			
				ib_t = utility_functions.Ub_calc(ia_t)
				xb_t = utility_functions.Ub_calc(xa_t)
				ub_t = utility_functions.Ub_calc(ua_t)
			
				ic_t = utility_functions.Uc_calc(ia_t)
				xc_t = utility_functions.Uc_calc(xa_t)
				uc_t = utility_functions.Uc_calc(ua_t)			
			
				self.ibR_t = ib_t.real 
				self.ibI_t = ib_t.imag 
				self.xbR_t = xb_t.real 
				self.xbI_t = xb_t.imag 
				self.ubR_t = ub_t.real
				self.ubI_t = ub_t.imag 

				self.icR_t = ic_t.real 
				self.icI_t = ic_t.imag 
				self.xcR_t = xc_t.real 
				self.xcI_t = xc_t.imag 
				self.ucR_t = uc_t.real 
				self.ucI_t = uc_t.imag
			
				self.Vdc_t = np.array([self.PV_model.y0[6]]) #DC link voltage variable
			
				self.xDC_t = np.array([self.PV_model.y0[7]]) #DC link voltage control variable
				self.xQ_t = np.array([self.PV_model.y0[8]]) #Reactive power control variable
			
				self.xPLL_t = np.array([self.PV_model.y0[9]]) #PLL variables			
				self.wte_t = np.array([self.PV_model.y0[10]]) #Frequency integration to get angle
		
			elif type(self.PV_model).__name__ == 'SolarPVDERSinglePhaseConstantVdc':		
				self.Vdc_t = np.array([self.PV_model.Vdc]) #Voltage is constant
			
				self.xP_t = np.array([self.PV_model.y0[6]])  #Active power control variable
				self.xQ_t = np.array([self.PV_model.y0[7]])  #Reactive power control variable
			
				self.xPLL_t = np.array([self.PV_model.y0[8]]) #PLL variables			
				self.wte_t = np.array([self.PV_model.y0[9]]) #Frequency integration to get angle
		except:
			LogUtil.exception_handler()


	def reset_stored_trajectories(self):
		"""Reset for plotting."""
		try:
			self._t_t = np.array(0.0)
			self.Vdc_t = self._Vdc_t = np.array(self.PV_model.Vdc)
		
			self.ia_t = self._ia_t = np.array(self.PV_model.ia)
			self.ma_t = self._ma_t = np.array(self.PV_model.ma)
			self.vta_t = self._vta_t = np.array(self.PV_model.vta)
			self.va_t = self._va_t = np.array(self.PV_model.va)
		
			self.ma_absolute_t = self._ma_absolute_t = np.array(abs(self.PV_model.ma))
			self.Varms_t  = self._Varms_t = np.array(abs(self.PV_model.va)/math.sqrt(2))
		
			if type(self.PV_model).__name__ in templates.three_phase_models:
				self.mb_absolute_t = self._mb_absolute_t = np.array(abs(self.PV_model.mb))
				self.mc_absolute_t = self._mc_absolute_t = np.array(abs(self.PV_model.mc))
			
				self.Vbrms_t  = self._Vbrms_t = np.array(abs(self.PV_model.vb)/math.sqrt(2))
				self.Vcrms_t  = self._Vcrms_t = np.array(abs(self.PV_model.vc)/math.sqrt(2))
			
				self.ib_t = self._ib_t = np.array(self.PV_model.ib)
				self.mb_t = self._mb_t = np.array(self.PV_model.mb)
				self.vtb_t = self._vtb_t = np.array(self.PV_model.vtb)
				self.vb_t = self._vb_t = np.array(self.PV_model.vb)
			
				self.ic_t = self._ic_t = np.array(self.PV_model.ic)
				self.mc_t = self._mc_t = np.array(self.PV_model.mc)
				self.vtc_t = self._vtc_t = np.array(self.PV_model.vtc)
				self.vc_t = self._vc_t = np.array(self.PV_model.vc)
		
			self.Irms_t = self._Irms_t = np.array(self.PV_model.Irms)
			self.Ppv_t = self._Ppv_t = np.array(self.PV_model.Ppv)
			self.S_PCC_t = self._S_PCC_t = np.array(self.PV_model.S_PCC)
			self.S_t = self._S_t = np.array(self.PV_model.S)
			self.Vtrms_t = self._Vtrms_t = np.array(self.PV_model.Vtrms)
			self.Vrms_t = self._Vrms_t = np.array(self.PV_model.Vrms)
		except:
			LogUtil.exception_handler()


	def ODE_model(self,y,t):
		""" Combine all derivatives."""
		try:
			y1 = y[0:self.PV_model.n_ODE]
			if self.PV_model.standAlone:
				self.grid_model.steady_state_model(t)

			y = self.PV_model.ODE_model(y1,t)
			if self.DEBUG_SIMULATION:
				self.debug_simulation(t)

			return y
		except:
			LogUtil.exception_handler()


	def jac_ODE_model(self,y,t):
		""" Combine all derivatives."""
		try:
			y1 = y[0:self.PV_model.n_ODE]
			if self.PV_model.standAlone:
				self.grid_model.steady_state_model(t)
			y = self.PV_model.jac_ODE_model(y1,t)

			return y
		except:
			LogUtil.exception_handler()


	def debug_simulation(self,t):
		""" Print to terminal for debugging."""
		try:
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
		except:
			LogUtil.exception_handler()


	def time_series_PCC_HV_side_voltage(self):
		"""Calculate time series PCC voltage."""
		try:
			assert len(self.ia_t) == len(self.vag_t) != None, "States must be available from simulation."
			self.vaHV_t = self.vag_t + (self.ia_t/self.PV_model.a)*self.grid_model.Z2 - (self.va_t/self.PV_model.a)*(self.grid_model.Z2/self.Zload1_t)

			if type(self.PV_model).__name__ in templates.single_phase_models:	
				self.vbHV_t = self.vbg_t 
				self.vcHV_t = self.vcg_t 
			elif type(self.PV_model).__name__ in templates.three_phase_models:
				self.vbHV_t = self.vbg_t + (self.ib_t/self.PV_model.a)*self.grid_model.Z2 - (self.vb_t/self.PV_model.a)*(self.grid_model.Z2/self.Zload1_t)
				self.vcHV_t = self.vcg_t + (self.ic_t/self.PV_model.a)*self.grid_model.Z2 - (self.vc_t/self.PV_model.a)*(self.grid_model.Z2/self.Zload1_t)
		except:
			LogUtil.exception_handler()


	def time_series_duty_cycle(self):
		"""Calculate time series PCC voltage."""
		try:
			self.ma_t =  utility_functions.m_time_series(self.ua_t,self.xa_t,self.PV_model.Kp_GCC)
			self.maR_t = self.ma_t.real 
			self.maI_t = self.ma_t.imag
		
			self.ma_absolute_t = utility_functions.Uabsolute_time_series(self.ma_t)
				
			if type(self.PV_model).__name__ in templates.three_phase_models:
				self.mb_t =  utility_functions.m_time_series(self.ub_t,self.xb_t,self.PV_model.Kp_GCC)
				self.mc_t =  utility_functions.m_time_series(self.uc_t,self.xc_t,self.PV_model.Kp_GCC)

				self.mbR_t = self.mb_t.real 
				self.mbI_t = self.mb_t.imag
				self.mcR_t = self.mc_t.real 
				self.mcI_t = self.mc_t.imag
			
				self.mb_absolute_t = utility_functions.Uabsolute_time_series(self.mb_t)
				self.mc_absolute_t = utility_functions.Uabsolute_time_series(self.mc_t)
		except:
			LogUtil.exception_handler()


	def time_series_PCC_LV_side_voltage(self):
		"""Calculate time series PCC voltage."""
		try:
			if self.PV_model.standAlone:
				self.va_t = ((self.vag_t+(self.ia_t/self.PV_model.a)*self.grid_model.Z2)/(self.PV_model.a) +self.ia_t*self.PV_model.Z1)*((self.Zload1_t*self.PV_model.a*self.PV_model.a)/((self.PV_model.a*self.PV_model.a*(self.PV_model.Z1+self.Zload1_t))+self.grid_model.Z2))
			
			else:
				self.va_t = np.repeat(self.PV_model.gridVoltagePhaseA,len(self.t))
			
			self.vaR_t = self.va_t.real
			self.vaI_t = self.va_t.imag
		
			if type(self.PV_model).__name__ in templates.three_phase_models:
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
		except:
			LogUtil.exception_handler()


	def time_series_inv_terminal_voltage(self):
		"""Calculate time series inverter terminal voltage."""
		try:
			self.vta_t =  utility_functions.Vinv_terminal_time_series(self.ma_t,self.Vdc_t)
			if type(self.PV_model).__name__ in templates.three_phase_models:
				self.vtb_t =  utility_functions.Vinv_terminal_time_series(self.mb_t,self.Vdc_t)
				self.vtc_t =  utility_functions.Vinv_terminal_time_series(self.mc_t,self.Vdc_t)
		except:
			LogUtil.exception_handler()


	def time_series_RMS(self):
		"""Calculate time series RMS quantities."""
		try:
			if type(self.PV_model).__name__ in templates.single_phase_models:
				self.Vtrms_t = utility_functions.Urms_time_series(self.vta_t,self.vta_t,self.vta_t)
				self.Vrms_t = utility_functions.Urms_time_series(self.va_t,self.va_t,self.va_t)
				self.Irms_t = utility_functions.Urms_time_series(self.ia_t,self.ia_t,self.ia_t)
				self.Varms_t = self.Vrms_t
			elif type(self.PV_model).__name__ in templates.three_phase_models:
				self.Vtrms_t = utility_functions.Urms_time_series(self.vta_t,self.vtb_t,self.vtc_t)
				self.Vrms_t = utility_functions.Urms_time_series(self.va_t,self.vb_t,self.vc_t)
				self.Irms_t = utility_functions.Urms_time_series(self.ia_t,self.ib_t,self.ic_t)	   
				self.Varms_t = utility_functions.Uphrms_time_series(self.va_t)
				self.Vbrms_t = utility_functions.Uphrms_time_series(self.vb_t)
				self.Vcrms_t = utility_functions.Uphrms_time_series(self.vc_t)

			if self.PV_model.standAlone:
				self.Vgrms_t = utility_functions.Urms_time_series(self.vag_t,self.vbg_t,self.vcg_t)
				self.Vhvrms_t= utility_functions.Urms_time_series(self.vaHV_t,self.vbHV_t,self.vcHV_t)
		except:
			LogUtil.exception_handler()


	def time_series_standalone_grid(self):
		"""Time series grid voltage and frequency for standalone model."""
		try:
			self.vag_t = []
			self.vbg_t = []
			self.vcg_t = []
			self.wgrid_t = []
			for i,t in enumerate(self.t):   #Loop through grid events and calculate wgrid at each time step

					Vagrid_new,self.grid_model.wgrid = self.simulation_events.grid_events(t)
				
					#Conversion of grid voltage setpoint
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
		
			if not self.LOOP_MODE:
				self.simulation_events.reset_event_counters() #reset event counters
		except:
			LogUtil.exception_handler()


	def time_series_PLL(self):
		"""Calculate time series PLL and d-q quantities."""
		try:
			self.vd_t = []
			self.vq_t = []
			self.we_t = []
		
			for i,t in enumerate(self.t):	 #Loop through time steps and calculate d-q values at each time step
			
				self.PV_model.wte = self.wte_t[i]
				self.PV_model.xPLL = self.xPLL_t[i]
			
				if type(self.PV_model).__name__ in templates.single_phase_models:
					self.PV_model.valpha = utility_functions.phasor_to_time_1phase(self.va_t[i],self.wgrid_t[i],t)
					self.PV_model.vbeta =utility_functions.phasor_to_time_1phase(self.va_t[i]*pow(math.e,-1j*(math.pi/2)),self.wgrid_t[i],t)
					self.PV_model.vd,self.PV_model.vq = utility_functions.alpha_beta_to_d_q(self.PV_model.valpha,self.PV_model.vbeta,self.PV_model.wte)
			
				elif type(self.PV_model).__name__ in templates.three_phase_models:
					self.PV_model.vat = utility_functions.phasor_to_time_1phase(self.va_t[i],self.wgrid_t[i],t)
					self.PV_model.vbt = utility_functions.phasor_to_time_1phase(self.vb_t[i],self.wgrid_t[i],t)
					self.PV_model.vct = utility_functions.phasor_to_time_1phase(self.vc_t[i],self.wgrid_t[i],t)
				
					self.PV_model.vd,self.PV_model.vq,self.PV_model.v0 = utility_functions.abc_to_dq0(self.PV_model.vat,self.PV_model.vbt,self.PV_model.vct,self.PV_model.wte)
				
				self.PV_model.we = self.PV_model.we_calc() #Calculate inverter frequency from PLL equation
			
				self.vd_t.append(self.PV_model.vd)
				self.vq_t.append(self.PV_model.vq)
				self.we_t.append(self.PV_model.we)
			
			self.vd_t = np.asarray(self.vd_t)
			self.vq_t = np.asarray(self.vq_t)
			self.we_t = np.asarray(self.we_t)
		except:
			LogUtil.exception_handler()


	def time_series_Zload1(self,tOverride=None):
		"""Calculate time series load impedance."""
		try:
			if tOverride:
				tAll=tOverride
			else:
				tAll=self.t

			self.Zload1_t = []
			for i,t in enumerate(tAll):	   #Loop through load events
				_Zload1_actual =  self.simulation_events.load_events(t)  
				_Zload1 = _Zload1_actual/self.PV_model.Zbase
				self.Zload1_t.append(_Zload1)

			self.Zload1_t = np.asarray(self.Zload1_t)
			if not self.LOOP_MODE:
				self.simulation_events.reset_event_counters() #reset event counters
		except:
			LogUtil.exception_handler()


	def time_series_Ppv(self):
		"""Calculate time series Solar PV power output."""
		try:
			self.Ppv_t = []
			self.Sinsol_t = []
			for i,t in enumerate(self.t):	   #Loop through solar events and calculate Ppv at each time step
				self.PV_model.Sinsol,self.PV_model.Tactual =  self.simulation_events.solar_events(t)  #Parse through solar events
				self.PV_model.Vdc = self.Vdc_t[i]
				self.PV_model.Iph = self.PV_model.Iph_calc()
				self.Sinsol_t.append(self.PV_model.Sinsol)
				self.Ppv_t.append(self.PV_model.Ppv_calc(self.PV_model.Vdc_actual))
			self.Sinsol_t = np.asarray(self.Sinsol_t)
			self.Ppv_t = np.asarray(self.Ppv_t)
			if not self.LOOP_MODE:
				self.simulation_events.reset_event_counters() #reset event counters
		except:
			LogUtil.exception_handler()


	def time_series_S(self):
		"""Calculate time series apparent power."""
		try:
			if type(self.PV_model).__name__ in templates.single_phase_models:
				self.S_t = (1/2)*(self.vta_t*self.ia_t.conjugate())
				self.S_PCC_t = (1/2)*(self.va_t*self.ia_t.conjugate())			
			
				if self.PV_model.standAlone:
					self.S_G_t = (1/2)*(-(self.ia_t - (self.va_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vag_t
					self.S_load1_t = (1/2)*(self.va_t*(-(self.va_t/self.Zload1_t)).conjugate())
		
			elif type(self.PV_model).__name__ in templates.three_phase_models:
				self.S_t = (1/2)*(self.vta_t*self.ia_t.conjugate() + self.vtb_t*self.ib_t.conjugate() + self.vtc_t*self.ic_t.conjugate())
				self.S_PCC_t = (1/2)*(self.va_t*self.ia_t.conjugate()+self.vb_t*self.ib_t.conjugate()+self.vc_t*self.ic_t.conjugate())
				 
				if self.PV_model.standAlone:
					self.S_G_t = (1/2)*(-(self.ia_t - (self.va_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vag_t + (-(self.ib_t -(self.vb_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vbg_t + (-(self.ic_t -(self.vc_t/self.Zload1_t))/self.PV_model.a).conjugate()*self.vcg_t
					self.S_load1_t = (1/2)*(self.va_t*(-(self.va_t/self.Zload1_t)).conjugate() + self.vb_t*(-(self.vb_t/self.Zload1_t)).conjugate() + self.vc_t*(-(self.vc_t/self.Zload1_t)).conjugate())
		except:
			LogUtil.exception_handler()


	def time_series_phase_angle(self):
		"""Calculate time series phase angle between grid and inverter voltage."""
		try:
			assert len(self.va_t) == len(self.vta_t) != None, "States must be available from simulation."
			self.phi_at_t = np.angle(self.vta_t)
			self.phi_a_t = np.angle(self.va_t)

			if self.PV_model.standAlone:
				self.phi_aHV_t = np.angle(self.vaHV_t) 
				self.phi_ag_t = np.angle(self.vag_t)		

			if type(self.PV_model).__name__ in templates.three_phase_models:
				self.phi_bt_t = np.angle(self.vtb_t)
				self.phi_ct_t = np.angle(self.vtc_t)
				self.phi_b_t = np.angle(self.vb_t)
				self.phi_c_t = np.angle(self.vc_t)
				if self.PV_model.standAlone:
					self.phi_bHV_t = np.angle(self.vbHV_t)
					self.phi_cHV_t = np.angle(self.vcHV_t)
					self.phi_bg_t = np.angle(self.vbg_t)
					self.phi_cg_t = np.angle(self.vcg_t)
		except:
			LogUtil.exception_handler()


	def time_series_power_transfer(self):
		"""Calculate time series power transfer between grid and PCC."""
		try:
			self.P_transfer = (self.Vhvrms_t*self.Vgrms_t*np.sin(self.phi_aHV_t-self.phi_ag_t))/np.abs(self.grid_model.Z2)
		except:
			LogUtil.exception_handler()


	def collect_solution(self,solution,t=None):
		"""Collect solution."""
		try:
			if self.LOOP_MODE:
				self.t = t
				self.collect_full_trajectory(solution)
				self.collect_last_states()
			else:
				self.collect_full_trajectory(solution)
			LogUtil.logger.debug('{}:Stored solution for {} time points starting at {:.3f} s and ending at {:.3f} s!'.format(self.name,len(self.t),self.t[0],self.t[-1]))
		except:
			LogUtil.exception_handler()


	def collect_last_states(self):
		"""Collect states at last time step."""
		try:
			self._t_t = np.append(self._t_t,self.t[1:])
			self._Vdc_t = np.append(self._Vdc_t,self.Vdc_t[1:])
			self._ia_t = np.append(self._ia_t,self.ia_t[1:])
			self._ma_t = np.append(self._ma_t,self.ma_t[1:])
			self._vta_t = np.append(self._vta_t,self.vta_t[1:])
			self._va_t = np.append(self._va_t,self.va_t[1:])
		
			self._ma_absolute_t = np.append(self._ma_absolute_t,self.ma_absolute_t[1:])
			self._Varms_t = np.append(self._Varms_t,self.Varms_t[1:])
		
			if type(self.PV_model).__name__ in templates.three_phase_models:
				self._mb_absolute_t = np.append(self._mb_absolute_t,self.mb_absolute_t[1:])
				self._mc_absolute_t = np.append(self._mc_absolute_t,self.mc_absolute_t[1:])
			
				self._Vbrms_t = np.append(self._Vbrms_t,self.Vbrms_t[1:])
				self._Vcrms_t = np.append(self._Vcrms_t,self.Vcrms_t[1:])
			
				self._ib_t = np.append(self._ib_t,self.ib_t[1:])
				self._mb_t = np.append(self._mb_t,self.mb_t[1:])
				self._vtb_t = np.append(self._vtb_t,self.vtb_t[1:])
				self._vb_t = np.append(self._vb_t,self.vb_t[1:])
			
				self._ic_t = np.append(self._ic_t,self.ic_t[1:])
				self._mc_t = np.append(self._mc_t,self.mc_t[1:])
				self._vtc_t = np.append(self._vtc_t,self.vtc_t[1:])
				self._vc_t = np.append(self._vc_t,self.vc_t[1:])
		
			self._Irms_t = np.append(self._Irms_t,self.Irms_t[1:])		

			self._Vtrms_t = np.append(self._Vtrms_t,self.Vtrms_t[1:])
			self._Vrms_t = np.append(self._Vrms_t,self.Vrms_t[1:])		
		
			self._Ppv_t = np.append(self._Ppv_t,self.Ppv_t[1:])
			self._S_t = np.append(self._S_t,self.S_t[1:])
			self._S_PCC_t = np.append(self._S_PCC_t,self.S_PCC_t[1:])
		except:
			LogUtil.exception_handler()


	def collect_full_trajectory(self,solution):
		"""Collect full solution from solver."""
		try:
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
				
			if type(self.PV_model).__name__ in templates.three_phase_models:
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
				self.time_series_standalone_grid() #Grid voltage
				self.time_series_Zload1()		
		
			self.time_series_inv_terminal_voltage()
			self.time_series_PCC_LV_side_voltage()
		
			if self.PV_model.standAlone:
				self.time_series_PCC_HV_side_voltage()
		
			self.time_series_S()
			self.time_series_RMS()
			self.time_series_Ppv()
		
			self.time_series_phase_angle()
			#self.time_series_power_transfer()
			#self.time_series_wgrid()
			if self.PV_model.standAlone:
				self.time_series_PLL()
		
			LogUtil.logger.debug("All states collected!")
		except:
			LogUtil.exception_handler()


	def collect_states(self,solution):
		"""Collect states from ode solution."""
		try:
			 #Phase a states
			self.iaR_t = solution[:,0]
			self.iaI_t = solution[:,1]
			self.xaR_t = solution[:,2]
			self.xaI_t = solution[:,3]
			self.uaR_t = solution[:,4]
			self.uaI_t = solution[:,5]
		
			if self.DER_model_type == 'SolarPVDERSinglePhase':
				#DC link voltage variables
				self.Vdc_t = solution[:,6]
				self.xDC_t = solution[:,7]
				self.xQ_t = solution[:,8]
				#PLL variables
				self.xPLL_t = solution[:,9]
				#Frequency integration to get angle
				self.wte_t = solution[:,10]

			elif self.DER_model_type == 'SolarPVDERThreePhase' or self.DER_model_type == 'SolarPVDERThreePhaseNumba':
			
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
		
			elif self.DER_model_type  == 'SolarPVDERSinglePhaseConstantVdc':
				self.Vdc_t = np.array([self.PV_model.Vdc]*len(self.iaR_t))
						
				self.xP_t = solution[:,6] #Active power control variable
				self.xQ_t = solution[:,7] #Reactive power control variable
			
				self.xPLL_t = solution[:,8] #PLL variables			
				self.wte_t = solution[:,9] #Frequency integration to get angle
		
			elif self.DER_model_type == 'SolarPVDERThreePhaseConstantVdc':
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
			
				self.Vdc_t = np.array([self.PV_model.Vdc]*len(self.iaR_t))
			
				self.xP_t = solution[:,18]  #Active power control variable
				self.xQ_t = solution[:,19]  #Reactive power control variable
			
				self.xPLL_t = solution[:,20] #PLL variables			
				self.wte_t = solution[:,21] #Frequency integration to get angle
			
			elif self.DER_model_type == 'SolarPVDERThreePhaseBalanced':
			
				ia_t = self.iaR_t+self.iaI_t*1j
				xa_t = self.xaR_t+self.xaI_t*1j
				ua_t = self.uaR_t+self.uaI_t*1j
			
				ib_t = utility_functions.Ub_calc(ia_t)
				xb_t = utility_functions.Ub_calc(xa_t)
				ub_t = utility_functions.Ub_calc(ua_t)
			
				ic_t = utility_functions.Uc_calc(ia_t)
				xc_t = utility_functions.Uc_calc(xa_t)
				uc_t = utility_functions.Uc_calc(ua_t)			
			
				self.ibR_t = ib_t.real 
				self.ibI_t = ib_t.imag 
				self.xbR_t = xb_t.real 
				self.xbI_t = xb_t.imag 
				self.ubR_t = ub_t.real
				self.ubI_t = ub_t.imag 

				self.icR_t = ic_t.real 
				self.icI_t = ic_t.imag 
				self.xcR_t = xc_t.real 
				self.xcI_t = xc_t.imag 
				self.ucR_t = uc_t.real 
				self.ucI_t = uc_t.imag
			
				#DC link voltage variables
				self.Vdc_t = solution[:,6]
				self.xDC_t = solution[:,7]
				self.xQ_t = solution[:,8]
				#PLL variables
				self.xPLL_t = solution[:,9]
				#Frequency integration to get angle
				self.wte_t = solution[:,10]
		except:
			LogUtil.exception_handler()


	def invert_arrays(self):
		"""Inverter arrays before usage by plots."""
		try:
			self.t_t = self._t_t
			self.Vdc_t = self._Vdc_t
		
			self.ia_t = self._ia_t
			self.ma_t = self._ma_t
			self.vta_t = self._vta_t
			self.va_t = self._va_t
		
			self.ma_absolute_t = self._ma_absolute_t
			self.Varms_t = self._Varms_t
		
			if type(self.PV_model).__name__ in templates.three_phase_models:
				self.mb_absolute_t = self._mb_absolute_t
				self.mc_absolute_t = self._mc_absolute_t
			
				self.Vbrms_t = self._Vbrms_t
				self.Vcrms_t = self._Vcrms_t
			
				self.ib_t = self._ib_t
				self.mb_t = self._mb_t
				self.vtb_t = self._vtb_t
				self.vb_t = self._vb_t
			
				self.ic_t = self._ic_t
				self.mc_t = self._mc_t
				self.vtc_t = self._vtc_t
				self.vc_t = self._vc_t
				
			self.Vtrms_t = self._Vtrms_t
			self.Vrms_t = self._Vrms_t
				
			self.Irms_t = self._Irms_t
		
			self.Ppv_t  = self._Ppv_t 
			self.S_t =  self._S_t 
			self.S_PCC_t = self._S_PCC_t
		except:
			LogUtil.exception_handler()


	def run_simulation(self,gridVoltagePhaseA=None,gridVoltagePhaseB=None,gridVoltagePhaseC=None,y0=None,t=None):
		"""Call the ODE solver and collect states."""
		try:
			self.solution_time = None #Always reset simulation time to None
			if self.LOOP_MODE:
			
				if isinstance(gridVoltagePhaseA,complex) and isinstance(y0,list) and isinstance(t,list):
					self.update_grid_measurements(gridVoltagePhaseA=gridVoltagePhaseA, gridVoltagePhaseB=gridVoltagePhaseB, gridVoltagePhaseC=gridVoltagePhaseC)
				else:
					raise ValueError('Grid voltage (complex scalar), initial states (y0), and time steps (t) should be provided in loop mode.')
			
				if t[0] == 0.0:
					self.t = t
					six.print_("{}:Simulation started in loop mode with a step size of {:.4f} s!".format(self.name,self.t[-1]-self.t[0]))
					if self.jacFlag:
						LogUtil.logger.debug("{}:Analytical Jacobian will be provided to ODE solver.".format(self.name))
				solution,_,_  = self.call_ODE_solver(self.ODE_model,self.jac_ODE_model,y0,t)
			
			else:
				self.t = self.t_calc()
				self.initialize_y0_t()
			
				timer_start = time.time()
				six.print_("{}:Simulation started at {} s and will end at {} s".format(self.name,self.tStart,self.tStop))
				LogUtil.logger.debug('{}:{} solver will be used.'.format(self.name,self.solver_type))
				if self.jacFlag:
					LogUtil.logger.debug("{}:Analytical Jacobian will be provided to ODE solver.".format(self.name))
				
				if self.solver_type != 'odeint':
					self.ode_solver.set_initial_value(self.y0,self.tStart) 
					LogUtil.logger.debug("{}:Resetting {} internal time step to {} and states to:\n{}.".format(self.name,self.solver_type,self.tStart,self.y0))
		
				solution,_,_  = self.call_ODE_solver(self.ODE_model,self.jac_ODE_model,self.y0,self.t)
			
				self.solution_time = time.time() - timer_start
			
				self.show_simulation_time()
				
				self.simulation_events.reset_event_counters() #Reset grid and solar event counters
				self.PV_model.reset_reference_counters() #Reset counters for reference change events
		
			if self.COLLECT_SOLUTION:
				self.collect_solution(solution,t=t)  #Collect full solution
			else:
				self.collect_states(solution)  #Atleast states must be collected	
		except:
			LogUtil.exception_handler()


	def update_grid_measurements(self,gridVoltagePhaseA, gridVoltagePhaseB, gridVoltagePhaseC):
		"""Update grid voltage and frequency in non-standalone model.
		Args:
			 gridVoltagePhaseA (complex): Value of gridVoltagePhaseA
			 gridVoltagePhaseB (complex): Value of gridVoltagePhaseB
			 gridVoltagePhaseC (complex): Value of gridVoltagePhaseC			 
		"""
		try:
			self.PV_model.gridVoltagePhaseA = gridVoltagePhaseA
			if type(self.PV_model).__name__ == 'SolarPVDERThreePhase':
				self.PV_model.gridVoltagePhaseB = gridVoltagePhaseB
				self.PV_model.gridVoltagePhaseC = gridVoltagePhaseC
			elif type(self.PV_model).__name__ == 'SolarPVDERThreePhaseBalanced':
				self.PV_model.gridVoltagePhaseB = utility_functions.Ub_calc(gridVoltagePhaseA)
				self.PV_model.gridVoltagePhaseC = utility_functions.Uc_calc(gridVoltagePhaseA)
		except:
			LogUtil.exception_handler()


	def show_simulation_time(self):
		"""Show simulation time."""
		try:
			print('{}:Simulation was completed in {}'.format(self.name,time.strftime("%H:%M:%S", time.gmtime(self.solution_time))))
		except:
			LogUtil.exception_handler()


	def get_trajectories(self):
		"""Return trajectories as a dictionary."""
		try:
			if self.LOOP_MODE:
				self.invert_arrays()
			trajectory_dictionary = {'ia_t':self.ia_t,'ma_t':self.ma_t,
									 'Vdc_t':self.Vdc_t,'vta_t':self.vta_t,'va_t':self.va_t,
									 'Irms_t':self.Irms_t,'Vtrms_t':self.Vtrms_t,'Vrms_t':self.Vrms_t,
									 'Ppv_t':self.Ppv_t,'S_t':self.S_t,'S_PCC_t':self.S_PCC_t
									}
			if type(self.PV_model).__name__ in templates.three_phase_models:
				trajectory_dictionary.update({'ib_t':self.ib_t,'ic_t':self.ic_t,'mb_t':self.mb_t,'mc_t':self.mc_t,
											  'vtb_t':self.vtb_t,'vtc_t':self.vtc_t,
											  'vb_t':self.vb_t,'vc_t':self.vc_t})	 

			for trajectory in trajectory_dictionary.values():
				assert len(self.t_t) == len(trajectory), 'The length of each trajectory should be equal to the time_steps.'

			return trajectory_dictionary
		except:
			LogUtil.exception_handler()


