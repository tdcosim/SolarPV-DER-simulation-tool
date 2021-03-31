""" Manage simulation results and ODE solver."""

from __future__ import division
import sys
import six
import time

import numpy as np
import math
import cmath
from scipy.integrate import odeint,ode

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from pvder.utility_classes import Utilities
from pvder import defaults
from pvder.logutil import LogUtil


class SimulationResults(Utilities):
	""" Utility class for simulation results."""
	
	count = 0
	SAVE_PLOT_JPEG = False
	SAVE_PLOT_SVG = False
	figure_DPI = defaults.FIGURE_DPI
	
	parameters= {"figure":{"height": defaults.FIGURE_HEIGHT,
						 "width": defaults.FIGURE_WIDTH
						}
				}
	
	available_plot_types = ['power','active_power','active_power_Ppv_Pac_PCC','active_power_Pac_PCC','reactive_power','reactive_power_Q_PCC','reactive_power_Q_PCC_smoothed','voltage','voltage_Vdc','voltage_HV','voltage_HV_imbalance','voltage_LV','voltage_Vpcclv','voltage_Vpcclv_smoothed','current','phase_angle','frequency','duty_cycle','voltage_Vpcclv_all_phases']
	
	def __init__(self,simulation,figure_index=1,PER_UNIT=True,font_size=18,PLOT_TITLE=True,SOLUTION_TIME=True,verbosity='INFO',identifier=None):
		"""Creates an instance of `SimulationResults`.
		Args:
		simulation: An instance of `GridSimulation`.
		figure_index: An integer specifying the figure index.
		font_size: An integer spcifying the font size to be used withing plots.
		PLOT_TITLE: A boolean specifying whether the title will be displayed in plots.
		verbosity: A string specifying the verbosity level (DEBUG,INFO,WARNING,ERROR).
		"""
		try:
			# do a lazy import. This is against PEP style guidelines. However,
			# matplotlib has a ~25mb memory foot print. Hence, if each opendssapi.worker
			# imports this then memory consumption becomes a burden very quickly. This
			# way only when simulationresults class is called, matplotlib will be
			# imported.
			#import matplotlib.pyplot as plt
			#import matplotlib.ticker as ticker

			#Increment count to keep track of number of simulation results instances
			SimulationResults.count = SimulationResults.count + 1
		
			#Generate a name for the instance
			self.name_instance(identifier)
			self.figure_index =1
			self.simulation = simulation
			self.PER_UNIT = PER_UNIT
			self.PLOT_TITLE = PLOT_TITLE
			self.SOLUTION_TIME = SOLUTION_TIME
			self.font_size = font_size
		except:
			LogUtil.exception_handler()


	def change_units(self):
		"""Change units from per unit to S.I. or vice versa."""
		try:
			if self.PER_UNIT:
				self.V_multiplier = 1.0
				self.I_multiplier = 1.0
				self.S_multiplier = 1.0
				LogUtil.logger.debug('Using per unit quantities!')
			else:
				self.V_multiplier = self.simulation.Vbase
				self.I_multiplier = self.simulation.Ibase
				self.S_multiplier = self.simulation.Sbase
				LogUtil.logger.debug('Using SI quantities!')
		except:
			LogUtil.exception_handler()


	def group_quantities_for_plotting(self,plot_type):
		"""Plot power from simulation.
		Args:
		plot_type (string): The quantities to be plotted.
		
		Raises:
			 ValueError: If `plot_type` is not available. 
		"""
		try:
			if plot_type not in SimulationResults.available_plot_types:
	 #{'power','active_power','active_power_Ppv_Pac_PCC','active_power_Pac_PCC','reactive_power','reactive_power_Q_PCC','reactive_power_Q_PCC_smoothed','voltage','voltage_Vdc','voltage_HV','voltage_HV_imbalance','voltage_LV','voltage_Vpcclv','voltage_Vpcclv_smoothed','current','phase_angle','frequency','duty_cycle','voltage_Vpcclv_all_phases'}:
				raise ValueError('Unknown plot type: ' + str(plot_type))
		
			if self.simulation.LOOP_MODE:
				self.simulation.invert_arrays()
			time = self.simulation.t_t	
			self.change_units()
		
			if self.PER_UNIT:
			 _voltage_label = 'V (p.u.)'
			 _current_label = 'A (p.u.)'
			 _power_label = 'W/VAR (p.u.)'
			 _active_power_label = 'W (p.u.)'
			 _reactive_power_label = 'VAR (p.u.)'
			else:
			 _voltage_label = 'V'
			 _current_label = 'A'
			 _power_label = 'W/VAR'
			 _active_power_label = 'Active Power (W)'
			 _reactive_power_label = 'Reactive Power (VAR)'
			 
			if plot_type == 'voltage':
				plot_values = [self.simulation.Vdc_t*self.V_multiplier,self.simulation.Vtrms_t*self.V_multiplier,self.simulation.Vrms_t*self.V_multiplier,self.simulation.vdg_t*self.V_multiplier]
				if self.simulation.PV_model.standAlone:
				 plot_values = plot_values + [self.simulation.Vhvrms_t*self.V_multiplier,self.simulation.Vgrms_t*self.V_multiplier]
			
				#legends=['DC link voltage','Inverter terminal voltage','PCC LV side voltage','PCC HV side voltage','Grid voltage','Grid voltage - d component']
				legends=[r"$V_{DC}$",r"$v^{inv}_{rms}$",r"$v^{PCC-LV}_{rms}$",r"$v^{PCC-HV}_{rms}$",r"$v^{grid}_{rms}$",r"$v^{grid}_{d}$"]
			
				plot_title='All available voltage quantities!'
				y_labels=_voltage_label
		
			elif plot_type == 'voltage_Vdc':
				plot_values = [self.simulation.Vdc_t*self.V_multiplier]
				legends=[r"$V_{DC}$"]
				#legends=['Vdc']
				plot_title='DC link capacitor: Voltage'
				y_labels=_voltage_label
		
			elif plot_type == 'voltage_LV':
				plot_values = [self.simulation.Vtrms_t*self.V_multiplier,self.simulation.Vrms_t*self.V_multiplier]
				#legends=['Vt','Vpcclv']
			
				legends=[r"$v^{inv}_{rms}$",r"$v^{PCC-LV}_{rms}$"] 
				plot_title='LV side: L-G RMS voltage'
				y_labels=_voltage_label
		
			elif plot_type == 'voltage_Vpcclv':
				plot_values = [self.simulation.Vrms_t*self.V_multiplier]
				#legends=['Vpcclv']
			
				legends=[r"$v^{PCC-LV}_{rms}$"] 
				plot_title='PCC-LV side: L-G RMS voltage'
				y_labels=_voltage_label
		
			elif plot_type == 'voltage_Vpcclv_all_phases':
				plot_values = [self.simulation.Varms_t*self.V_multiplier]
				legends=[r"$v^{PCC-LV}_{arms}$"] 
				#legends=['Vpcclv-a']
			
				if type(self.simulation.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
					plot_values = plot_values + [self.simulation.Vbrms_t*self.V_multiplier,self.simulation.Vcrms_t*self.V_multiplier]
					legends = legends + [r"$v^{PCC-LV}_{brms}$",r"$v^{PCC-LV}_{crms}$"] 
			
				plot_title='PCC-LV side: L-G RMS voltage - all phases'
				y_labels=_voltage_label		
		
			elif plot_type == 'voltage_HV':
				if not self.simulation.PV_model.standAlone:
					raise ValueError('Plot '+str(plot_type)+' is only available in Stand Alone mode!')
				
				plot_values = [self.simulation.Vhvrms_t*self.V_multiplier,self.simulation.Vgrms_t*self.V_multiplier]
				#legends=['Vpcchv','Vgrid']
			
				legends=[r"$v^{PCC-HV}_{rms}$",r"$v^{grid}_{rms}$"] 
				plot_title='PCC HV side and Grid voltage source:L-G RMS voltage'
				y_labels=_voltage_label
		
			elif plot_type == 'current':
				plot_values = [self.simulation.Irms_t*self.I_multiplier]
				#legends=['Inverter RMS current']
				legends=[r"$i^{inv}_{rms}$"] 
				plot_title='Inverter terminal: RMS current'
				y_labels = _current_label
		
			elif plot_type == 'power':
				plot_values = [self.simulation.Ppv_t*self.S_multiplier,self.simulation.S_t.real*self.S_multiplier,self.simulation.S_t.imag*self.S_multiplier,
	 self.simulation.S_PCC_t.real*self.S_multiplier,self.simulation.S_PCC_t.imag*self.S_multiplier]
			
				if self.simulation.PV_model.standAlone:
					plot_values = plot_values + [self.simulation.S_load1_t.real*self.S_multiplier,self.simulation.S_G_t.real*self.S_multiplier]
			
				#legends=['Ppv','Pac','Q','Pac_PCC','Q_PCC','Pac_load1','Pac_G']
			
				legends=[r"$P^{PV}_{DC}$",r"$P^{inv}_{AC}$",r"$Q^{inv}$",r"$P^{PCC-LV}_{AC}$",r"$Q^{PCC-LV}$",r"$P^{load1}_{AC}$",r"$P^{G}_{AC}$"] 
				plot_title='All available power quantities!'
				y_labels=_power_label
		
			elif plot_type == 'active_power':
				plot_values = [self.simulation.Ppv_t*self.S_multiplier,self.simulation.S_t.real*self.S_multiplier,self.simulation.S_PCC_t.real*self.S_multiplier]
				if self.simulation.PV_model.standAlone:
					plot_values = plot_values + [self.simulation.S_load1_t.real*self.S_multiplier,self.simulation.S_G_t.real*self.S_multiplier]
				#legends=['Ppv','Pac','Pac_PCC','Pac_load1','Pac_G']
			
				legends=[r"$P^{PV}_{DC}$",r"$P^{inv}_{AC}$",r"$P^{PCC-LV}_{AC}$",r"$P^{load1}_{AC}$",r"$P^{G}_{AC}$"] 
				plot_title='PV,Inverter,PCC-LV,Load,& Grid voltage:Active power'
				y_labels=_active_power_label		
		
			elif plot_type == 'active_power_Ppv_Pac_PCC':
				plot_values = [self.simulation.Ppv_t*self.S_multiplier,self.simulation.S_PCC_t.real*self.S_multiplier]
			
				#legends=['Ppv','Pac_PCC']
			
				legends=[r"$P^{PV}_{DC}$",r"$P^{PCC-LV}_{ac}$"] 
				plot_title='PV,PCC-LV:Active power'
				y_labels=_active_power_label							 
		
			elif plot_type == 'active_power_Pac_PCC':
				plot_values = [self.simulation.S_PCC_t.real*self.S_multiplier]
			
				#legends=['Pac_PCC']
			
				legends=[r"$P^{PCC-LV}_{ac}$"] 
				plot_title='PCC-LV:Active power'
				y_labels=_active_power_label 
		
			elif plot_type == 'reactive_power':
				plot_values = [self.simulation.S_t.imag*self.S_multiplier,self.simulation.S_PCC_t.imag*self.S_multiplier]
				legends=[r"$Q^{inv}$",r"$Q^{PCC-LV}$"]
				plot_title='Inverter,PCC-LV:Reactive power'
				if self.simulation.PV_model.standAlone:
					plot_values = plot_values + [self.simulation.S_load1_t.imag*self.S_multiplier,self.simulation.S_G_t.imag*self.S_multiplier]
					legends= legends + [r"$Q^{load1}$",r"$Q^{G}$"] 
					plot_title = 'Inverter,PCC-LV,Load,& Grid voltage:Reactive power'
			
				#legends=['Q','Q_PCC','Q_load1','Q_G']						
				y_labels=_reactive_power_label
		
			elif plot_type == 'reactive_power_Q_PCC':
				plot_values = [self.simulation.S_PCC_t.imag*self.S_multiplier]
				#legends=['Q_PCC']
			
				legends=[r"$Q^{PCC-LV}$"]
				plot_title='PCC-LV:Reactive power'
				y_labels=_reactive_power_label
	
			elif plot_type == 'phase_angle':
				if not self.simulation.PV_model.standAlone:
					raise ValueError('Plot {} is only available in Stand Alone mode!'.format(plot_type))
			
				plot_values = [self.simulation.phi_at_t,self.simulation.phi_a_t]
				if self.simulation.PV_model.standAlone:
					plot_values = plot_values + [self.simulation.phi_ag_t,self.simulation.phi_a_t-self.simulation.phi_ag_t]
				
				#legends=['theta_vat','theta_va','theta_vag','theta_va-theta_vag']
			
				legends=[r"$\phi^{inv}_{a}$",r"$\phi^{PCC-LV}_{a}$",r"$\phi^{grid}_{a}$",r"$\phi^{PCC-LV}_{a}-\phi^{grid}_{a}$"]
				plot_title='All available phase angle quantities!'
				y_labels='radians'
		
			elif plot_type == 'frequency':
				if not self.simulation.PV_model.standAlone:
					raise ValueError('Plot {} is only available in Stand Alone mode!'.format(plot_type))
				
				plot_values = [self.simulation.wgrid_t,self.simulation.we_t] #legends=['Grid frequency','PLL frequency']
			
				legends=[r"$\omega_{grid}$",r"$\omega_{PLL}$"]
				plot_title='Grid voltage source and PLL: Frequency'
				y_labels='radians/s'
			 
			elif plot_type == 'duty_cycle':
				plot_values = [self.simulation.ma_absolute_t]
				#legends = [r"$\m^{absolute}_{a}$"]
				legends = ["ma"]
				if type(self.simulation.PV_model).__name__ == 'SolarPV_DER_ThreePhase':
					plot_values = plot_values + [self.simulation.mb_absolute_t,self.simulation.mc_absolute_t]
					#legends = legends + [r"$\m^{absolute}_{b}$",r"$\mc^{absolute}_{c}$"]
					legends = legends + ["mb","mc"]
					#legends=['m absolute']
						
				plot_title = 'Inverter: Duty cycle'
				y_labels = 'p.u.'
		
			plot_title = self.simulation.PV_model.name + ' -- ' + plot_title
		
			if self.simulation.solution_time is not None and self.SOLUTION_TIME:
				plot_title = plot_title + '\nSolution time:{:.4f} seconds'.format(self.simulation.solution_time) 
		
			for plot_value in plot_values:
				assert (plot_value.size == time.size) and time.size > 1, 'Number of time points should be greater than one and equal to number of value points!'
		
			return time,plot_values,legends,plot_title,y_labels
		except:
			LogUtil.exception_handler()


	def plot_DER_simulation(self,plot_type = 'power'):
		"""Plot desired electrical quantity from simulation.
		Args:
		plot_type (str): Specify type of plot (to see available plot types use SimulationResults.available_plot_types).
		"""
		try:
			time,plot_values,legends,plot_title,y_labels = self.group_quantities_for_plotting(plot_type)
			time_values = [time]*len(plot_values)
			self.plot_multiple(time_values,plot_values,legends,plot_title,y_labels)
			self.save_plot(plt,plot_type)
			self.figure_index =self.figure_index+1
		except:
			LogUtil.exception_handler()


	def plot_multiple(self,time_values,plot_values,legends,plot_title,y_labels):
		"""Function to plot multiple time series on same plot."""
		try:
			assert len(plot_values) == len(time_values) == len(legends)," The number of legends should be equal to the number of quantities"
			fig = plt.figure(self.figure_index, figsize=(self.parameters["figure"]["width"], self.parameters["figure"]["height"]))
			for i,item in enumerate(plot_values):
				plt.plot(time_values[i],plot_values[i],label=legends[i])
			plt.ylabel(y_labels,weight = "bold", fontsize=self.font_size)
			plt.xlabel('Time (s)',weight = "bold", fontsize=self.font_size)

			if self.PLOT_TITLE:
				plt.title(plot_title,weight = "bold", fontsize=self.font_size) 
			plt.xlim(0, time_values[0][-1])
			plt.legend(fontsize =self.font_size)

			ax = plt.axes()
			ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
			ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
		
			ax.tick_params(axis = 'both', which = 'major',labelsize = self.font_size)
			ax.tick_params(axis = 'both', which = 'minor',labelsize = self.font_size)

			plt.tight_layout()
			plt.show()
		except:
			LogUtil.exception_handler()


	def save_plot(self,plot_object,plot_name='results'):
		"""Save the plots."""
		try:
			if self.SAVE_PLOT_JPEG:
				plot_object.savefig(plot_name+".jpg", dpi=self.figure_DPI)
			if self.SAVE_PLOT_SVG:
				plot_object.savefig(plot_name+".svg", dpi=self.figure_DPI)
		except:
			LogUtil.exception_handler()


	def compare_with_external(self,external_time_values,external_plot_values,external_plot_legends,plot_type='power'):
		"""Compare with external simulation."""
		try:
			assert external_plot_values is not None, 'There should be external values to plot'
			assert len(external_plot_values) == len(external_plot_legends),'Legends should be equal to number of plots'

			time,internal_plot_values,internal_legends,plot_title,y_labels = self.group_quantities_for_plotting(plot_type)
			internal_time_values = [time]*len(internal_plot_values)
			time_values = internal_time_values+external_time_values
			plot_values = internal_plot_values + external_plot_values
			legends = internal_legends + external_plot_legends
			self.plot_multiple(time_values,plot_values,legends,plot_title,y_labels)
			self.figure_index =self.figure_index+1
		except:
			LogUtil.exception_handler()


class SimulationUtilities():
	""" Utility class for dynamic simulations."""

	max_steps = defaults.max_steps #Max steps to be used by solver before producing error
	solver_list = ['odeint','ode-vode-bdf']

	def call_ODE_solver(self,derivatives,jacobian,y,t):
		"""Call the SciPy ODE solver."""
		try:
			if self.solver_type == 'odeint':
				solution,infodict = self.call_odeint_solver(derivatives,jacobian,y,t)
				self.check_simulation(infodict,t) #Check whether solver successful for all time intervals
			elif self.solver_type == 'ode-vode-bdf':
				solution,infodict = self.call_ode_solver()			
			else:
				LogUtil.logger.debug('Solver not found!')			
			return solution,infodict,self.SOLVER_CONVERGENCE
		except:
			LogUtil.exception_handler()


	def call_odeint_solver(self,derivatives,jacobian,y,t):
		"""Use the SciPy odeint solver."""
		try:
			if self.jacFlag:
				solution,infodict = odeint(derivatives,y,t,Dfun=jacobian,full_output=1,printmessg=True,\
										 hmax = 1/120.,mxstep=self.max_steps,atol=1e-4,rtol=1e-4)
			else:
				solution,infodict =odeint(derivatives,y,t,full_output=1,printmessg=True,\
										 hmax = 1/120.,mxstep=self.max_steps,atol=1e-4,rtol=1e-4)
			return solution,infodict
		except:
			LogUtil.exception_handler()


	def call_ode_solver(self):
		"""Use the SciPy ode solver."""
		try:
			solution = self.ode_solver.y
			self.t =np.array([self.ode_solver.t]) #np.array([0.0])
			info = []
			#print('t:',self.ode_solver.t,'y1:',solution[0])
			#while self.ode_solver.successful() and self.ode_solver.t < self.tStop:
			while self.ode_solver.successful() and(self.tStop + 1e-6 - self.ode_solver.t) >= self.tInc:		
				y = self.ode_solver.integrate(self.ode_solver.t+self.tInc)
				return_code = self.ode_solver.get_return_code()
				info.append(return_code)
				#print('t:',self.ode_solver.t,'y1:',y[0])
				solution = np.vstack((solution, y))
				self.t = np.hstack((self.t, np.array([self.ode_solver.t])))
			#print('Solution shape:',solution.shape)
			#print('Time steps shape:',self.t.shape)
			return solution,info
		except:
			LogUtil.exception_handler()


	def initialize_solver(self,solver_type,t0=0.0):
		"""Initialize an integrator."""
		try:
			assert solver_type in self.solver_list, 'Solver type {} is not available!'.format(self.solver_type)
			self.solver_type = solver_type

			if self.solver_type == 'ode-vode-bdf':
				self.ode_solver = ode(self.ODE_model_ode,self.jac_ODE_model_ode).set_integrator('vode',method='bdf',rtol=1e-4,atol=1e-4)
				self.ode_solver.set_initial_value(self.y0,t0)		
			
			elif self.solver_type == 'odeint':
				pass #Initialization not required if using odeint

			LogUtil.logger.debug('{}:{} solver initialized.'.format(self.name,self.solver_type))
		except:
			LogUtil.exception_handler()


	def ODE_model_ode(self,t,y):
		"""" Combine all derivatives when using ode method"""
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


	def jac_ODE_model_ode(self,t,y):
		""" Combine all derivatives."""
		try:
			y1 = y[0:self.PV_model.n_ODE]
			if self.PV_model.standAlone:
				self.grid_model.steady_state_model(t)
			y = self.PV_model.jac_ODE_model(y1,t)

			return y
		except:
			LogUtil.exception_handler()


	def check_simulation(self,infodict,t):
		"""Check whether the ODE solver failed at any time step."""
		try:
			infodict_mused = infodict['mused']
			if any(status == 0 or status > 2 for status in infodict_mused): #Alway check whether solver is failing to converge at any step
				assert len(infodict_mused) == len(t)-1," The number of methods used should be equal to the number of time steps."
				failure_time_point = list(map(lambda status: status == 0 or status > 2,infodict_mused)).index(True)
				self.SOLVER_CONVERGENCE = False
				self.convergence_failure_list.append({'Model':self.PV_model.name,
													'Simulation':self.name,
													'failure_time_point':t[failure_time_point],
													'failure_code':infodict_mused[failure_time_point],
													'S':self.PV_model.S*self.PV_model.Sbase})
					 
				error_message = '{}:ODE solver failed at {:.6f} s for {} with failure code:{}!'.format(self.name,t[failure_time_point],self.PV_model.name,infodict_mused[failure_time_point])
			
				states_at_failure = '\n___States at failure___\nVdc:{:.4f},Vta:{:.4f},Vpcca:{:.4f}\nia:{:.4f}\nPpv:{:.4f},S:{:.4f},\nma:{:.4f}'.format(self.PV_model.Vdc*self.PV_model.Vdcbase,
																																					 self.PV_model.vta*self.PV_model.Vbase,self.PV_model.va*self.PV_model.Vbase,self.PV_model.ia*self.PV_model.Ibase,self.PV_model.Ppv*self.PV_model.Sbase,self.PV_model.S*self.PV_model.Sbase,self.PV_model.ma)
				primary_controller_states = '\nxa:{:.4f},ua:{:.4f}'.format(self.PV_model.xa,self.PV_model.ua) 
			
				if self.DER_model_type in ['SolarPVDERSinglePhase', 'SolarPVDERThreePhase','SolarPVDERThreePhaseNumba','SolarPVDERThreePhaseBalanced']:
					secondary_controller_states = '\nxDC:{:.4f},xQ:{:.4f}'.format(self.PV_model.xDC,self.PV_model.xQ)
				elif self.DER_model_type == 'SolarPVDERSinglePhaseConstantVdc' or self.DER_model_type == 'SolarPVDERThreePhaseConstantVdc':
					secondary_controller_states = '\nxP:{:.4f},xQ:{:.4f}'.format(self.PV_model.xP,self.PV_model.xQ)
				else:
					raise ValueError('{} is not a valid DER model type!'.format(self.DER_model_type))
				raise ValueError(error_message+states_at_failure+primary_controller_states+secondary_controller_states)
			else:
				if self.DEBUG_SOLVER:
					infodict_time_steps = infodict['hu']
					assert len(infodict_mused) == len(infodict_time_steps) == len(t)-1," The number of methods used should be equal to the number of time steps."
					six.print_('{}:Time:{},Methods:{},Time steps:{}'.format(self.PV_model.name,t,infodict_mused,infodict_time_steps))
					if all(status == 1 or status == 2 for status in infodict_mused):
						six.print_('{}:Simulation successful for all time steps!'.format(self.PV_model.name))
				self.SOLVER_CONVERGENCE = True
		except:
			LogUtil.exception_handler()


