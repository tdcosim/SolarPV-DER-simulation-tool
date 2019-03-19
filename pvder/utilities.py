from __future__ import division
import numpy as np
import math
import cmath
import sys
import six
from scipy.integrate import odeint
import time
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

class SimulationResults():
    """ Utility class for simulation results."""
    SAVE_PLOT_JPEG = False
    SAVE_PLOT_SVG = False
    figure_DPI = 1200

    def __init__(self,simulation,figure_index=1,PER_UNIT=True,font_size=18,PLOT_TITLE=True):
        # do a lazy import. This is against PEP style guidelines. However,
        # matplotlib has a ~25mb memory foot print. Hence, if each opendssapi.worker
        # imports this then memory consumption becomes a burden very quickly. This
        # way only when simulationresults class is called, matplotlib will be
        # imported.
        #import matplotlib.pyplot as plt
        #import matplotlib.ticker as ticker

        self.figure_index =1
        self.simulation = simulation
        
        self.PER_UNIT = PER_UNIT
        self.PLOT_TITLE = PLOT_TITLE
        self.font_size = font_size

    def change_units(self):
        """Change units from per unit to S.I. or vice versa."""
        
        if self.PER_UNIT == True:
           self.V_multiplier = 1.0
           self.I_multiplier = 1.0
           self.S_multiplier = 1.0
           six.print_('Using per unit quantities!')
        else:
           self.V_multiplier = self.simulation.grid_model.Vbase
           self.I_multiplier = self.simulation.grid_model.Ibase
           self.S_multiplier = self.simulation.grid_model.Sbase
           six.print_('Using SI quantities!')
    
    def group_quantities_for_plotting(self,plot_type):
        """Plot power from simulation."""
        
        if plot_type not in {'power','active_power','active_power_Ppv_Pac_PCC','active_power_Pac_PCC','reactive_power','reactive_power_Q_PCC','reactive_power_Q_PCC_smoothed','voltage','voltage_Vdc','voltage_HV','voltage_HV_imbalance','voltage_LV','voltage_Vpcclv','voltage_Vpcclv_smoothed','current','phase_angle','frequency'}:
            raise ValueError('Unknown plot type: ' + str(plot_type))
        
        time = self.simulation.t
        self.change_units()
        if self.PER_UNIT == True:
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
            plot_values = [self.simulation.Vdc_t*self.V_multiplier,self.simulation.Vtrms_t*self.V_multiplier,self.simulation.Vrms_t*self.V_multiplier,self.simulation.Vhvrms_t*self.V_multiplier,self.simulation.Vgrms_t*self.V_multiplier,self.simulation.vdg_t*self.V_multiplier]
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
        
        elif plot_type == 'voltage_HV': 
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
 self.simulation.S_PCC_t.real*self.S_multiplier,self.simulation.S_PCC_t.imag*self.S_multiplier,self.simulation.S_load1_t.real*self.S_multiplier,self.simulation.S_G_t.real*self.S_multiplier]
            #legends=['Ppv','Pac','Q','Pac_PCC','Q_PCC','Pac_load1','Pac_G']
            
            legends=[r"$P^{PV}_{DC}$",r"$P^{inv}_{AC}$",r"$Q^{inv}$",r"$P^{PCC-LV}_{AC}$",r"$Q^{PCC-LV}$",r"$P^{load1}_{AC}$",r"$P^{G}_{AC}$"] 
            plot_title='All available power quantities!'
            y_labels=_power_label
        
        elif plot_type == 'active_power':
            plot_values = [self.simulation.Ppv_t*self.S_multiplier,self.simulation.S_t.real*self.S_multiplier,self.simulation.S_PCC_t.real*self.S_multiplier,self.simulation.S_load1_t.real*self.S_multiplier,self.simulation.S_G_t.real*self.S_multiplier]
            
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
            plot_values = [self.simulation.S_t.imag*self.S_multiplier,self.simulation.S_PCC_t.imag*self.S_multiplier,self.simulation.S_load1_t.imag*self.S_multiplier,self.simulation.S_G_t.imag*self.S_multiplier]
            
            #legends=['Q','Q_PCC','Q_load1','Q_G']
            
            legends=[r"$Q^{inv}$",r"$Q^{PCC-LV}$",r"$Q^{load1}$",r"$Q^{G}$"] 
            plot_title='Inverter,PCC LV,Load,& Grid voltage:Reactive power'
            y_labels=_reactive_power_label
        
        elif plot_type == 'reactive_power_Q_PCC':
            plot_values = [self.simulation.S_PCC_t.imag*self.S_multiplier]
            #legends=['Q_PCC']
            
            legends=[r"$Q^{PCC-LV}$"]
            plot_title='PCC-LV side:Reactive power'
            y_labels=_reactive_power_label
    
        elif plot_type == 'phase_angle':
            plot_values = [self.simulation.phi_at_t,self.simulation.phi_a_t,self.simulation.phi_ag_t,self.simulation.phi_a_t-self.simulation.phi_ag_t]
            #legends=['theta_vat','theta_va','theta_vag','theta_va-theta_vag']
            
            legends=[r"$\phi^{inv}_{a}$",r"$\phi^{PCC-LV}_{a}$",r"$\phi^{grid}_{a}$",r"$\phi^{PCC-LV}_{a}-\phi^{grid}_{a}$"]
            plot_title='All available phase angle quantities!'
            y_labels='radians'
        
        elif plot_type == 'frequency':
            plot_values = [self.simulation.wgrid_t,self.simulation.we_t]
            #legends=['Grid frequency','PLL frequency']
            
            legends=[r"$\omega_{grid}$",r"$\omega_{PLL}$"]
            plot_title='Grid voltage source and PLL: Frequency'
            y_labels='radians/s'
        
        return time,plot_values,legends,plot_title,y_labels

    def plot_DER_simulation(self,plot_type = 'power'):
        """Plot desired electrical quantity from simulation."""
        time,plot_values,legends,plot_title,y_labels = self.group_quantities_for_plotting(plot_type)
        time_values = [time]*len(plot_values)
        
        self.plot_multiple(time_values,plot_values,legends,plot_title,y_labels)
        self.save_plot(plt,plot_type)
        self.figure_index =  self.figure_index+1
        
    def plot_multiple(self,time_values,plot_values,legends,plot_title,y_labels):
        """Function to plot multiple time series on same plot."""
        
        assert len(plot_values) == len(time_values) == len(legends),  " The number of legends should be equal to the number of quantities"
        fig = plt.figure(self.figure_index, figsize=(8,8))
        
        for i in range(len(plot_values)):
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
    
    def save_plot(self,plot_object,plot_name='results'):
        """Save the plots."""
        
        if self.SAVE_PLOT_JPEG == True:
           plot_object.savefig(plot_name+".jpg", dpi=self.figure_DPI)
        if self.SAVE_PLOT_SVG == True:
           plot_object.savefig(plot_name+".svg", dpi=self.figure_DPI)
    
    def compare_with_external(self,external_time_values,external_plot_values,external_plot_legends,plot_type='power'):
        """Compare with external simulation."""
        
        assert external_plot_values != None, 'There should be external values to plot'
        assert len(external_plot_values) == len(external_plot_legends),'Legends should be equal to number of plots'
        
        time,internal_plot_values,internal_legends,plot_title,y_labels = self.group_quantities_for_plotting(plot_type)
        internal_time_values = [time]*len(internal_plot_values)
                
        time_values = internal_time_values+external_time_values
        plot_values = internal_plot_values + external_plot_values
        legends = internal_legends + external_plot_legends
        self.plot_multiple(time_values,plot_values,legends,plot_title,y_labels)
        self.figure_index =  self.figure_index+1
    
class SimulationUtilities():
    """ Utility class for dynamic simulations."""
    
    def call_ODE_solver(self,derivatives,jacobian,y,t):
        """Call the SciPy ODE solver."""
        
        if self.jacFlag:
            solution,infodict = odeint(derivatives,y,t,Dfun=jacobian,full_output=1,printmessg=True,\
                                       hmax = 1/120.,mxstep=50,atol=1e-4,rtol=1e-4)
        else:
            solution,infodict =  odeint(derivatives,y,t,full_output=1,printmessg=True,\
                                       hmax = 1/120.,mxstep=50,atol=1e-4,rtol=1e-4)
        #return odeint(derivatives,y,t,full_output=1,printmessg=True,\
        #hmax =1/120.,mxstep=1000)#atol=1.49012e-6
        
        self.check_simulation(infodict,t) #Check whether solver successful for all time intervals
        
        return solution,infodict
    
    def solve_ODE(self):
        """Solve given system of ODE's."""
        
        timer_start = time.time()
        six.print_("Simulation for {} started at {} s and will end at {} s".format(self.name,self.tStart,self.tStop))
        if self.jacFlag:
           six.print_("Analytical Jacobian will be provided to ODE solver.")
        
        solution,infodict  = self.call_ODE_solver(self.ODE_model,self.jac_ODE_model,self.y0,self.t)
        six.print_("Simulation Finished")
        six.print_('Total simulation time is {}'.format(time.strftime("%H:%M:%S", time.gmtime(time.time()-timer_start))))
        
        return solution
    
    def check_simulation(self,infodict,t):
        """Check whether the ODE solver failed at any time step."""
        infodict_mused = infodict['mused']
        
        if any(status == 0 or status > 2 for status in infodict_mused): #Alway check whether solver is failing to converge at any step
                 assert len(infodict_mused) == len(t)-1,  " The number of methods used should be equal to the number of time steps."
                 failure_time_point = list(map(lambda status: status == 0 or  status > 2,  infodict_mused)).index(True)
                 
                 raise ValueError('ODE solver failed at {:.4f} s for {} with failure code:{}!\n___States at failure___\nVdc:{:.4f},Vta:{:.4f},Vpcca:{:.4f}\nPpv:{:.4f},S:{:.4f}'
                                  .format(t[failure_time_point],self.PV_model.name,infodict_mused[failure_time_point],self.PV_model.Vdc*self.PV_model.Vdcbase,
                                  self.PV_model.vta*self.PV_model.Vbase,self.PV_model.va*self.PV_model.Vbase,self.PV_model.Ppv*self.PV_model.Sbase,self.PV_model.S*self.PV_model.Sbase))
        
        elif self.DEBUG_SOLVER:
            infodict_time_steps = infodict['hu']
            assert len(infodict_mused) == len(infodict_time_steps) == len(t)-1,  " The number of methods used should be equal to the number of time steps."
            six.print_('Time:{},Methods:{},Time steps:{}'.format(t,infodict_mused,infodict_time_steps))
            
            if all(status == 1 or status == 2 for status in infodict_mused):
               six.print_('Simulation successful for all time steps for {}!'.format(self.PV_model.name))
