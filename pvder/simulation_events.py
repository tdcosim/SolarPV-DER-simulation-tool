"""Manage simulation events."""

from __future__ import division
import math
import operator
import six
import numpy as np
from pvder import utility_functions

class SimulationEvents():
    """ Utility class for events."""
    
    events_count = 0
    Vgrms_default = 1.0  #This is a fraction and not per unit value
    fgrid_default = 60.0 
    Sinsol_default = 100.0
    Tactual_default = 298.15
    Zload1_actual_default = 10e6+0j
    
    def __init__(self,SOLAR_EVENT_ENABLE = True,GRID_EVENT_ENABLE = True, LOAD_EVENT_ENABLE = True):
        """Creates an instance of `SimulationEvents`.
        
        Args:
          SOLAR_EVENT_ENABLE: A boolean to enable solar insolation events.
          GRID_EVENT_ENABLE: A boolean to enable grid voltage or frequency events.
          LOAD_EVENT_ENABLE: A boolean to enable load change events at PCC-LV side.
        """
                
        #Increment count to keep track of number of simulation events instances
        SimulationEvents.events_count = SimulationEvents.events_count + 1
        self.events_ID =SimulationEvents.events_count
        #Object name
        self.name = 'events_'+str(self.events_ID)
        
        self.SOLAR_EVENT_ENABLE = SOLAR_EVENT_ENABLE
        self.GRID_EVENT_ENABLE = GRID_EVENT_ENABLE
        self.LOAD_EVENT_ENABLE = LOAD_EVENT_ENABLE
        
        self.solar_events_list = [{'T':3.0,'Sinsol':self.Sinsol_default,'Tactual':self.Tactual_default}]
        self.load_events_list = [{'T':4.0,'Zload1_actual':self.Zload1_actual_default}]
        self.grid_events_list = [{'T':5.0,'Vgrms':self.Vgrms_default,'fgrid':self.fgrid_default}]
        
        self.update_event_totals()
        self.reset_event_counters()
    
    def solar_events(self,t):
        """List of all simulation events."""
        
        if self.SOLAR_EVENT_ENABLE == True and self.solar_events_list: #Check whether list is empty
            
            if t<self.solar_events_list[0]['T']: 
                Sinsol = self.Sinsol_default
                Tactual = self.Tactual_default
                
            elif t<self.solar_events_list[self.solar_event_counter]['T']  and self.solar_event_counter >=1:
                Sinsol = self.solar_events_list[self.solar_event_counter-1]['Sinsol']
                Tactual = self.solar_events_list[self.solar_event_counter-1]['Tactual']
            elif t>=self.solar_events_list[self.solar_event_counter]['T']:
                
                Sinsol = self.solar_events_list[self.solar_event_counter]['Sinsol']
                Tactual = self.solar_events_list[self.solar_event_counter]['Tactual']
                self.solar_event_counter = min(self.solar_events_total-1,self.solar_event_counter+1)
        else:
            Sinsol = self.Sinsol_default
            Tactual = self.Tactual_default
        return Sinsol,Tactual
    
    def grid_events(self,t):
        """Generate events during simulation."""
        
        if self.GRID_EVENT_ENABLE == True  and self.grid_events_list: #Check whether list is empty
            
            if t<self.grid_events_list[self.grid_event_counter]['T'] and self.grid_event_counter ==0:
                Vgrms = self.Vgrms_default
                fgrid = self.fgrid_default
            elif t<self.grid_events_list[self.grid_event_counter]['T']  and self.grid_event_counter >=1:
                Vgrms = self.grid_events_list[self.grid_event_counter-1]['Vgrms']
                fgrid = self.grid_events_list[self.grid_event_counter-1]['fgrid']
            elif t>=self.grid_events_list[self.grid_event_counter]['T']:
                Vgrms = self.grid_events_list[self.grid_event_counter]['Vgrms']
                fgrid = self.grid_events_list[self.grid_event_counter]['fgrid']
                self.grid_event_counter = min(self.grid_events_total-1,self.grid_event_counter+1)
        else:
            Vgrms = self.Vgrms_default
            fgrid = self.fgrid_default
        
        Vphase_angle_a = 0.0
        return Vgrms*pow(math.e,(1j*math.radians(Vphase_angle_a))),2.0*math.pi*fgrid
    
    def load_events(self,t):
        """Generate load events at PCC LV side during simulation."""
        
        if self.LOAD_EVENT_ENABLE == True and self.load_events_list: #Check whether list is empty
            
            if t<self.load_events_list[0]['T']: 
                Zload1_actual = self.Zload1_actual_default
                
            elif t<self.load_events_list[self.load_event_counter]['T']  and self.load_event_counter >=1:
                Zload1_actual = self.load_events_list[self.load_event_counter-1]['Zload1_actual']
                
            elif t>=self.load_events_list[self.load_event_counter]['T']:
                Zload1_actual = self.load_events_list[self.load_event_counter]['Zload1_actual']
                self.load_event_counter = min(self.load_events_total-1,self.load_event_counter+1)
        else:
            Zload1_actual = self.Zload1_actual_default
           
        return Zload1_actual
    
    def add_solar_event(self,T,Sinsol=100.0,Tactual=298.15):
        """Add new solar event.
        
        Args:
           T: A scalar specifying start time of solar event in seconds.
           Sinsol: A scalar specifying solar insolation in percentage.
           Tactual: A scalar specifying module temperature in Kelvin.
        """        
        
        T = float(T)
        Sinsol = float(Sinsol)
        Tactual = float(Tactual)
        
        if Sinsol >150.0 or Sinsol < 0.0:
           raise ValueError('{} W/m2 is not a valid value for solar insolation!'.format(Sinsol))
        if Tactual >375.0 or Tactual < 250.0:
           raise ValueError('{} K is not a valid value for temperature!'.format(Tactual))
        
        for event in self.solar_events_list:
            if T==event['T']:
                six.print_('Removing existing solar event at {:.2f}'.format(event['T']))
                self.solar_events_list.remove(event)   # in {}Remove exi,self.events_IDsting event at same time stamp
        
        six.print_('Adding new solar event at {:.2f} s'.format(T))
        self.solar_events_list.append({'T':T,'Sinsol':Sinsol,'Tactual':Tactual})  #Append new event to existing event list
        self.solar_events_list.sort(key=operator.itemgetter('T'))  #Sort new events list
        self.update_event_totals()
    
    def add_grid_event(self,T,Vgrms=1.0,fgrid=60.0):
        """Add new solar event.
        
        Args:
           T: A scalar specifying start time of grid event in seconds.
           Vgrms: A scalar specifying grid voltage in fraction.
           fgrid: A scalar specifying grid frequency in Hz.
        """
        
        assert Vgrms >=0.0 and Vgrms <=1.2 and fgrid >= 0.0 and fgrid <= 100.0, 'Grid event {} is not within feasible limits'.format(Vgrms,fgrid)
        
        T = float(T)
        Vgrms = float(Vgrms)
        fgrid = float(fgrid)
        for event in self.grid_events_list:
            if T==event['T']:
                six.print_('Removing existing grid event at {:.2f}'.format(event['T']))
                self.grid_events_list.remove(event)   #Remove existing event at same time stamp
        
        six.print_('Adding new grid event at {:.2f} s'.format(T))
        self.grid_events_list.append({'T':T,'Vgrms':Vgrms,'fgrid':fgrid})  #Append new event to existing event list
        self.grid_events_list.sort(key=operator.itemgetter('T'))  #Sort new events list
        self.update_event_totals()
    
    def add_load_event(self,T,Zload1_actual=10e6+0j):
        """Add new load event.
        
        Args:
           T: A scalar specifying start time of load event in seconds.
           Vgrms: A complex scalar specifying load in ohm.
        """
        
        T = float(T)
        if type(Zload1_actual) != complex:
            Zload1_actual = complex(Zload1_actual)
        
        if Zload1_actual.real < 0.0:
            raise ValueError('{} ohm is not a valid value for load resistance!'.format(Zload1_actual.real))
        
        for event in self.load_events_list:
            if T==event['T']:
                six.print_('Removing existing load event at {:.2f}'.format(event['T']))
                self.load_events_list.remove(event)   #Remove existing event at same time stamp
        
        six.print_('Adding new load event at {:.2f} s'.format(T))
        self.load_events_list.append({'T':T,'Zload1_actual':Zload1_actual})  #Append new event to existing event list
        self.load_events_list.sort(key=operator.itemgetter('T'))  #Sort new events list
        self.update_event_totals()
    
    def remove_solar_event(self,T=None,REMOVE_ALL=False):
        """Remove solar event at 'T'."""
        
        if not REMOVE_ALL and T!=None:
            T = float(T)
            REMOVE_FLAG = False
            for event in self.solar_events_list:
                if event["T"] == T:
                    six.print_('Solar event at {:.2f} s removed'.format(T))
                    self.solar_events_list.remove(event)
                    REMOVE_FLAG = True
            if REMOVE_FLAG == False:
                six.print_('No solar event at {:.2f} s'.format(T))
        else:
            six.print_('Removing all events in solar events list and replacing with default event')
            self.solar_events_list.clear()
            self.solar_events_list.append({'T':3.0,'Sinsol':self.Sinsol_default,'Tactual':self.Tactual_default})  #Defaut solar event
            
        self.update_event_totals() #Update total events
    
    def remove_grid_event(self,T):
        """Remove grid event at 'T'."""
        
        T = float(T)
        REMOVE_FLAG = False
        for event in self.grid_events_list:
           if event["T"] == T:
              six.print_('Grid event at {:.2f} s removed'.format(T))
              self.grid_events_list.remove(event)
              REMOVE_FLAG = True
        if REMOVE_FLAG == False:
           six.print_('No grid event at {:.2f} s'.format(T))
        self.update_event_totals()
    
    def remove_load_event(self,T=None,REMOVE_ALL=False):
        """Remove solar event at 'T'"""
        
        if not REMOVE_ALL and T!=None:
            
            T = float(T)
            REMOVE_FLAG = False
            for event in self.load_events_list:
                if event["T"] == T:
                    six.print_('Load event at {:.2f} s removed'.format(event["T"]))
                    self.load_events_list.remove(event)
                    REMOVE_FLAG = True
            if REMOVE_FLAG == False:
                six.print_('No load event at {:.2f} s'.format(T)) 
        else:
            six.print_('Removing all events in load events list and replacing with default event')
            self.load_events_list.clear()
            self.load_events_list.append({'T':5.0,'Zload1_actual':self.Zload1_default})  #Defaut load event
            
        self.update_event_totals()
    
    def insolation_ramp(self,tstart,tstop,Sinsol_target,tstep=0.25):
        """Create a ramp signal.
        
        Args:
           tstart: A scalar specifying start time of solar insolation event in seconds.
           tstop: A scalar specifying stop time of solar insolation event in seconds.
           Sinsol_target: A scalar specifying target solar insolation at 'tstop' in percentage.
           tstep: A scalar specifying time step size.
        """
        
        Sinsol,_ = self.solar_events(t=tstart)
        
        trange = np.arange(tstart,tstop,tstep)
        Sinsol_range = np.linspace(Sinsol,Sinsol_target,len(trange))
        for i,Sinsol in enumerate(Sinsol_range):
            self.add_solar_event(T=trange[i],Sinsol=Sinsol)
    
    def voltage_ramp(self,tstart,tstop,vg_target,tstep=0.5):
        """Create a ramp signal for grid voltage.
        
        Args:
           tstart: A scalar specifying start time of grid voltage event in seconds.
           tstop: A scalar specifying stop time of grid voltage event in seconds.
           vg_target: A scalar specifying target grid voltage at 'tstop' in fraction.
           tstep: A scalar specifying time step size (optional).
        """
        
        vg,_ = self.grid_events(t=tstart)
        
        trange = np.arange(tstart,tstop,tstep)
        vg_range = np.linspace(vg,vg_target,len(trange))
        for i,vg in enumerate(vg_range):
            self.add_grid_event(T=trange[i],Vgrms=vg) 
    
    def show_events(self):
        """Print all the simulation events."""
        
        print('Showing all event in events instance {}'.format(self.name ))
        print('Total solar events:{}\nTotal grid events:{}'.format(len(self.solar_events_list),len(self.grid_events_list)))
        if self.simulation_events_list:
           
           for event in self.simulation_events_list:
                if 'S' and 'Tactual' in event.keys():
                   six.print_('t:{:.3f},Solar event, Solar insolation is {:.2f} W/cm2, Temperature is {:.2f}'.format(event['T'],event['Sinsol'],event['Tactual']))
                if 'Vgrms' and 'fgrid' in event.keys():
                   six.print_('t:{:.3f},Grid event, Grid voltage is {:.2f} V, Frequency is {:.2f}'.format(event['T'],event['Vgrms'],event['fgrid']))
                if 'Zload1_actual' in event.keys():
                    six.print_('t:{:.3f},Load event, Impedance is {:.2f} ohm'.format(event['T'],event['Zload1_actual']))
        else:
            six.print_("No simulation events!!!")
    
    def update_event_totals(self):
        """Update event counts. """
        
        self.solar_events_total = len(self.solar_events_list)
        self.grid_events_total = len(self.grid_events_list)
        self.load_events_total = len(self.load_events_list)
        self.events_total =self.solar_events_total +self.grid_events_total  +  self.load_events_total 
    
    def reset_event_counters(self):
        """Reset event counts. """
        
        self.solar_event_counter = 0
        self.grid_event_counter = 0
        self.load_event_counter = 0
        
        six.print_('Simulation event counters reset!')
        
    @property
    def simulation_events_list(self):
        """List of all simulation events."""
        
        return  sorted(self.solar_events_list + self.grid_events_list + self.load_events_list, key=operator.itemgetter('T'))  
   