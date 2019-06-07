"""Manage simulation events."""

from __future__ import division
import operator
import six
import logging

import random
import math
import numpy as np

from pvder.utility_classes import Logging
from pvder import utility_functions

class SimulationEvents(Logging):
    """ Utility class for events."""
    
    events_count = 0
    Tactual_default = 298.15
    Zload1_actual_default = 10e6+0j
    
    _events_spec = {'insolation':{'default':100.0,'min':25.0,'max':100.0},
                   'voltage':{'default':1.0,'min':0.1,'max':1.1},  #This is a fraction and not per unit value
                   'frequency':{'default':60.0,'min':58.5,'max':61.5}}  #Time delay between events
    
    def __init__(self,SOLAR_EVENT_ENABLE = True,GRID_EVENT_ENABLE = True, LOAD_EVENT_ENABLE = True,verbosity='INFO'):
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
        
        self.initialize_logger()
        #Set logging level - {DEBUG,INFO,WARNING,ERROR}
        self.verbosity = verbosity
        
        self.SOLAR_EVENT_ENABLE = SOLAR_EVENT_ENABLE
        self.GRID_EVENT_ENABLE = GRID_EVENT_ENABLE
        self.LOAD_EVENT_ENABLE = LOAD_EVENT_ENABLE
        
        self.solar_events_list = []#{'T':3.0,'Sinsol':self.Sinsol_default,'Tactual':self.Tactual_default}]
        self.load_events_list = [] #{'T':4.0,'Zload1_actual':self.Zload1_actual_default}
        self.grid_events_list = [] #{'T':5.0,'Vgrms':self.Vgrms_default,'fgrid':self.fgrid_default}
        
        self.update_event_totals()
        self.reset_event_counters()
    
    def solar_events(self,t):
        """Generate solar event during simulation at specified time.
        
        Args:
           t (float): A scalar specifying the time (s).
        """
        
        if self.SOLAR_EVENT_ENABLE == True and self.solar_events_list: #Check whether list is empty
            
            if t<self.solar_events_list[0]['T']: 
                Sinsol = self._events_spec['insolation']['default']
                Tactual = self.Tactual_default
                
            elif t<self.solar_events_list[self.solar_event_counter]['T']  and self.solar_event_counter >=1:
                Sinsol = self.solar_events_list[self.solar_event_counter-1]['Sinsol']
                Tactual = self.solar_events_list[self.solar_event_counter-1]['Tactual']
            elif t>=self.solar_events_list[self.solar_event_counter]['T']:
                
                Sinsol = self.solar_events_list[self.solar_event_counter]['Sinsol']
                Tactual = self.solar_events_list[self.solar_event_counter]['Tactual']
                self.solar_event_counter = min(self.solar_events_total-1,self.solar_event_counter+1)
        else:
            Sinsol = self._events_spec['insolation']['default']
            Tactual = self.Tactual_default
        return Sinsol,Tactual
    
    def grid_events(self,t):
        """Generate grid event during simulation at specified time.
        
        Args:
           t (float): A scalar specifying the time (s).
        """
        
        if self.GRID_EVENT_ENABLE == True  and self.grid_events_list: #Check whether list is empty
            
            if t<self.grid_events_list[self.grid_event_counter]['T'] and self.grid_event_counter ==0:
                Vgrid = self._events_spec['voltage']['default']
                fgrid = self._events_spec['frequency']['default']
            elif t<self.grid_events_list[self.grid_event_counter]['T']  and self.grid_event_counter >=1:
                Vgrid = self.grid_events_list[self.grid_event_counter-1]['Vgrid']
                fgrid = self.grid_events_list[self.grid_event_counter-1]['fgrid']
            elif t>=self.grid_events_list[self.grid_event_counter]['T']:
                Vgrid = self.grid_events_list[self.grid_event_counter]['Vgrid']
                fgrid = self.grid_events_list[self.grid_event_counter]['fgrid']
                self.grid_event_counter = min(self.grid_events_total-1,self.grid_event_counter+1)
        else:
            Vgrid =self._events_spec['voltage']['default']
            fgrid =self._events_spec['frequency']['default']
        
        Vphase_angle_a = 0.0
        return Vgrid*pow(math.e,(1j*math.radians(Vphase_angle_a))),2.0*math.pi*fgrid
    
    def load_events(self,t):
        """Generate load event at PCC LV side during simulation.
        
        Args:
           t (float): A scalar specifying the time (s).
        """
        
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
           T (float): A scalar specifying start time of solar event in seconds.
           Sinsol (float): A scalar specifying solar insolation in percentage.
           Tactual (float): A scalar specifying module temperature in Kelvin.
        """        
        
        T = float(T)
        Sinsol = float(Sinsol)
        Tactual = float(Tactual)
        
        if Sinsol < self._events_spec['insolation']['min'] or Sinsol >  self._events_spec['insolation']['max']:
            raise ValueError('{} W/m2 is not a valid value for solar insolation! - Min:{},Max:{}'.format(Sinsol,self._events_spec['insolation']['min'],self._events_spec['insolation']['max']))
        
        if Tactual >375.0 or Tactual < 250.0:
            raise ValueError('{} K is not a valid value for temperature!'.format(Tactual))
        
        for event in self.solar_events_list:
            if T==event['T']:
                self.logger.debug('{}:Removing existing solar event at {:.2f}'.format(self.name,event['T']))
                self.solar_events_list.remove(event)   # in {}Remove exi,self.events_IDsting event at same time stamp
        
        self.logger.debug('{}:Adding new solar event at {:.2f} s'.format(self.name,T))
        self.solar_events_list.append({'T':T,'Sinsol':Sinsol,'Tactual':Tactual})  #Append new event to existing event list
        self.solar_events_list.sort(key=operator.itemgetter('T'))  #Sort new events list
        self.update_event_totals()
    
    def add_grid_event(self,T,Vgrid=1.0,fgrid=60.0):
        """Add new grid event.
        
        Args:
           T (float): A scalar specifying start time of grid event in seconds.
           Vgrid (float): A scalar specifying grid voltage in fraction.
           fgrid (float): A scalar specifying grid frequency in Hz.
        """
        
        if Vgrid < self._events_spec['voltage']['min'] or Vgrid >  self._events_spec['voltage']['max']:
            raise ValueError('{} V p.u. is not a valid value for grid voltage! - Min:{},Max:{}'.format(Vgrid,self._events_spec['voltage']['min'],self._events_spec['voltage']['max']))
        
        if fgrid < self._events_spec['frequency']['min'] or fgrid >  self._events_spec['frequency']['max']:
            raise ValueError('{} Hz is not a valid value for grid frequency! - Min:{},Max:{}'.format(fgrid,self._events_spec['frequency']['min'],self._events_spec['frequency']['max']))
            
        #assert Vgrms >=0.0 and Vgrms <=1.2 and fgrid >= 0.0 and fgrid <= 100.0, 'Grid event {} is not within feasible limits'.format(Vgrms,fgrid)
        
        T = float(T)
        Vgrid = float(Vgrid)
        fgrid = float(fgrid)
        
        for event in self.grid_events_list:
            if T==event['T']:
                self.logger.debug('{}:Removing existing grid event at {:.2f}'.format(self.name,event['T']))
                self.grid_events_list.remove(event)   #Remove existing event at same time stamp
        
        self.logger.debug('{}:Adding new grid event at {:.2f} s'.format(self.name,T))
        self.grid_events_list.append({'T':T,'Vgrid':Vgrid,'fgrid':fgrid})  #Append new event to existing event list
        self.grid_events_list.sort(key=operator.itemgetter('T'))  #Sort new events list
        self.update_event_totals()
    
    def add_load_event(self,T,Zload1_actual=10e6+0j):
        """Add new load event.
        
        Args:
           T: A scalar specifying start time of load event in seconds.
           Zload1_actual: A complex scalar specifying load in ohm.
        """
        
        T = float(T)
        if type(Zload1_actual) != complex:
            Zload1_actual = complex(Zload1_actual)
        
        if Zload1_actual.real < 0.0:
            raise ValueError('{} ohm is not a valid value for load resistance!'.format(Zload1_actual.real))
        
        for event in self.load_events_list:
            if T==event['T']:
                self.logger.debug('{}:Removing existing load event at {:.2f}'.format(self.name,event['T']))
                self.load_events_list.remove(event)   #Remove existing event at same time stamp
        
        self.logger.debug('{}:Adding new load event at {:.2f} s'.format(self.name,T))
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
                    logging.debug('{}:Solar event at {:.2f} s removed'.format(self.name,T))
                    self.solar_events_list.remove(event)
                    REMOVE_FLAG = True
            if not REMOVE_FLAG:
                logging.debug('{}:No solar event at {:.2f} s'.format(self.name,T))
        else:
            self.logger.debug('{}:Removing all events in solar events list and replacing with default event'.format(self.name))
            self.solar_events_list.clear()
            self.solar_events_list.append({'T':3.0,'Sinsol':self.Sinsol_default,'Tactual':self.Tactual_default})  #Defaut solar event
            
        self.update_event_totals() #Update total events
    
    def remove_grid_event(self,T):
        """Remove grid event at 'T'."""
        
        T = float(T)
        REMOVE_FLAG = False
        for event in self.grid_events_list:
            if event["T"] == T:
                self.logger.debug('{}:Grid event at {:.2f} s removed'.format(self.name,T))
                self.grid_events_list.remove(event)
                REMOVE_FLAG = True
        if not REMOVE_FLAG:
            self.logger.debug('{}:No grid event at {:.2f} s'.format(self.name,T))
        self.update_event_totals()
    
    def remove_load_event(self,T=None,REMOVE_ALL=False):
        """Remove solar event at 'T'"""
        
        if not REMOVE_ALL and T!=None:
            
            T = float(T)
            REMOVE_FLAG = False
            for event in self.load_events_list:
                if event["T"] == T:
                    self.logger.debug('{}:Load event at {:.2f} s removed'.format(self.name,event["T"]))
                    self.load_events_list.remove(event)
                    REMOVE_FLAG = True
            if not REMOVE_FLAG:
                self.logger.debug('No load event at {:.2f} s'.format(T)) 
        else:
            self.logger.debug('{}:Removing all events in load events list and replacing with default event'.format(self.name))
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
                if 'Vgrid' and 'fgrid' in event.keys():
                    six.print_('t:{:.3f},Grid event, Grid voltage is {:.2f} V, Frequency is {:.2f}'.format(event['T'],event['Vgrid'],event['fgrid']))
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
        
        logging.debug('{}:Simulation event counters reset!'.format(self.name))
    
    def create_random_events(self,t_event_start,t_event_end,t_event_step,events_type=['insolation','voltage']):
        """Create random events of specified types."""
        
        t_events = np.arange(t_event_start,t_event_end,t_event_step)
        
        for t_event in t_events:
            event_choice = random.choice(events_type)
            
            if event_choice == 'insolation':
                self.create_random_insolation_events(t_event)
            elif event_choice == 'voltage':
                self.create_random_voltage_events(t_event)
            else:
                print('{}:{} is not a valid event choice!'.format(self.name,event_choice))
    
    def create_random_insolation_events(self,t_event):
        """Create random voltage event at specified time."""
        
        insolation = self._events_spec['insolation']['min'] + random.random()*(self._events_spec['insolation']['max']-self._events_spec['insolation']['min'])
        self.add_solar_event(t_event,insolation)
    
    def create_random_voltage_events(self,t_event):
        """Create random voltage event at specified time."""
        
        voltage = self._events_spec['voltage']['min'] + random.random()*(self._events_spec['voltage']['max']-self._events_spec['voltage']['min']) 
        self.add_grid_event(t_event,voltage)    
    
    @property
    def simulation_events_list(self):
        """List of all simulation events."""
        
        return  sorted(self.solar_events_list + self.grid_events_list + self.load_events_list, key=operator.itemgetter('T'))  
   