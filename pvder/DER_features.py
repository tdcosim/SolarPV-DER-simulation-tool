"""Code for features inside PV-DER model instances."""

from __future__ import division
import six

import math

from pvder import utility_functions
from pvder import defaults,templates,specifications

class PVDER_SmartFeatures():
    """Class for describing smart inverter inverter features of PV-DER."""
    
    #Flags
    VOLT_VAR_ENABLE = False
    VOLT_VAR_FLAG = False
    
    VOLT_WATT_ENABLE = False
    VOLT_WATT_FLAG = False    
    
    V_threshold_high_limit = defaults.Vthreshold_high_limit #percentage
    t_threshold_high_limit = defaults.tthreshold_high_limit #seconds
    t_disconnect_low_limit = defaults.tdisconnect_low_limit #seconds #Minimum time to initiate DER disconnection (output cessation)
    t_reconnect_low_limit = defaults.treconnect_low_limit #seconds #Minimum time to initiate DER reconnection (output restoration)
    
    f_ref = 60.0
    DER_CONNECTED = True
    DER_MOMENTARY_CESSATION = False
    DER_TRIP     = False        
    
    def initialize_Volt_VAR(self):
        """Initialize the Volt-VAR controller settings."""
        
        self.Volt_VAR_dict = {'1':{'V':0.92,
                                   'Q':0.44},
                              '2':{'V':0.98,
                                   'Q':0.0},
                              '3':{'V':1.02,
                                   'Q':0.0},
                              '4':{'V':1.08,
                                   'Q':-0.44}
                              }  #From IEEE 1547-2018 Catergy B (Table 8 - page 39)        

        self.Volt_VAR_dict['m_inject'] = (self.Volt_VAR_dict['2']['Q'] - self.Volt_VAR_dict['1']['Q'])/(self.Volt_VAR_dict['2']['V'] - self.Volt_VAR_dict['1']['V'])
        
        self.Volt_VAR_dict['m_absorb'] = (self.Volt_VAR_dict['4']['Q'] - self.Volt_VAR_dict['3']['Q'])/(self.Volt_VAR_dict['4']['V'] - self.Volt_VAR_dict['3']['V'])
        
        self.Volt_VAR_dict['c_inject'] = self.Volt_VAR_dict['2']['Q'] - self.Volt_VAR_dict['2']['V']*self.Volt_VAR_dict['m_inject']
        self.Volt_VAR_dict['c_absorb'] = self.Volt_VAR_dict['4']['Q'] - self.Volt_VAR_dict['4']['V']*self.Volt_VAR_dict['m_absorb']
        
        self.VOLT_VAR_ACTIVE = False
        
        self.Qlimit = self.Qlimit_calc()
        self.Volt_VAR_logic_t = 0.0
        
    def Volt_VAR_logic(self,t):
        """ Volt-VAR."""
        
        #Select RMS voltage source
        Vrms_measured = self.Vrms
        self.Qlimit = self.Qlimit_calc()
        
        deadband_V = 0.01
        
        if t>self.Volt_VAR_logic_t: #Update volt-var logic only if solver time is greater that current time
            self.Volt_VAR_logic_t = t  
       
            if self.VOLT_VAR_ACTIVE:
                
                if self.Volt_VAR_inject_range(Vrms_measured,deadband_V):
                    Qref = self.Q_Volt_VAR_inject(Vrms_measured)                    
                elif self.Volt_VAR_absorb_range(Vrms_measured,deadband_V):
                    Qref = self.Q_Volt_VAR_absorb(Vrms_measured)       
                else:
                    utility_functions.print_to_terminal("Volt-VAR control with reference voltage {:.3f} is deactivated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,Vrms_measured*self.Vbase))
                    Qref = 0.0
                    self.VOLT_VAR_ACTIVE = False                    
            else:
                if self.Volt_VAR_inject_range(Vrms_measured,deadband_V):
                    Qref = self.Q_Volt_VAR_inject(Vrms_measured)                    
                    utility_functions.print_to_terminal("Volt-VAR control (injection) with reference voltage {:.3f} is activated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,Vrms_measured*self.Vbase))
                    self.VOLT_VAR_ACTIVE = True
                elif self.Volt_VAR_absorb_range(Vrms_measured,deadband_V):
                    Qref = self.Q_Volt_VAR_absorb(Vrms_measured)
                    utility_functions.print_to_terminal("Volt-VAR control (absorbption) with reference voltage {:.3f} is activated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,Vrms_measured*self.Vbase))            
                    self.VOLT_VAR_ACTIVE = True

                else: #Don't supply/absorb reactive power outside operating range
                    Qref = 0.0
        else:
            self.Volt_VAR_logic_t = self.Volt_VAR_logic_t
            Qref = self.Q_ref            
                    
        return Qref
    
    def Qlimit_calc(self):
        """Calculate maximum Q reference."""        
        
        Qmax = math.sqrt(math.pow(self.Sinverter_nominal,2)-math.pow(max(-self.Sinverter_nominal,min(self.S.real,self.Sinverter_nominal)),2))    #Qmax = sqrt(S^2 - P^2)
        
        Qlimit = Qmax - self.n_phases*self.Xf*(abs(self.ia)/math.sqrt(2))**2  #Qlimit = Qmax - Qfilter
        Qlimit = max(0.0,Qlimit)
        
        return Qlimit
    
    def Volt_VAR_inject_range(self,Vrms_measured,del_V=0.0):
        """Check if voltage in volt-VAR operating range."""
        
        return Vrms_measured >= (self.Volt_VAR_dict['1']['V'] - del_V)*self.Vrms_ref and \
               Vrms_measured <= (self.Volt_VAR_dict['2']['V'] + del_V)*self.Vrms_ref   
    
    def Volt_VAR_absorb_range(self,Vrms_measured,del_V=0.0):
        """Check if voltage in volt-VAR operating range."""
        
        return Vrms_measured >= (self.Volt_VAR_dict['3']['V'] - del_V)*self.Vrms_ref and \
               Vrms_measured <= (self.Volt_VAR_dict['4']['V'] + del_V)*self.Vrms_ref       
    
    def Q_Volt_VAR_inject(self,Vrms_measured):
        """Calculate reactive power for Volt-VAR control."""
        
        if self.Volt_VAR_inject_range(Vrms_measured,del_V=0.0):
            Vrms_measured = (Vrms_measured/self.Vrms_ref)
            Q = max(0.0,((Vrms_measured)*self.Volt_VAR_dict['m_inject'] + self.Volt_VAR_dict['c_inject'])*self.Sinverter_nominal)
            Q = min(self.Qlimit,Q)
        else:
            Q = self.Q_ref            
        
        return Q
    
    def Q_Volt_VAR_absorb(self,Vrms_measured):
        """Calculate reactive power for Volt-VAR control."""
        
        if self.Volt_VAR_absorb_range(Vrms_measured,del_V=0.0):
            Vrms_measured = (Vrms_measured/self.Vrms_ref)
            Q = min(0.0,((Vrms_measured)*self.Volt_VAR_dict['m_absorb'] + self.Volt_VAR_dict['c_absorb'])*self.Sinverter_nominal)
            Q = max(-self.Qlimit,Q)        
        else:
            Q = self.Q_ref            
        
        return Q
    
    def update_ridethrough_flags(self,t):
        """Check VRT and FRT logic."""
        
        if t < self.RT_logic_t_lock:
            self.logger.debug('{}:{:.4f}:Time going backwards in VRT/FRT logic - previous time was {:.4f} s.'.format(self.name,t,self.RT_logic_t_lock))
        
        elif t > self.RT_logic_t_lock: #Go through logic only if time is greater than previous time
            self.RT_logic_t_lock = t
            #LVRT logic        
            if self.LVRT_ENABLE:
                self.LVRT(t)

            #HVRT logic
            if self.HVRT_ENABLE:
                self.HVRT(t)

            #LFRT trip logic
            if self.LFRT_ENABLE:
                self.FRT(t)

            self.DER_TRIP = self.LVRT_TRIP or self.HVRT_TRIP# or self.LFRT_TRIP
            self.DER_MOMENTARY_CESSATION = self.LVRT_MOMENTARY_CESSATION or self.HVRT_MOMENTARY_CESSATION   
        
    def disconnect_or_reconnect(self,t):
        """Check flags and ither disconnect or reconnect DER."""
        
        if t < self.connect_logic_t_lock:
            self.logger.debug('{}:{:.4f}:Time going backwards in disconnect/reconnect logic - previous time was {:.4f} s .'.format(self.name,t,self.connect_logic_t_lock))
            
        elif t > self.connect_logic_t_lock: #Go through logic only if time is greater than previous time
            self.connect_logic_t_lock = t
            if self.DER_CONNECTED:
                if self.DER_TRIP or self.DER_MOMENTARY_CESSATION:
                    self.DER_disconnect_logic(t)       #Update status of DER_CONNECTED
            else:
                if not self.DER_TRIP:
                    self.DER_reconnect_logic(t)
        
        if not self.DER_CONNECTED: #Keep DER connected only if DER_CONNECTED is True
            self.DER_disconnect()
        
    def DER_disconnect_logic(self,t):
        """Logic for disconnecting/deenergizing DER from grid."""
        
        assert self.DER_CONNECTED, 'Disconnection logic can only be used if DER is already connected.'
        
        if self.t_disconnect_start == 0.0: #Start disconnect timer
            text_string = '{}:{:.4f}:DER disconnect timer started.'.format(self.name,t)
            print(text_string)
            self.print_event(text_string,False) 
            self.t_disconnect_start = t
        elif t-self.t_disconnect_start < self.t_disconnect_delay: #Disconnect DER only after disconnect time delay has elapsed
            text_string = '{}:{:.4f}:DER is in disconnect timer zone.'.format(self.name,t)
            self.logger.debug(text_string)
        elif t-self.t_disconnect_start >= self.t_disconnect_delay: #Disconnect DER only after disconnect time delay has elapsed
            text_string = '{}:{:.4f}:DER will be disconnected.'.format(self.name,t)
            print(text_string)
            self.print_event(text_string,False) 
            self.t_disconnect_start = 0.0
            self.DER_CONNECTED = False
    
    def DER_disconnect(self):
        """Function to disconnect PV-DER from grid."""
        
        #self.logger.debug('{}:DER disconnected'.format(self.name))
       
        self.VOLT_VAR_ENABLE = False
        self.VOLT_WATT_ENABLE = False
        
        self.Q_ref = 0.0  #Set reactive power reference to zero
        self.Vdc_ref = self.Vdc  #Maintain DC link voltage
        self.Ppv = 0.0     #Disconnect PV panel
        
        self.ia_ref = 0.0 + 1j*0.0
        self.ib_ref = 0.0 + 1j*0.0
        self.ic_ref = 0.0 + 1j*0.0        
    
    def DER_reconnect_logic(self,t):
        """Logic used to decide reconnection."""
        
        #Select RMS voltage source
        Vrms_measured = self.Vrms   #Select PCC - LV side voltage   
        fgrid = self.we/(2.0*math.pi)  #Use grid frequency as estimated by PLL
        
        assert not self.DER_CONNECTED, 'Reconnection logic can only be used if DER is disconnected.'
        
        if self.DER_MOMENTARY_CESSATION:
            if self.t_reconnect_start > 0.0:
                self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='reconnect_reset')
                self.t_reconnect_start = 0.0
            elif self.t_reconnect_start == 0.0:
                self.print_reconnect_events(t, Vrms_measured,fgrid,event_name='DER_tripped')    #Inverter output remains zero
                
        elif t > 0.0:
            if self.t_reconnect_start > 0.0:
                if t-self.t_reconnect_start < self.t_reconnect_delay:
                    self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='reconnect_zone')
            
                if t-self.t_reconnect_start > self.t_reconnect_delay:
                    self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='DER_reconnection')
                    self.DER_CONNECTED = True
                    self.t_reconnect_start = 0.0
                    if self.RESTORE_Vdc: #Whether to restore Vdc to pre-cessation setpoint
                        self.Vdc_ref_ramp(t,self.set_Vdc_ref()*self.Vbase)
                        self.logger.info('Ramping Vdc to pre-anomaly setpoint with Vdc reference list:{}'.format(self.Vdc_ref_list))
            
            elif self.t_reconnect_start == 0.0:
                self.t_reconnect_start = t
                self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='reconnect_start')
    
    def RT_initialize(self,DER_arguments):
        """Initialize VRT and FRT settings."""      
                
        self.VRT_initialize()
        self.FRT_initialize()        
    
    def VRT_initialize(self):
        """Initialize LVRT and HVRT settings."""
        
        #Common ride through variables
        self.RT_logic_t_lock = 0.0
        self.connect_logic_t_lock = 0.0
        self.t_disconnect_start = 0.0
        self.t_reconnect_start = 0.0
        
        #LVRT flags    
        self.LVRT_ENABLE = True
        self.LVRT_TRIP = False
        self.LVRT_MOMENTARY_CESSATION = False
        
        #HVRT flags
        self.HVRT_ENABLE = True
        self.HVRT_TRIP = False
        self.HVRT_MOMENTARY_CESSATION = False
               
        
        
                
        
        self.LVRT_dict = self.RT_config['LVRT']
        self.HVRT_dict = self.RT_config['HVRT']
        
        self.t_disconnect_delay = self.RT_config['VRT_delays']['output_cessation_delay']#(1/120.0)
        self.t_reconnect_delay = self.RT_config['VRT_delays']['output_restore_delay'] 
        self.RESTORE_Vdc = self.RT_config['VRT_delays']['restore_Vdc'] #Restore Vdc to nominal reference by using a ramp
                          
        self.check_VRT_settings()
    
    def LVRT(self,t):
        """Function to implement LVRT ridethrough and trip logic."""
        
        #Select RMS voltage source
        Vrms_measured = self.Vrms   #Select PCC - LV side voltage             
        
        if t > self.t_stable: #Go through logic only after a short time delay
            for LVRT_key,LVRT_values in self.LVRT_dict.items():
                zone_name = 'LV'+str(LVRT_key)
                
                V_threshold = LVRT_values['V_threshold']*self.Vrms_ref #Convert % thresholds into p.u.
                t_threshold = LVRT_values['t_threshold']
                
                if Vrms_measured < V_threshold and not LVRT_values['threshold_breach']: #Check if voltage below threshold
                    if LVRT_values['t_start'] == 0.0: #Start timer if voltage goes below threshold
                        LVRT_values['t_start']  = t
                        self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],event_name='zone_entered')
                        
                        if LVRT_values['mode'] == 'momentary_cessation': #Go into momentary cessation
                            self.LVRT_MOMENTARY_CESSATION = True
                            self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],event_name='momentary_cessation')
                        
                    elif t-LVRT_values['t_start'] <= t_threshold: #Remain in LV zone and monitor
                        self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],event_name='zone_continue')
                      
                    elif t-LVRT_values['t_start'] >= t_threshold: #Trip DER if timer exceeds threshold
                        self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],event_name='trip')
                        LVRT_values['threshold_breach'] = True                           
                        self.LVRT_TRIP = True
                        LVRT_values['t_start'] = 0.0
                    
                elif  Vrms_measured > V_threshold: #Check if voltage above threshold
                    if LVRT_values['t_start'] > 0.0: #Reset timer if voltage goes above  threshold
                        self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],event_name='zone_reset')
                        LVRT_values['t_start']  = 0.0 
                        self.LVRT_MOMENTARY_CESSATION = False #Reset momentary cessation flags
                    else: #Do nothing
                        pass

    def HVRT(self,t):
        """Function to implement HVRT ridethrough and trip logic."""
        
        #Select RMS voltage source
        Vrms_measured = self.Vrms   #Select PCC - LV side voltage     
        
        if t > self.t_stable: #Go through logic only after a short time delay
            for HVRT_key,HVRT_values in self.HVRT_dict.items():
                zone_name = 'HV'+str(HVRT_key)
                V_threshold = HVRT_values['V_threshold']*self.Vrms_ref #Convert % thresholds into p.u.
                t_threshold = HVRT_values['t_threshold']
                
                if Vrms_measured > V_threshold and not HVRT_values['threshold_breach']: #Check if voltage above threshold
                    if HVRT_values['t_start'] == 0.0: #Start timer if voltage goes above threshold
                        HVRT_values['t_start']  = t
                        self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],event_name='zone_entered')
                        
                        if HVRT_values['mode'] == 'momentary_cessation': #Go into momentary cessation
                            self.HVRT_MOMENTARY_CESSATION = True
                            self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],event_name='momentary_cessation')
                    
                    elif t-HVRT_values['t_start'] <= t_threshold: #Remain in LV zone and monitor
                        self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],event_name='zone_continue')
                      
                    elif t-HVRT_values['t_start'] >= t_threshold: #Trip DER if timer exceeds threshold
                        self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],event_name='trip')
                        HVRT_values['threshold_breach'] = True                           
                        self.HVRT_TRIP = True
                        HVRT_values['t_start'] = 0.0
                    
                elif  Vrms_measured < V_threshold: #Check if voltage below threshold
                    if HVRT_values['t_start'] > 0.0: #Reset timer if voltage goes below threshold
                        self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],event_name='zone_reset')
                        HVRT_values['t_start']  = 0.0 
                        self.HVRT_MOMENTARY_CESSATION = False #Reset momentary cessation flags
                    else: #Do nothing
                        pass
    
    def show_RT_flags(self):
        """Show all RT flags."""
        
        print('DER_MOMENTARY_CESSATION:{},LVRT_MOMENTARY_CESSATION:{},HVRT_MOMENTARY_CESSATION:{}'.format(self.DER_MOMENTARY_CESSATION,self.LVRT_MOMENTARY_CESSATION,self.HVRT_MOMENTARY_CESSATION))
        print('DER_CONNECTED:{},DER_TRIP:{}'.format(self.DER_CONNECTED,self.DER_TRIP))
    
    def check_VRT_settings(self):
        """Sanity check for VRT settings."""
        
        for RT_type in ['LV','HV']:
            V_threshold_list = []
            t_threshold_list = []
            
            if RT_type == 'LV':
                VRT_dict = self.LVRT_dict
            elif RT_type == 'HV':
                VRT_dict = self.HVRT_dict
            
            for VRT_key,VRT_values in VRT_dict.items():
                zone_name = RT_type+str(VRT_key)
                V_threshold_list.append(VRT_values['V_threshold'])
                t_threshold_list.append(VRT_values['t_threshold'])
                voltage_error_text = '{} is not a valid voltage threshold for zone {}!'.format(VRT_values['V_threshold'],zone_name) 
                time_error_text = '{} is not a valid time threshold for zone {}!'.format(VRT_values['t_threshold'],zone_name) 
                
                if RT_type == 'LV':
                    if VRT_values['V_threshold'] >= 1.0:
                        raise ValueError(voltage_error_text) 
                    if not six.PY2: #Since dictionaries don't preserve insertion order in Python 2
                        if not V_threshold_list == sorted(V_threshold_list): #Check if voltage values are in ascending order
                            raise ValueError('LVRT voltage limits - {} are infeasible!'.format(V_threshold_list))

                elif RT_type == 'HV':
                    if VRT_values['V_threshold'] <= 1.0:
                        raise ValueError(voltage_error_text) 
                    if not six.PY2: #Since dictionaries don't preserve insertion order in Python 2
                        if not V_threshold_list == sorted(V_threshold_list,reverse=True): #Check if voltage values are in descending order
                            raise ValueError('HVRT voltage limits - {} are infeasible!'.format(V_threshold_list))
                    
                if VRT_values['V_threshold'] <= 0.0:
                    raise ValueError(voltage_error_text) 
                if VRT_values['V_threshold'] >=self.V_threshold_high_limit:
                    raise ValueError(voltage_error_text) 
                
                if VRT_values['t_threshold'] <= 0.0:
                    raise ValueError(time_error_text) 
                if VRT_values['t_threshold'] >= self.t_threshold_high_limit:
                    raise ValueError(time_error_text) 

                if VRT_values['mode'] not in {'mandatory_operation','momentary_cessation'}:
                    raise ValueError('{} is not a valid mode for zone {}'.format(VRT_values['mode'],zone_name))    
            
            if self.t_disconnect_delay < self.t_disconnect_low_limit:
                raise ValueError('DER output cessation time delay {} s is infeasible!'.format(self.t_disconnect_delay))   
            
            if self.t_reconnect_delay < self.t_reconnect_low_limit:
                raise ValueError('DER output restore time delay {} s after momentary cessation is infeasible!'.format(self.t_reconnect_delay))   
    
    
    def print_VRT_events(self,simulation_time,voltage,zone_name,timer_start=0.0,event_name='',print_inline = True,verbose = False):
        """Print logs for VRT events."""
        
        voltage = (voltage*self.Vbase)/(self.Vrms_ref*self.Vbase)#175

        if event_name == 'zone_entered':
            text_string = '{}:{:.4f}:{} zone entered at {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V)'                           .format(self.name,simulation_time,zone_name,timer_start,voltage,self.Vrms_ref*self.Vbase)
        
        elif event_name == 'zone_reset':
            text_string = '{}:{:.4f}:{} flag reset at {:.4f} s after {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V)'.format(self.name,simulation_time,zone_name,simulation_time,simulation_time-timer_start,voltage,self.Vrms_ref*self.Vbase)
        
        elif event_name == 'zone_continue' and verbose:
            text_string = '{}:{:.4f}:{} zone entered at:{:.4f} s and continuing for {:.4f} s'\
                            .format(self.name,simulation_time,zone_name,timer_start,simulation_time-timer_start)
        
        elif event_name == 'momentary_cessation':
            text_string = '{}:{:.4f}:{} zone - momentary cessation at {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V)'                           .format(self.name,simulation_time,zone_name,timer_start,voltage,self.Vrms_ref*self.Vbase)
            six.print_(text_string)
        
        elif event_name == 'trip':
            text_string = '{}:{:.4f}:{} violation at {:.4f}s after {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V) - DER will be tripped'\
                            .format(self.name,simulation_time,zone_name,simulation_time,simulation_time-timer_start,voltage,self.Vrms_ref*self.Vbase)
            six.print_(text_string)
        
        else:
            text_string =''

        self.print_event(text_string,print_inline)           
    
    def print_reconnect_events(self,simulation_time,voltage,frequency,timer_start=0.0,event_name='',print_inline = True,verbose = False):
        """Print logs for VRT events."""
        
        voltage = (voltage*self.Vbase)/(self.Vrms_ref*self.Vbase)#175
        
        if event_name == 'reconnect_start':
            text_string = '{}:{:.4f}:Reconnect timer started at {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V)'\
                            .format(self.name,simulation_time,timer_start,voltage,self.Vrms_ref*self.Vbase)

        elif event_name == 'reconnect_reset':
            text_string = '{}:{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                           .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'reconnect_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
                           .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'DER_reconnection':
            text_string = '{}:{time_stamp:.4f}:DER reconnecting after momentary cessation at {time_stamp:.4f}s after {time_elasped:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
            six.print_(text_string)

        elif event_name == 'DER_tripped' and verbose: 
            text_string = '{}:{time_stamp:.4f}:Inverter in tripped condition for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'.format(self.name,time_stamp=simulation_time,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        else:
            text_string =''
        
        self.print_event(text_string,print_inline)
    
    def check_anomaly(self):
        """Check if voltage anomaly was detected."""
        
        for LVRT_key,LVRT_values in self.LVRT_dict.items():
            if LVRT_values['t_start'] > 0.0: #Check if any timer started
                print('Voltage anomaly started at {:.4f} due to zone {} with threshold {:.2f} V p.u.'.format(LVRT_values['t_start'],LVRT_key,LVRT_values['V_threshold']))
      
    def update_RT_config_old(self,derConfig):
        """Check whether the config file is good."""        
                   
        for item in templates.RT_config_template.keys():
            if item in derConfig:
                self.RT_config[item] = derConfig[item]
            elif item in self.DER_config:
                self.RT_config[item] = self.DER_config[item]
                self.logger.debug('{}:{} updated with value from config file {}.'.format(self.name,item,self.DER_config[item]))
            elif item in self.default_RT_config:
                self.logger.debug('{}:{} updated with default value {}.'.format(self.name,item,self.default_RT_config[item]))    
                self.RT_config[item] = self.default_RT_config[item]
            else:
                raise KeyError('{}:Ridethrough setting {} could not be found!'.format(self.name,item))
    
    def update_RT_config(self,DER_config,DER_arguments,config_dict):
        """Check whether the config file is good."""             
                     
        for RT in list(templates.VRT_config_template.keys()) +  list(templates.FRT_config_template.keys()):
            if RT in DER_arguments['derConfig']:
                self.RT_config[RT] = DER_arguments['derConfig'][RT]
                self.logger.debug('{}:{} updated with settings from derConfig.'.format(self.name,RT))
            elif RT in self.DER_config:              
                if 'config_id' in self.DER_config[RT]:
                    self.RT_config[RT] = config_dict[self.DER_config[RT]['config_id']]['config']
                    self.logger.debug('{}:{} updated with settings from config id {}.'.format(self.name,RT,self.DER_config[RT]['config_id']))
                elif 'config' in self.DER_config[RT]:
                    self.RT_config[RT] = self.DER_config[RT]['config'] 
                    self.logger.debug('{}:{} updated with settings from DER config file.'.format(self.name,RT))
            else:
                if RT in list(templates.VRT_config_template.keys()):
                    self.RT_config[RT] = templates.VRT_config_template[RT]['config']
                if RT in list(templates.FRT_config_template.keys()):
                    self.RT_config[RT] = templates.FRT_config_template[RT]['config']
                self.logger.debug('{}:{} updated with settings from template.'.format(self.name,RT))
   
        for RT in ['LVRT','HVRT']:
            for RT_level,RT_config in self.RT_config[RT].items():
               
                assert isinstance(RT_config,dict), 'Ridethrough level must be specified as a dictionary!'
                if 't_start' not in RT_config.keys():
                    RT_config.update({'t_start':0.0})
                if 'threshold_breach' not in RT_config.keys():
                    RT_config.update({'threshold_breach':False})    
    
    def check_RT_config(self):
        """Check whether the config file is good."""        
        
        for RT in list(templates.VRT_config_template.keys()):
            if RT != 'VRT_delays':
                for RT_setting in templates.VRT_config_template[RT]['config']['0']:                    
                    for level,setting in self.RT_config[RT].items():                        
                        if RT_setting not in setting:
                            raise ValueError('{} not found for level {} in {}!'.format(RT_setting,level,RT))
            
        for RT in list(templates.FRT_config_template.keys()):
            if RT != 'FRT_delays':
                for RT_setting in templates.FRT_config_template[RT]['config']['1']:
                    for level,setting in self.RT_config[RT].items():                        
                        if RT_setting not in setting:
                            raise ValueError('{} not found for level {} in {}!'.format(RT_setting,level,RT))        
    
    def show_RT_settings(self,settings_type='LVRT',PER_UNIT=True):
        """Method to show LVRT settings."""
        
        if settings_type not in {'LVRT','HVRT','LFRT'}:#
            raise ValueError('Unknown quantity: ' + str(settings_type))
        
        if PER_UNIT:
            V_multiplier = (1/self.Vrms_ref)
        else:
            V_multiplier = self.Vbase
        
        print('\n______{} - {}_____'.format(self.name,settings_type))
        
        if settings_type ==  'LVRT':
            print('______Flags______')
            print('LVRT_ENABLE:{}\nLVRT_TRIP:{},LVRT_MOMENTARY_CESSATION:{} '.format(self.LVRT_ENABLE,self.LVRT_TRIP,self.LVRT_MOMENTARY_CESSATION))    
            print('______Thresholds______')
            print('Vrms_ref:{:.2f} V'.format(self.Vrms_ref*self.Vbase))
            for LVRT_key,LVRT_values in self.LVRT_dict.items():
                print('Zone:{},Vthreshold:{:.2f},tthreshold:{:.2f},mode:{}'.format(LVRT_key,LVRT_values['V_threshold'],LVRT_values['t_threshold'],LVRT_values['mode']))
            
        if settings_type ==  'HVRT':
            print('______Flags______')
            print('HVRT_ENABLE:{}\nHVRT_TRIP:{},HVRT_MOMENTARY_CESSATION:{}'.format(self.HVRT_ENABLE,self.HVRT_TRIP,self.HVRT_MOMENTARY_CESSATION))
            print('______Thresholds______')
            print('Vrms_ref:{:.2f} V'.format(self.Vrms_ref*self.Vbase))
            for HVRT_key,HVRT_values in self.HVRT_dict.items():
                print('Zone:{},Vthreshold:{:.2f},tthreshold:{:.2f},mode:{}'.format(HVRT_key,HVRT_values['V_threshold'],HVRT_values['t_threshold'],HVRT_values['mode']))
            
        
        if settings_type ==  'LFRT':
            print('f_ref:{:.2f} Hz\nF_LF1:{:.2f} Hz\nF_LF2:{:.2f} Hz'.format(self.f_ref,self.LFRT_dict['1']['F_LF'],self.LFRT_dict['2']['F_LF']))
            print('t_LF1:{:.2f} s\nt_LF2:{:.2f} s'.format(self.LFRT_dict['1']['t_LF_limit'],self.LFRT_dict['2']['t_LF_limit']))
            print('______Flags______')
            print('LFRT_ENABLE:{}\nLFRT_TRIP:{} '.format(self.LFRT_ENABLE,self.LFRT_TRIP))        
        
        if settings_type ==  'LFRT' or settings_type ==  'HFRT':
            print('FRT_INSTANTANEOUS_TRIP:{}'.format(self.FRT_INSTANTANEOUS_TRIP))        
        print('OUTPUT_CESSATION_DELAY:{},OUTPUT_RESTORE_DELAY:{}'.format(self.t_disconnect_delay,self.t_reconnect_delay))
     
    def FRT_initialize(self):
        """Initialize LFRT and HFRT settings."""
        
        #LVRT flags    
        self.LFRT_ENABLE = False
        self.LFRT_TRIP = False        
        
        #HVRT flags
        self.HFRT_ENABLE = False
        self.HFRT_TRIP = False
        
        self.del_f =0.02

        self.LFRT_dict = self.RT_config['LFRT']
        self.HFRT_dict = self.RT_config['HFRT']
        self.FRT_INSTANTANEOUS_TRIP = self.RT_config['FRT_delays']['FRT_INSTANTANEOUS_TRIP'] #Disconnects PV-DER within one cycle for frequency anomaly
    
        self.check_LFRT_settings()    
        self.check_HFRT_settings()
    
    def check_LFRT_settings(self):
        """Sanity check for LFRT settings."""        
        
        if not self.LFRT_dict['2']['F_LF'] > self.LFRT_dict['1']['F_LF']:
            
            raise ValueError('LFRT frequency settings - F_LF1:{:.2f},F_LF2:{:.2f} are infeasible!'.format(self.LFRT_dict['1']['F_LF'],self.LFRT_dict['2']['F_LF']))
            
    def check_HFRT_settings(self):
        """Sanity check for HFRT settings."""
        
        if not self.HFRT_dict['2']['F_HF'] > self.HFRT_dict['1']['F_HF']:
            
            raise ValueError('HFRT frequency settings - F_HF1:{:.2f},F_HF2:{:.2f} are infeasible!'.format(self.HFRT_dict['1']['F_HF'],self.HFRT_dict['2']['F_HF']))
    
    def FRT(self,t):
        """Frequency ride through and trip logic. """
        
        #Select frequency source
        fgrid = self.we/(2.0*math.pi)  #Use grid frequency as estimated by PLL
        
        if t > self.t_stable: #Go through logic only after a short time delay
            
            if not self.LFRT_TRIP: #Logic before tripping/momentary cessation
                
                if any(self.LFRT_dict[key]['t_LFstart']==0.0 for key in self.LFRT_dict.keys()):
                    
                    for LFRT_key,LFRT_values in self.LFRT_dict.items():
                        if LFRT_values['t_LFstart'] == 0.0 and fgrid < LFRT_values['F_LF']:
                            LFRT_values['t_LFstart']  = t
                            self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='LF'+LFRT_key+'_start')
                            #print(fgrid, LFRT_values['F_LF'],LFRT_values['t_LFstart'])

                if any(self.LFRT_dict[key]['t_LFstart'] > 0.0 for key in self.LFRT_dict.keys()):    
  
                    for LFRT_key,LFRT_values in self.LFRT_dict.items():
                        if fgrid > LFRT_values['F_LF']  + self.del_f and LFRT_values['t_LFstart'] > 0.0: #Reset timer if freqency goes above
                           self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='LF'+LFRT_key+'_reset')
                           LFRT_values['t_LFstart']  = 0.0

                    for LFRT_key,LFRT_values in self.LFRT_dict.items():
                        if fgrid <= LFRT_values['F_LF'] and t-LFRT_values['t_LFstart'] >= LFRT_values['t_LF_limit'] and  LFRT_values['t_LFstart'] > 0.0: #Set HVRT_TRIP flag if timer exeeds limit
                            self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='inverter_trip_LF'+LFRT_key)
                            self.LFRT_TRIP = True
                            LFRT_values['t_LFstart'] = 0.0
                            self.t_LF_reconnect = 0.0

                        elif fgrid <= LFRT_values['F_LF'] and t-LFRT_values['t_LFstart'] < LFRT_values['t_LF_limit'] and  LFRT_values['t_LFstart'] > 0.0: #Remain in LF zone and monitor
                            self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='LF'+LFRT_key+'_zone')

            elif self.LFRT_TRIP: #Logic after tripping/momentary cessation
                self.DER_reconnect_logic(t)
 
        else:
            self.LFRT_TRIP = False
            for LFRT_key,LFRT_values in self.LFRT_dict.items():
                 LFRT_values['t_LFstart'] = 0.0
    
    @property
    def FRT_INSTANTANEOUS_TRIP(self):
        return self.__FRT_INSTANTANEOUS_TRIP
    
    @FRT_INSTANTANEOUS_TRIP.setter
    def FRT_INSTANTANEOUS_TRIP(self,FRT_INSTANTANEOUS_TRIP):
        
        self.__FRT_INSTANTANEOUS_TRIP = FRT_INSTANTANEOUS_TRIP
        
        if FRT_INSTANTANEOUS_TRIP:            
            self.LFRT_dict['1']['t_LF_limit'] = self.LFRT_dict['2']['t_LF_limit']  = 30*(1/60) #Disconnect within n cycles                     
            self.HFRT_dict['1']['t_LF_limit'] = self.HFRT_dict['2']['t_LF_limit']  = 30*(1/60) #Disconnect within n cycles  
            
        else:
            self.LFRT_dict['1']['t_LF_limit'] = self.RT_config['LFRT']['1']['t_LF_limit'] 
            self.LFRT_dict['2']['t_LF_limit'] = self.RT_config['LFRT']['2']['t_LF_limit'] 
            self.HFRT_dict['1']['t_HF_limit']  = self.RT_config['HFRT']['1']['t_HF_limit'] 
            self.HFRT_dict['2']['t_HF_limit']  = self.RT_config['HFRT']['2']['t_HF_limit'] 
                    
        return self.__FRT_INSTANTANEOUS_TRIP          
   
    
    def print_LFRT_events(self,simulation_time,frequency,timer_start=0.0,event_name='',print_inline = False,verbose = False):
        """Print LFRT events."""    
        #print(event_name)
        if event_name == 'LF1_start':
            text_string = '{time_stamp:.4f}:LF1 zone entered at {timer_start:.4f}s for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)

        elif event_name == 'LF2_start':
            text_string = '{time_stamp:.4f}:LF2 zone entered at {timer_start:.4f}s for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)

        elif event_name == 'LF3_start':
            text_string = '{time_stamp:.4f}:LF3 zone entered at {timer_start:.4f}s for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)

        elif event_name == 'LF1_reset':
            text_string = '{time_stamp:.4f}:LF1 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LF1 zone for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)

        elif event_name == 'LF2_reset':
            text_string = '{time_stamp:.4f}:LF2 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LF2 zone for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)

        elif event_name == 'LF3_reset':
            text_string = '{time_stamp:.4f}:LF3 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LF3 zone for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)

        elif event_name == 'LF1_zone' and verbose:
            text_string = '{time_stamp:.4f}:LF1 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'LF2_zone' and verbose:
            text_string = '{time_stamp:.4f}:LF2 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'LF3_zone' and verbose:
            text_string = '{time_stamp:.4f}:LF3 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)


        elif event_name == 'inverter_trip_LF1':
            text_string = '{time_stamp:.4f}:LF1 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
                            .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
            six.print_(text_string)

        elif event_name == 'inverter_trip_LF2':
            text_string = '{time_stamp:.4f}:LF2 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
                            .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)    
            six.print_(text_string)

        elif event_name == 'inverter_trip_LF3':
            text_string = '{time_stamp:.4f}:LF3 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
                            .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)    
            six.print_(text_string)

        elif event_name == 'reconnect_start':
            text_string = '{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)

        elif event_name == 'reconnect_reset':
            text_string = '{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {frequency:.3f} Hz'\
                           .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)

        elif event_name == 'reconnect_zone' and verbose:
            text_string = '{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
                           .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'inverter_reconnection':
            text_string = '{time_stamp:.4f}:Inverter reconnecting after LF trip at {time_stamp:.4f}s after {time_elasped:.4f}s for {frequency:.3f} Hz'\
                            .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
            six.print_(text_string)

        elif event_name == 'inverter_tripped' and verbose: 
            text_string = '{time_stamp:.4f}:Inverter in tripped condition for {frequency:.3f} Hz'.format(time_stamp=simulation_time,frequency=frequency)

        else:
            text_string =''    

        self.print_event(text_string,print_inline)    

   
     
    def print_event(self,text_string,print_inline):
        """Print information about ride through events."""
        
        if text_string != '':
            if print_inline:  #Print log in notebook window
                self.logger.info(text_string)

            else: #Print in console window
                utility_functions.print_to_terminal(text_string)
        else:
            pass
    