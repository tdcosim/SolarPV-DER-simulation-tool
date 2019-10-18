"""Code for features inside PV-DER model instances."""

from __future__ import division
import six

import math

from pvder import utility_functions
from pvder import config

class PVDER_SmartFeatures():
    """Class for describing smart inverter inverter features of PV-DER."""
    
    #Flags
    VOLT_VAR_ENABLE = False
    VOLT_VAR_FLAG = False
    
    VOLT_WATT_ENABLE = False
    VOLT_WATT_FLAG = False
    
    #LVRT_ENABLE = False
    LFRT_ENABLE = False
    
    #LFRT variables
    LFRT_DEBUG = False
    t_LF1start = 0.0
    t_LF2start = 0.0
    t_LF3start = 0.0
    t_LF_reconnect = 0.0
    LFRT_TRIP = False
    LFRT_RECONNECT = False    
    
    config_default = {'V_LV0':0.50,'V_LV1':0.70,'V_LV2':0.88,
                   't_LV0_limit':0.1,'t_LV1_limit':1.0,'t_LV2_limit':2.0,
                   'V_HV1':1.06,'V_HV2':1.12,
                   't_HV1_limit':2.0,'t_HV2_limit':1/60.0,
                   'OUTPUT_RESTORE_DELAY':0.5,
                   'VRT_INSTANTANEOUS_TRIP':False,'VRT_MOMENTARY_CESSATION':True}
    
   
    
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
        
        #LVRT logic        
        if self.LVRT_ENABLE:
            self.LVRT(t)
            
        #HVRT logic
        if self.HVRT_ENABLE:
            self.HVRT(t)
        
        #LFRT trip logic
        if self.LFRT_ENABLE:
            self.FRT(t)
            
    def check_and_trip(self):
        """Check whether any trip flags are true and trip DER."""
        
        #LVRT trip logic
        if self.LVRT_TRIP and not self.LVRT_RECONNECT:
            self.PV_DER_disconnect()
        
        #HVRT trip logic
        if self.HVRT_TRIP and not self.HVRT_RECONNECT:
            self.PV_DER_disconnect()        
        
        #LFRT trip logic
        if self.LFRT_TRIP and not self.LFRT_RECONNECT:
            self.PV_DER_disconnect()         
    
    def VRT_initialize(self):
        """Initialize LVRT settings."""
        
        #LVRT variables
        self.t_LV0start = 0.0
        self.t_LV1start = 0.0
        self.t_LV2start = 0.0
        
        #LVRT flags    
        self.LVRT_ENABLE = True
        self.LVRT_TRIP = False
        self.LVRT_RECONNECT = False
        
        #HVRT flags
        self.HVRT_ENABLE = True
        self.HVRT_TRIP = False
        self.HVRT_RECONNECT = False        
        
        #Common ride through variables
        self.t_reconnect = 0.0
        self.Vreconnect_LV = config.DEFAULT_Vreconnect_LV*self.Vrms_ref
        self.Vreconnect_HV = config.DEFAULT_Vreconnect_HV*self.Vrms_ref 
        
        if self.pvderConfig is None:
            self.pvderConfig = {}
        
        self.update_config() #Checks and updates pvderConfig if any entries are missing
        
        #IEEE 1547-2018 standards        
        #V1 to V2 - zone 2,V1 < - zone 1 
        self.V_LV0 = self.pvderConfig['V_LV0']*self.Vrms_ref  #From IEEE 1557-2018 Category III (Table 16)   
        self.V_LV1 = self.pvderConfig['V_LV1']*self.Vrms_ref  #From IEEE 1557-2018 Category III (Table 16)   
        self.V_LV2 = self.pvderConfig['V_LV2']*self.Vrms_ref  #From IEEE 1557-2018  Category III (Table 16)
        self.t_LV0_limit = self.pvderConfig['t_LV0_limit'] #Time limit for LV0 zone from IEEE 1557-2018 Category III  (Table 16, page 48)
        self.t_LV1_limit = self.pvderConfig['t_LV1_limit'] #Time limit for LV1 zone from IEEE 1557-2018 Category III  (Table 16, page 48)
        self.t_LV2_limit = self.pvderConfig['t_LV2_limit'] #Time limit for LV2 zone from IEEE 1557-2018 Category III (Table 16, page 48) 
        
        self.HVRT_dict = {'1':{'V_HV':self.pvderConfig['V_HV1']*self.Vrms_ref,
                               't_HV_limit':self.pvderConfig['t_HV1_limit'],
                               't_HVstart':0.0
                              },
                          '2':{'V_HV':self.pvderConfig['V_HV2']*self.Vrms_ref,
                               't_HV_limit':self.pvderConfig['t_HV2_limit'],
                               't_HVstart':0.0
                              }
                         }        
       
        """
        self.LVRT_dict = {'1':{'V_LV':self.V_HV,
                               't_LV_limit':self.t_HV_limit,
                               't_LVstart':self.t_HVstart
                              }
                          '2':{'V_LV':self.V_HV,
                               't_LV_limit':self.t_HV_limit,
                               't_LVstart':self.t_HVstart
                              }
                         }
        """
        
        self.VRT_INSTANTANEOUS_TRIP = self.pvderConfig['VRT_INSTANTANEOUS_TRIP'] #Disconnects PV-DER within one cycle for voltage anomaly
        self.VRT_MOMENTARY_CESSATION = self.pvderConfig['VRT_MOMENTARY_CESSATION'] #Reconnect PV-DER after voltage anomaly
        self.check_LVRT_settings()
        
    def update_config(self):
        """Check whether the config file is good."""        
                          
        for item in self.config_default.keys():
            if item not in self.pvderConfig:
                self.logger.info('{}:{} was not in pvderconfig, updating with default value {}.'.format(self.name,item,self.config_default[item]))    
                self.pvderConfig[item] = self.config_default[item]

    def check_LVRT_settings(self):
        """Sanity check for LVRT settings."""
        
        V_LV0_low_limit = 0.45
        V_LV0_high_limit = 0.6
        V_LV1_low_limit = 0.65
        V_LV1_high_limit = 0.75
        V_LV2_low_limit = 0.8
        V_LV2_high_limit = 0.92      
        
        t_LV0_low_limit = 1/120.0
        t_LV0_high_limit = 1.0
        t_LV1_low_limit = 1.0
        t_LV1_high_limit = 10.0
        t_LV2_low_limit = 2.0
        t_LV2_high_limit = 20.0      
        
        #if (self.V_LV0 > self.V_LV1 or self.V_LV1 > self.V_LV2) :
        if not self.V_LV2 > self.V_LV1 > self.V_LV0:
            
            raise ValueError('LVRT voltage limits - V_LV0:{:.2f},V_LV1:{:.2f},V_LV2:{:.2f} are infeasible!'.format(self.V_LV0,self.V_LV1,self.V_LV2))
        
        assert (self.V_LV0 >= V_LV0_low_limit*self.Vrms_ref and self.V_LV0 <= V_LV0_high_limit*self.Vrms_ref), 'V_LV0 should be between {} and {}'.format(V_LV0_low_limit,V_LV0_high_limit)
        assert (self.V_LV1 >= V_LV1_low_limit*self.Vrms_ref and self.V_LV1 <= V_LV1_high_limit*self.Vrms_ref), 'V_LV1 should be between {} and {}'.format(V_LV1_low_limit,V_LV1_high_limit)
        assert (self.V_LV2 >= V_LV2_low_limit*self.Vrms_ref and self.V_LV2 <= V_LV2_high_limit*self.Vrms_ref), 'V_LV2 should be between {} and {}'.format(V_LV2_low_limit,V_LV2_high_limit)        
        
        if (self.t_LV0_limit > self.t_LV1_limit or self.t_LV1_limit > self.t_LV2_limit) and  not self.VRT_INSTANTANEOUS_TRIP:
            raise ValueError('LVRT ridethrough times - t_LV0:{:.2f},t_LV1:{:.2f},t_LV2:{:.2f} are infeasible!'.format(self.t_LV0_limit,self.t_LV1_limit,self.t_LV2_limit))
        
        assert self.t_LV0_limit >= t_LV0_low_limit and self.t_LV0_limit <= t_LV0_high_limit, 't_LV0 should be between {} and {}'.format(t_LV0_low_limit,t_LV0_high_limit)
        assert self.t_LV1_limit >= t_LV1_low_limit and self.t_LV1_limit <= t_LV1_high_limit, 't_LV1 should be between {} and {}'.format(t_LV1_low_limit,t_LV1_high_limit)
        assert self.t_LV2_limit >= t_LV2_low_limit and self.t_LV2_limit <= t_LV2_high_limit, 't_LV2 should be between {} and {}'.format(t_LV2_low_limit,t_LV2_high_limit)
        
        if self.t_reconnect_delay < 0.4:
            raise ValueError('Reconnect time delay after momentary cessation {} is infeasible!'.format(self.t_reconnect_delay))
            
    def check_HVRT_settings(self):
        """Sanity check for HVRT settings."""
        
        V_HV1_low_limit = 1.06
        V_HV1_high_limit = 1.1
        V_HV2_low_limit = 1.12
        V_HV2_high_limit = 1.2      
        
        t_HV1_low_limit = 2.0
        t_HV1_high_limit = 12.0
        t_HV2_low_limit = 1/120.0
        t_HV2_high_limit = 1/60.0     
        
        if self.HVRT_dict['1']['V_HV'] > self.HVRT_dict['2']['V_HV']:
            
            raise ValueError('HVRT voltage limits - V_HV1:{:.2f},V_HV2:{:.2f} are infeasible!'.format(self.HVRT_dict['1']['V_HV'],self.HVRT_dict['2']['V_HV'] ))
        
        assert (self.HVRT_dict['1']['V_HV'] >= V_HV1_low_limit*self.Vrms_ref and self.HVRT_dict['1']['V_HV'] <= V_HV1_high_limit*self.Vrms_ref), 'V_HV1 should be between {} and {}'.format(V_HV1_low_limit,V_HV1_high_limit)
        assert (self.HVRT_dict['2']['V_HV'] >= V_HV2_low_limit*self.Vrms_ref and self.HVRT_dict['2']['V_HV'] <= V_HV2_high_limit*self.Vrms_ref), 'V_HV2 should be between {} and {}'.format(V_HV2_low_limit,V_HV2_high_limit)        
        
        if self.HVRT_dict['2']['t_HV_limit'] > self.HVRT_dict['1']['t_HV_limit'] and  not self.VRT_INSTANTANEOUS_TRIP:
            raise ValueError('HVRT ridethrough times - t_HV1:{:.2f},t_HV2:{:.2f} are infeasible!'.format(self.HVRT_dict['1']['t_HV_limit'],self.HVRT_dict['2']['t_HV_limit']))
        
        assert self.HVRT_dict['1']['t_HV_limit'] >= t_HV1_low_limit and self.HVRT_dict['1']['t_HV_limit']  <= t_HV1_high_limit, 't_HV1 should be between {} and {}'.format(t_LV1_low_limit,t_LV1_high_limit)
        assert self.HVRT_dict['2']['t_HV_limit'] >= t_HV2_low_limit and self.HVRT_dict['2']['t_HV_limit']  <= t_HV2_high_limit, 't_HV2 should be between {} and {}'.format(t_LV2_low_limit,t_LV2_high_limit)
    
    def show_RT_settings(self,settings_type='LVRT',PER_UNIT=True):
        """Method to show LVRT settings."""
        
        if settings_type not in {'LVRT','HVRT'}:#,'LFRT'
            raise ValueError('Unknown quantity: ' + str(settings_type))
        
        print('\n______{} - {}_____'.format(self.name,settings_type))
        
        if settings_type ==  'LVRT':
            print('Vrms_ref:{:.2f} V\nV_LV0:{:.2f} V\nV_LV1:{:.2f} V\nV_LV2:{:.2f} V'.format(self.Vrms_ref,self.V_LV0,self.V_LV1,self.V_LV2))
            print('t_LV0:{:.2f} s\nt_LV1:{:.2f} s\nt_LV2:{:.2f} s'.format(self.t_LV0_limit,self.t_LV1_limit,self.t_LV2_limit))
            print('______Flags______')
            print('LVRT_ENABLE:{}\nLVRT_TRIP:{} '.format(self.LVRT_ENABLE,self.LVRT_TRIP))    
            
        if settings_type ==  'HVRT':
            print('Vrms_ref:{:.2f} V\nV_HV1:{:.2f} V\nV_HV2:{:.2f} V'.format(self.Vrms_ref,self.HVRT_dict['1']['V_HV'],self.HVRT_dict['2']['V_HV']))
            print('t_HV1:{:.2f} s\nt_HV2:{:.2f} s'.format(self.HVRT_dict['1']['t_HV_limit'],self.HVRT_dict['2']['t_HV_limit']))
            print('______Flags______')
            print('HVRT_ENABLE:{}\nHVRT_TRIP:{} '.format(self.HVRT_ENABLE,self.HVRT_TRIP))
        
        print('VRT_INSTANTANEOUS_TRIP:{},VRT_MOMENTARY_CESSATION:{},OUTPUT_RESTORE_DELAY:{}'.format(self.VRT_INSTANTANEOUS_TRIP,self.VRT_MOMENTARY_CESSATION,self.t_reconnect_delay))
    
    def LVRT(self,t):
        """Function to implement LVRT trip and reconnect logic."""
        
        #Select RMS voltage source
        _Vrms_measured = self.Vrms   #Select PCC - LV side voltage     

        if self.LVRT_ENABLE and t > self.t_stable and not self.LVRT_TRIP:

            if self.t_LV0start == 0.0 or self.t_LV1start == 0.0 or self.t_LV2start == 0.0:
                if self.t_LV0start == 0.0 and _Vrms_measured < self.V_LV0:
                    self.t_LV0start = t
                    self.print_LVRT_events(t,_Vrms_measured,self.t_LV0start,event_name='LV0_start')

                if self.t_LV1start == 0.0 and _Vrms_measured < self.V_LV1:
                    self.t_LV1start = t
                    self.print_LVRT_events(t,_Vrms_measured,self.t_LV1start,event_name='LV1_start')

                if self.t_LV2start == 0.0 and self.Vrms < self.V_LV2:
                    self.t_LV2start = t
                    self.print_LVRT_events(t,_Vrms_measured,self.t_LV2start,event_name='LV2_start')

            if self.t_LV0start > 0.0 or self.t_LV1start > 0.0 or self.t_LV2start > 0.0:

                if  _Vrms_measured >= self.V_LV0 or _Vrms_measured >= self.V_LV1  or  _Vrms_measured >= self.V_LV2:
                    
                    if  _Vrms_measured >= self.V_LV0 and self.t_LV0start > 0.0:
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV0start,event_name='LV0_reset')
                        self.t_LV0start = 0.0
                    if  _Vrms_measured >= self.V_LV1 and self.t_LV1start > 0.0:
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV1start,event_name='LV1_reset')
                        self.t_LV1start = 0.0
                    if  _Vrms_measured >= self.V_LV2 and self.t_LV2start > 0.0:
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV2start,event_name='LV2_reset')   
                        self.t_LV2start = 0.0

                if t-self.t_LV0start >= self.t_LV0_limit and self.t_LV0start > 0.0:   #Trip inverter if any limit breached
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV0start,event_name='inverter_trip_LV1')
                        self.LVRT_trip_signals()
                        
                if t-self.t_LV1start >= self.t_LV1_limit and self.t_LV1start > 0.0:   #Trip inverter if any limit breached
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV1start,event_name='inverter_trip_LV1')
                        #six.print_(t-self.t_LV1start,self.t_LV1_limit)  #debug script
                        self.LVRT_trip_signals()
                        
                elif t-self.t_LV2start >= self.t_LV2_limit and self.t_LV2start > 0.0:
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV2start,event_name='inverter_trip_LV2')
                        self.LVRT_trip_signals()
                        
                if _Vrms_measured < self.V_LV0 and self.t_LV0start > 0.0:
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV0start,event_name='LV0_zone') 
                elif _Vrms_measured < self.V_LV1 and self.t_LV1start > 0.0:
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV1start,event_name='LV1_zone') #,verbose = False
                elif _Vrms_measured < self.V_LV2 and self.t_LV2start > 0.0:
                        self.print_LVRT_events(t,_Vrms_measured,self.t_LV2start,event_name='LV2_zone')
                        #six.print_(self.t_LV2start)

        elif self.LVRT_ENABLE and t > self.t_stable and self.LVRT_TRIP:
            
            #Select RMS voltage source
            _Vrms_measured = self.Vrms   #Select PCC - LV side voltage     

            if self.t_reconnect > 0.0:
                if  _Vrms_measured < self.V_LV2:
                    self.print_LVRT_events(t, _Vrms_measured,self.t_reconnect,event_name='reconnect_reset')
                    self.t_reconnect = 0.0
                elif  _Vrms_measured >= self.V_LV2 and t-self.t_reconnect >= self.t_reconnect_delay:
                    self.print_LVRT_events(t, _Vrms_measured,self.t_reconnect,event_name='inverter_reconnection')
                    self.LVRT_TRIP = False             #Reset trip flag
                    self.LVRT_RECONNECT = True        #Set reconnect flag
                    self.t_reconnect = 0.0    #Reset timer
                elif _Vrms_measured >= self.V_LV2 and t-self.t_reconnect < self.t_reconnect_delay:
                    self.print_LVRT_events(t, _Vrms_measured,self.t_reconnect,event_name='reconnect_zone')

            elif self.t_reconnect == 0.0:
                if  _Vrms_measured >= self.V_LV2:
                    self.t_reconnect = t
                    self.print_LVRT_events(t, _Vrms_measured,self.t_reconnect,event_name='reconnect_start')
                else:
                    self.print_LVRT_events(t, _Vrms_measured,event_name='inverter_tripped')

        else:
            self.t_LV0start = 0.0
            self.t_LV1start = 0.0
            self.t_LV2start = 0.0
            self.LVRT_TRIP = False
        
    def LVRT_trip_signals(self):
        """ Trip Signals. """
        
        self.LVRT_TRIP = True
        self.t_LV0start = 0.0
        self.t_LV1start = 0.0
        self.t_LV2start = 0.0
        self.t_reconnect = 0.0
    
    def HVRT(self,t):
        """Function to implement HVRT trip and reconnect logic."""
        
        #Select RMS voltage source
        _Vrms_measured = self.Vrms   #Select PCC - LV side voltage     
        
        if t > self.t_stable: #Go through logic only after a short time delay
            
            if not self.HVRT_TRIP: #Logic before tripping/momentary cessation
                
                if any(self.HVRT_dict[key]['t_HVstart']==0.0 for key in self.HVRT_dict.keys()):
                    
                    for HVRT_key,HVRT_values in self.HVRT_dict.items():
                        if HVRT_values['t_HVstart'] == 0.0 and _Vrms_measured > HVRT_values['V_HV']:
                            HVRT_values['t_HVstart']  = t
                            self.print_HVRT_events(t,_Vrms_measured,HVRT_values['t_HVstart'],event_name='HV'+HVRT_key+'_start')

                if any(self.HVRT_dict[key]['t_HVstart'] > 0.0 for key in self.HVRT_dict.keys()):    
  
                    for HVRT_key,HVRT_values in self.HVRT_dict.items():
                        if _Vrms_measured < HVRT_values['V_HV'] and HVRT_values['t_HVstart'] > 0.0: #Reset timer if voltage goes below
                           self.print_HVRT_events(t,_Vrms_measured,HVRT_values['t_HVstart'],event_name='HV'+HVRT_key+'_reset')
                           HVRT_values['t_HVstart']  = 0.0

                    for HVRT_key,HVRT_values in self.HVRT_dict.items():
                        if _Vrms_measured >= HVRT_values['V_HV'] and t-HVRT_values['t_HVstart'] >= HVRT_values['t_HV_limit'] and  HVRT_values['t_HVstart'] > 0.0: #Set HVRT_TRIP flag if timer exeeds limit
                            self.print_HVRT_events(t,_Vrms_measured,HVRT_values['t_HVstart'],event_name='inverter_trip_HV'+HVRT_key)
                            self.HVRT_TRIP = True
                            HVRT_values['t_HVstart'] = 0.0
                            self.t_reconnect = 0.0

                        elif _Vrms_measured >= HVRT_values['V_HV'] and t-HVRT_values['t_HVstart'] < HVRT_values['t_HV_limit'] and  HVRT_values['t_HVstart'] > 0.0: #Remain in HV zone and monitor
                            self.print_HVRT_events(t,_Vrms_measured,self.t_LV1start,event_name='HV'+HVRT_key+'_zone')

            elif self.HVRT_TRIP: #Logic after tripping/momentary cessation
                self.DER_reconnect_logic(t)
 
        else:
            self.HVRT_TRIP = False
            for HVRT_key,HVRT_values in self.HVRT_dict.items():
                 HVRT_values['t_HVstart'] = 0.0
    
    def PV_DER_disconnect(self):
        """Function to disconnect PV DER from grid."""
        
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

        if self.t_reconnect > 0.0:
            if  Vrms_measured < self.Vreconnect_LV or Vrms_measured > self.Vreconnect_HV: #Reset reconnect timer
                self.print_reconnect_events(t, Vrms_measured,self.t_reconnect,event_name='reconnect_reset')
                self.t_reconnect = 0.0
            elif Vrms_measured >= self.Vreconnect_LV and Vrms_measured <= self.Vreconnect_HV and t-self.t_reconnect >= self.t_reconnect_delay:  #Reset HVRT_TRIP flag and set LVRT_Reconnect flag
                self.print_reconnect_events(t, Vrms_measured,self.t_reconnect,event_name='inverter_reconnection')
                self.HVRT_TRIP = False             #Reset trip flag
                self.HVRT_RECONNECT = True        #Set reconnect flag
                self.t_reconnect = 0.0    #Reset timer
            elif Vrms_measured >= self.Vreconnect_LV and Vrms_measured <= self.Vreconnect_HV and t-self.t_reconnect < self.t_reconnect_delay:
                self.print_reconnect_events(t, Vrms_measured,self.t_reconnect,event_name='reconnect_zone')

        elif self.t_reconnect == 0.0:
            if  Vrms_measured >= self.Vreconnect_LV and Vrms_measured <= self.Vreconnect_HV: #Start reconnect timer if voltage is nominal
                self.t_reconnect = t
                self.print_reconnect_events(t, Vrms_measured,self.t_reconnect,event_name='reconnect_start')
            else:
                self.print_reconnect_events(t, Vrms_measured,event_name='inverter_tripped')    
    
    def FRT(self,t):
        """Frequency ride through and trip logic. """
        
        t_reconnect_limit = 3.0   #Time lag before reconnecting

        #IEEE 1547-2018 standards
        F_LF1 = 57.0  #From IEEE 1557-2018 Category III    
        F_LF2 = 58.8  #From IEEE 1557-2018 Category III    
        t_LF1_limit = 0.16      #Time limit for LF1 zone From IEEE 1557-2018 Category III (Table 19 and figure H10)   
        t_LF2_limit = 299.0      #Time limit for LF2 zone From IEEE 1557-2018 Category III (Table 19 and figure H10)
        
        #Dummy variables to preserve logic
        F_LF3 = 0.0
        t_LF3_limit = 0.0
               
        del_f =0.02
        assert (F_LF1 < F_LF2) and (t_LF1_limit < t_LF2_limit), "Frequency level 2 should be greater than frequency level 1"

        #Use grid frequency as estimated by PLL
        fgrid = self.we/(2.0*math.pi)
        
        if self.LFRT_ENABLE and t >  self.t_stable and not self.LFRT_TRIP:
            
            if self.LFRT_DEBUG:
                text_string = '{time_stamp:.4f} -- fgrid:{f1:.3f} Hz, t_LF1start:{t1:.3f} s, t_LF2start:{t2:.3f} s, t_LF3start:{t3:.3f} s'\
                               .format(time_stamp=t,f1 = fgrid,t1=self.t_LF1start,t2=self.t_LF2start,t3=self.t_LF3start)
                utility_functions.print_to_terminal(text_string)
            
            if self.t_LF1start == 0.0 or self.t_LF2start == 0.0 or self.t_LF3start == 0.0:
            
                if self.t_LF1start == 0.0 and fgrid < F_LF1:
                    self.t_LF1start = t
                    utility_functions.print_LFRT_events(t,fgrid,self.t_LF1start,event_name='LF1_start')
                
                if self.t_LF2start == 0.0 and fgrid < F_LF2:
                    self.t_LF2start = t
                    utility_functions.print_LFRT_events(t,fgrid,self.t_LF2start,event_name='LF2_start')
                
                if self.t_LF3start == 0.0 and fgrid < F_LF3:
                    self.t_LF3start = t
                    utility_functions.print_LFRT_events(t,fgrid,self.t_LF3start,event_name='LF3_start')
                
            if self.t_LF1start > 0.0 or self.t_LF2start > 0.0 or self.t_LF3start > 0.0:
                
                if fgrid >= F_LF1+del_f  or fgrid >= F_LF2+del_f or fgrid >= F_LF3+del_f:
                    if fgrid >= F_LF1+del_f and self.t_LF1start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF1start,event_name='LF1_reset')
                        self.t_LF1start = 0.0
                    if fgrid >= F_LF2+del_f and self.t_LF2start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF2start,event_name='LF2_reset')   
                        self.t_LF2start = 0.0
                    if fgrid >= F_LF3+del_f and self.t_LF3start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF3start,event_name='LF3_reset')   
                        self.t_LF3start = 0.0
                
                if t-self.t_LF1start >= t_LF1_limit and self.t_LF1start > 0.0:   #Trip inverter if any limit breached
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF1start,event_name='inverter_trip_LF1')
                        #six.print_(t-self.t_LF1start)
                        self.LFRT_trip_signals()
                        
                elif t-self.t_LF2start >= t_LF2_limit and self.t_LF2start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF2start,event_name='inverter_trip_LF2')
                        self.LFRT_trip_signals()
                        
                elif t-self.t_LF3start >= t_LF3_limit and self.t_LF3start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF3start,event_name='inverter_trip_LF3')
                        self.LFRT_trip_signals()
                        
                if fgrid < F_LF1 and self.t_LF1start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF1start,event_name='LF1_zone') #,verbose = False
                if fgrid < F_LF2 and self.t_LF2start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF2start,event_name='LF2_zone')
                if fgrid < F_LF3 and self.t_LF3start > 0.0:
                        utility_functions.print_LFRT_events(t,fgrid,self.t_LF3start,event_name='LF3_zone')
                        #six.print_(self.t_LF2start)

        elif self.LFRT_ENABLE and t > self.t_stable and self.LFRT_TRIP:

            if self.t_LF_reconnect > 0.0:
                if fgrid < F_LF3:
                    utility_functions.print_LFRT_events(t,fgrid,self.t_LF_reconnect,event_name='reconnect_reset')
                    self.t_LF_reconnect = 0.0
                elif fgrid >= F_LF3 and t-self.t_LF_reconnect >= t_reconnect_limit:
                    utility_functions.print_LFRT_events(t,fgrid,self.t_LF_reconnect,event_name='inverter_reconnection')
                    self.LFRT_TRIP = False             #Reset trip flag
                    self.LFRT_RECONNECT = True        #Set reconnect flag
                    self.t_LF_reconnect = 0.0    #Reset timer
                elif fgrid >= F_LF3 and t-self.t_LF_reconnect < t_reconnect_limit:
                     utility_functions.print_LFRT_events(t,fgrid,self.t_LF_reconnect,event_name='reconnect_zone')
            elif self.t_LF_reconnect == 0.0:
                if fgrid >= F_LF3:
                    self.t_LF_reconnect = t
                    utility_functions.print_LFRT_events(t,fgrid,self.t_LF_reconnect,event_name='reconnect_start')
                else:
                    utility_functions.print_LFRT_events(t,fgrid,event_name='inverter_tripped')

        else:
            self.t_LF1start = 0.0
            self.t_LF2start = 0.0
            self.t_LF3start = 0.0
            self.LFRT_TRIP = False
    
    def LFRT_trip_signals(self):
        """ Trip Signals. """
        
        self.LFRT_TRIP = True
        self.t_LF1start = 0.0
        self.t_LF2start = 0.0
        self.t_LF3start = 0.0
        self.t_LF_reconnect = 0.0
        
    @property
    def VRT_INSTANTANEOUS_TRIP(self):
        return self.__VRT_INSTANTANEOUS_TRIP
    
    @VRT_INSTANTANEOUS_TRIP.setter
    def VRT_INSTANTANEOUS_TRIP(self,VRT_INSTANTANEOUS_TRIP):
        
        self.__VRT_INSTANTANEOUS_TRIP = VRT_INSTANTANEOUS_TRIP
        
        if VRT_INSTANTANEOUS_TRIP:            
            self.t_LV0_limit = self.t_LV1_limit = self.t_LV2_limit = 1/60 #Disconnect within one cycle
            self.HVRT_dict['1']['t_HV_limit'] = self.HVRT_dict['2']['t_HV_limit']  = 1/60 #Disconnect within one cycle                        
        else:
            self.t_LV0_limit = self.pvderConfig['t_LV0_limit']
            self.t_LV1_limit = self.pvderConfig['t_LV1_limit'] 
            self.t_LV2_limit = self.pvderConfig['t_LV2_limit']
            self.HVRT_dict['1']['t_HV_limit'] = self.pvderConfig['t_HV1_limit']
            self.HVRT_dict['2']['t_HV_limit']  = self.pvderConfig['t_HV2_limit']
                    
        return self.__VRT_INSTANTANEOUS_TRIP
    
    @property
    def VRT_MOMENTARY_CESSATION(self):
        return self.__VRT_MOMENTARY_CESSATION
    
    @VRT_MOMENTARY_CESSATION.setter
    def VRT_MOMENTARY_CESSATION(self,VRT_MOMENTARY_CESSATION):
        
        self.__VRT_MOMENTARY_CESSATION = VRT_MOMENTARY_CESSATION
        
        if VRT_MOMENTARY_CESSATION:
            self.t_reconnect_delay = self.pvderConfig['OUTPUT_RESTORE_DELAY']  #Delay before DER power output resumes after voltage anomaly (only valid if momentary cessation flage is true)
        else:
            self.t_reconnect_delay = 1000.0 # A large number to prevent DER from restoring power output after end of voltage anomaly      
                    
        return self.__VRT_MOMENTARY_CESSATION
    
    def print_LVRT_events(self,simulation_time,voltage,timer_start=0.0,event_name='',print_inline = False,verbose = False):
        """Print logs for LVRT events."""
        
        voltage = (voltage*self.Vbase)/(self.Vrms_ref*self.Vbase)#175

        if event_name == 'LV0_start':
            text_string = '{}:{time_stamp:.4f}:LV0 zone entered at {timer_start:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
        
        elif event_name == 'LV1_start':
            text_string = '{}:{time_stamp:.4f}:LV1 zone entered at {timer_start:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'LV2_start':
            text_string = '{}:{time_stamp:.4f}:LV2 zone entered at {timer_start:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'LV0_reset':
            text_string = '{}:{time_stamp:.4f}:LV0 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LV1 zone for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
            
        elif event_name == 'LV1_reset':
            text_string = '{}:{time_stamp:.4f}:LV1 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LV1 zone for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'LV2_reset':
            text_string = '{}:{time_stamp:.4f}:LV2 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LV2 zone for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'LV0_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:LV0 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'LV1_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:LV1 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'LV2_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:LV2 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'inverter_trip_LV0':
            text_string = '{}:{time_stamp:.4f}:LV0 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V) - Inverter will be tripped'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
            six.print_(text_string)
        
        elif event_name == 'inverter_trip_LV1':
            text_string = '{}:{time_stamp:.4f}:LV1 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V) - Inverter will be tripped'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
            six.print_(text_string)

        elif event_name == 'inverter_trip_LV2':
            text_string = '{}:{time_stamp:.4f}:LV2 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V) - Inverter will be tripped'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)    
            six.print_(text_string)

        elif event_name == 'reconnect_start':
            text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'reconnect_reset':
            text_string = '{}:{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                           .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'reconnect_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
                           .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'inverter_reconnection':
            text_string = '{}:{time_stamp:.4f}:Inverter reconnecting after LV trip at {time_stamp:.4f}s after {time_elasped:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
            six.print_(text_string)

        elif event_name == 'inverter_tripped' and verbose: 
            text_string = '{}:{time_stamp:.4f}:Inverter in tripped condition for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'.format(self.name,time_stamp=simulation_time,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        else:
            text_string =''

        self.print_event(text_string,print_inline)       
        
    def print_HVRT_events(self,simulation_time,voltage,timer_start=0.0,event_name='',print_inline = False,verbose = False):
        """Print logs for HVRT events."""
        
        voltage = (voltage*self.Vbase)/(self.Vrms_ref*self.Vbase)#175

        if event_name == 'HV1_start':
            text_string = '{}:{time_stamp:.4f}:HV1 zone entered at {timer_start:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'HV2_start':
            text_string = '{}:{time_stamp:.4f}:HV2 zone entered at {timer_start:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'HV1_reset':
            text_string = '{}:{time_stamp:.4f}:HV1 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in HV1 zone for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'HV2_reset':
            text_string = '{}:{time_stamp:.4f}:HV2 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in HV2 zone for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'HV1_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:HV1 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'HV2_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:HV2 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'inverter_trip_HV1':
            text_string = '{}:{time_stamp:.4f}:HV1 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V) - Inverter will be tripped'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
            six.print_(text_string)

        elif event_name == 'inverter_trip_HV2':
            text_string = '{}:{time_stamp:.4f}:HV2 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V) - Inverter will be tripped'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)    
            six.print_(text_string)

        else:
            text_string =''
        
        self.print_event(text_string,print_inline)
        
    def print_reconnect_events(self,simulation_time,voltage,timer_start=0.0,event_name='',print_inline = False,verbose = False):
        """Print logs for VRT events."""
        
        voltage = (voltage*self.Vbase)/(self.Vrms_ref*self.Vbase)#175
        
        if event_name == 'reconnect_start':
            text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'reconnect_reset':
            text_string = '{}:{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                           .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        elif event_name == 'reconnect_zone' and verbose:
            text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
                           .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'inverter_reconnection':
            text_string = '{}:{time_stamp:.4f}:Inverter reconnecting after momentary cessation at {time_stamp:.4f}s after {time_elasped:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)
            six.print_(text_string)

        elif event_name == 'inverter_tripped' and verbose: 
            text_string = '{}:{time_stamp:.4f}:Inverter in tripped condition for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'.format(self.name,time_stamp=simulation_time,voltage=voltage,Vref=self.Vrms_ref*self.Vbase)

        else:
            text_string =''
        
        self.print_event(text_string,print_inline)
     
    def print_event(self,text_string,print_inline):
        """Print information about ride through events."""
        
        if text_string != '':
            if print_inline:  #Print in notebook window
                logging.info(text_string)

            else: #Print in console window
                utility_functions.print_to_terminal(text_string)
        else:
            pass
    