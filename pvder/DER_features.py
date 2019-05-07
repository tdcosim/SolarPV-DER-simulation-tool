from __future__ import division

from pvder import utility_functions

class PVDER_SmartFeatures():
    """Class for describing smart inverter inverter features of PV-DER."""
    #Flags
    VOLT_VAR_ENABLE = False
    VOLT_VAR_FLAG = False
    
    VOLT_WATT_ENABLE = False
    VOLT_WATT_FLAG = False
    
    LVRT_ENABLE = False
    LFRT_ENABLE = False
    
    #LVRT variables
    t_LV1start = 0.0
    t_LV2start = 0.0
    t_reconnect = 0.0
    LVRT_TRIP = False
    LVRT_RECONNECT = False
    LVRT_INSTANTANEOUS_TRIP = False
    
    #LFRT variables
    LFRT_DEBUG = False
    t_LF1start = 0.0
    t_LF2start = 0.0
    t_LF3start = 0.0
    t_LF_reconnect = 0.0
    LFRT_TRIP = False
    LFRT_RECONNECT = False
    
    def __init__(self,pvderConfig=None):
        """Creates an instance of `PVDER_SmartFeatures`.
        
        Args:
          pvderConfig: A dictionary containing configuration parameters that may be supplied from an external program.
     
        Raises:
          ValueError: If LVRT settings don't make sense.
        """
        pass
    
    def update_ridethrough_flags(self,t):
        """Check VRT and FRT logic."""
        
        #LVRT trip logic        
        if self.LVRT_ENABLE == True:
            self.LVRT(t)
            if self.LVRT_TRIP == True and self.LVRT_RECONNECT == False:
                self.PV_DER_disconnect()
        
        #LFRT trip logic
        if self.LFRT_ENABLE == True:
            self.FRT(t)
            if self.LFRT_TRIP == True and self.LFRT_RECONNECT == False:
                self.PV_DER_disconnect() 
    
    def LVRT_initialize(self,pvderConfig=None):
        """Function to initialize LVRT settings."""
        if pvderConfig is None:
            pvderConfig = {}
            pvderConfig['scaling_factor']=10
            pvderConfig['V_LV1']=0.70
            pvderConfig['V_LV2']=0.88
            pvderConfig['t_LV1_limit']=10.0
            pvderConfig['t_LV2_limit']=20.0
            pvderConfig['LVRT_INSTANTANEOUS_TRIP']=False

        #IEEE 1547-2018 standards
        scaling_factor = pvderConfig['scaling_factor']
        self.V_LV0 = 0.50*self.Vrms_ref  #From IEEE 1557-2018 Category III (Table 16)   
        self.V_LV1 = pvderConfig['V_LV1']*self.Vrms_ref  #From IEEE 1557-2018 Category III (Table 16)   
        self.V_LV2 = pvderConfig['V_LV2']*self.Vrms_ref  #From IEEE 1557-2018  Category III (Table 16)
        self.t_LV0_limit = 1.0/scaling_factor      #Time limit for LV0 zone from IEEE 1557-2018 Category III  (Table 16, page 48)
        self.t_LV1_limit = pvderConfig['t_LV1_limit']/scaling_factor      #Time limit for LV2 zone from IEEE 1557-2018 Category III  (Table 16, page 48)
        self.t_LV2_limit = pvderConfig['t_LV2_limit']/scaling_factor       #Time limit for LV2 zone from IEEE 1557-2018 Category III (Table 16, page 48) 
        
        self.LVRT_INSTANTANEOUS_TRIP =  pvderConfig['LVRT_INSTANTANEOUS_TRIP']
        self.t_reconnect_limit = self.t_LV2_limit  #Time lag before reconnecting
        
        if self.LVRT_INSTANTANEOUS_TRIP:
            self.t_LV0_limit = self.t_LV1_limit = self.t_LV1_limit = 1/60 #Disconnect within one cycle

        assert (self.V_LV1 < self.V_LV2 and self.t_LV1_limit < self.t_LV2_limit) == True, "Voltage level 2 should be greater than Voltage level 1"
        self.check_LVRT_settings()

        #V1 to V2 - zone 2,V1 < - zone 1
    
    def check_LVRT_settings(self):
        """Method to do sanity check of LVRT settings."""
        
        if (self.V_LV0 > self.V_LV1 or self.V_LV1 > self.V_LV2) or \
            (self.V_LV0 < 0.45*self.Vrms_ref or self.V_LV0 > 0.6*self.Vrms_ref) or \
            (self.V_LV1 < 0.65*self.Vrms_ref or self.V_LV1 > 0.75*self.Vrms_ref) or\
            (self.V_LV2 < 0.8*self.Vrms_ref or self.V_LV2 > 0.92*self.Vrms_ref):
            raise ValueError('LVRT voltage limits:{:.2f},{:.2f},{:.2f} are infeasible!'.format(self.V_LV0,self.V_LV1,self.V_LV2))
        
        if (self.t_LV0_limit > self.t_LV1_limit or self.t_LV1_limit > self.t_LV2_limit or \
            self.t_LV0_limit < 1/60 or self.t_LV0_limit > 1 or \
            self.t_LV1_limit < 1 or self.t_LV1_limit > 10 or \
            self.t_LV2_limit < 2 or self.t_LV2_limit > 20) and not self.LVRT_INSTANTANEOUS_TRIP:
            raise ValueError('LVRT ridethrough times {},{},{} are infeasible!'.format(self.t_LV0_limit,self.t_LV1_limit,self.t_LV2_limit))
        
        if self.t_reconnect_limit < 0.4 or self.t_reconnect_limit > self.t_LV2_limit:
            raise ValueError('LVRT reconnect time limit {} is infeasible!'.format(self.t_reconnect_limit))
            
    def show_RT_settings(self,settings_type='LVRT',PER_UNIT=True):
        """Method to show LVRT settings."""
        
        if settings_type not in {'LVRT'}:#,'LFRT'
            raise ValueError('Unknown quantity: ' + str(settings_type))
        
        print('\n______{} - {}_____'.format(self.name,settings_type))
        
        if settings_type ==  'LVRT':
            print('Vrms_ref:{:.2f} V\nV_LV0:{:.2f} V\nV_LV1:{:.2f} V\nV_LV2:{:.2f} V'.format(self.Vrms_ref,self.V_LV0,self.V_LV1,self.V_LV2))
            print('t_LV0:{:.2f} s\nt_LV1:{:.2f} s\nt_LV2:{:.2f} s'.format(self.t_LV0_limit,self.t_LV1_limit,self.t_LV2_limit))
            print('______Flags______')
            print('LVRT_ENABLE:{}\nLVRT_INSTANTANEOUS_TRIP:{}\nLVRT_TRIP:{} '.format(self.LVRT_ENABLE,self.LVRT_INSTANTANEOUS_TRIP,self.LVRT_TRIP))    
    
    def Volt_VAR_logic(self,t):
        """Function for volt_var control."""
        
        assert (self.VOLT_VAR_ENABLE and self.VOLT_WATT_ENABLE) == False, "Volt-VAR and Volt-Watt cannot be active at the same time"
        
        #Volt-VAR logic
        V1_Volt_VAR = self.Vrms_ref - 0.04*self.Vrms_ref #From IEEE 1547-2018 Catergy B (Table 8 - page 39)
        V2_Volt_VAR = self.Vrms_ref - 0.01*self.Vrms_ref #From IEEE 1547-2018 Catergy B (Table 8 - page 39)
        V3_Volt_VAR = self.Vrms_ref + 0.01*self.Vrms_ref #From IEEE 1547-2018 Catergy B (Table 8 - page 39)
        V4_Volt_VAR = self.Vrms_ref + 0.04*self.Vrms_ref #From IEEE 1547-2018 Catergy B (Table 8 - page 39)       
        
        #Select RMS voltage source
        _Vrms_measured = self.Vrms
        
        _del_V = 0.02
        if self.VOLT_VAR_ENABLE == True and (_Vrms_measured < V2_Volt_VAR or _Vrms_measured > V3_Volt_VAR) and t>self.t_stable:
            if self.VOLT_VAR_FLAG == True:
                if (_Vrms_measured > V1_Volt_VAR - _del_V and _Vrms_measured < V4_Volt_VAR + _del_V):
                    
                    Qref = self.Qsetpoint_calc(t)                    
                    
                else:
                    utility_functions.print_to_terminal("Volt-VAR control with reference voltage {:.3f} is deactivated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,_Vrms_measured*self.Vbase))
                    Qref = 0.0
                    self.VOLT_VAR_FLAG = False
            else:
                if ( _Vrms_measured >= V1_Volt_VAR and  _Vrms_measured < V4_Volt_VAR):
                    utility_functions.print_to_terminal("Volt-VAR control with reference voltage {:.3f} is activated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,_Vrms_measured*self.Vbase))
                    self.VOLT_VAR_FLAG = True
                    Qref = self.Qsetpoint_calc(t)
                else:
                    utility_functions.print_to_terminal("Volt-VAR control not activated at {:.3f} s since {:.3f} V is outside range!".format(t,_Vrms_measured*self.Vbase))
                    Qref = 0.0
        else: 
            Qref = 0.0
        
        return Qref
    
    def Qsetpoint_calc(self,t):
        """Function to calculate Qsetpoint."""
        _del_V_to_Q = 200.0 #100.0
        
        _Vrms_measured = self.Vrms
        
        self.Qlimit = math.sqrt(math.pow(self.Sinverter_nominal,2)-math.pow(max(min(self.S.real,self.Sinverter_nominal),-self.Sinverter_nominal),2)) - self.Xf*(abs(self.ia)/math.sqrt(2))**2 #Find Qlimit based on current power output
        
        #utility_functions.print_to_terminal("Q is {:.3f} while Qlimit is {:.3f} at {:.3f} s".format(self.S_PCC.imag*self.Sbase,self.Qlimit*self.Sbase,t))
        
        if _Vrms_measured < self.Vrms_ref:
                Qref = min(self.Qlimit,(self.Vrms_ref -self.Vrms_ref*0.01-_Vrms_measured)*_del_V_to_Q)
                                
        elif _Vrms_measured > self.Vrms_ref:
                Qref = max(-self.Qlimit,(self.Vrms_ref + self.Vrms_ref*0.01-_Vrms_measured)*_del_V_to_Q)
        return Qref
    
    def LVRT(self,t):
        """Function to implement LVRT trip and reconnect logic."""    
        #Select RMS voltage source
        _Vrms_measured = self.Vrms   #Select PCC - LV side voltage     

        if self.LVRT_ENABLE == True and t > self.t_stable and self.LVRT_TRIP == False:

            if self.t_LV1start == 0.0 or self.t_LV2start == 0.0:
                if self.t_LV1start == 0.0 and _Vrms_measured < self.V_LV1:
                    self.t_LV1start = t
                    utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV1start,event_name='LV1_start')

                if self.t_LV2start == 0.0 and self.Vrms < self.V_LV2:
                    self.t_LV2start = t
                    utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV2start,event_name='LV2_start')

            if self.t_LV1start > 0.0 or self.t_LV2start > 0.0:

                if  _Vrms_measured >= self.V_LV1  or  _Vrms_measured >= self.V_LV2:
                    if  _Vrms_measured >= self.V_LV1 and self.t_LV1start > 0.0:
                        utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV1start,event_name='LV1_reset')
                        self.t_LV1start = 0.0
                    if  _Vrms_measured >= self.V_LV2 and self.t_LV2start > 0.0:
                        utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV2start,event_name='LV2_reset')   
                        self.t_LV2start = 0.0

                if t-self.t_LV1start >= self.t_LV1_limit and self.t_LV1start > 0.0:   #Trip inverter if any limit breached
                        utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV1start,event_name='inverter_trip_LV1')
                        #six.print_(t-self.t_LV1start,self.t_LV1_limit)  #debug script
                        self.LVRT_TRIP = True
                        self.t_LV1start = 0.0
                        self.t_LV2start = 0.0
                        self.t_reconnect = 0.0
                elif t-self.t_LV2start >= self.t_LV2_limit and self.t_LV2start > 0.0:
                        utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV2start,event_name='inverter_trip_LV2')
                        self.LVRT_TRIP = True
                        self.t_LV1start = 0.0
                        self.t_LV2start = 0.0
                        self.t_reconnect = 0.0

                if _Vrms_measured < self.V_LV1 and self.t_LV1start > 0.0:
                        utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV1start,event_name='LV1_zone') #,verbose = False
                elif _Vrms_measured < self.V_LV2 and self.t_LV2start > 0.0:
                        utility_functions.print_LVRT_events(t,_Vrms_measured/self.Vrms_ref,self.t_LV2start,event_name='LV2_zone')
                        #six.print_(self.t_LV2start)

        elif self.LVRT_ENABLE == True and t > self.t_stable and self.LVRT_TRIP == True:
            
            #Select RMS voltage source
            _Vrms_measured = self.Vrms   #Select PCC - LV side voltage     

            if self.t_reconnect > 0.0:
                if  _Vrms_measured < self.V_LV2:
                    utility_functions.print_LVRT_events(t, _Vrms_measured/self.Vrms_ref,self.t_reconnect,event_name='reconnect_reset')
                    self.t_reconnect = 0.0
                elif  _Vrms_measured >= self.V_LV2 and t-self.t_reconnect >= self.t_reconnect_limit:
                    utility_functions.print_LVRT_events(t, _Vrms_measured/self.Vrms_ref,self.t_reconnect,event_name='inverter_reconnection')
                    self.LVRT_TRIP = False             #Reset trip flag
                    self.LVRT_RECONNECT = True        #Set reconnect flag
                    self.t_reconnect = 0.0    #Reset timer
                elif _Vrms_measured >= self.V_LV2 and t-self.t_reconnect < self.t_reconnect_limit:
                    utility_functions.print_LVRT_events(t, _Vrms_measured/self.Vrms_ref,self.t_reconnect,event_name='reconnect_zone')


            elif self.t_reconnect == 0.0:
                if  _Vrms_measured >= self.V_LV2:
                    self.t_reconnect = t
                    utility_functions.print_LVRT_events(t, _Vrms_measured/self.Vrms_ref,self.t_reconnect,event_name='reconnect_start')
                else:
                    utility_functions.print_LVRT_events(t, _Vrms_measured/self.Vrms_ref,event_name='inverter_tripped')

        else:
            self.t_LV1start = 0.0
            self.t_LV2start = 0.0
            self.LVRT_TRIP = False
    
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

    
    def FRT(self,t):
        """Frequency ride through and trip logic. """
        t_reconnect_limit = 3.0   #Time lag before reconnecting
        """
        #ERCOT standards from NERC PRC - 024
        F_LF0 = 57.5  #From NERC PRC - 024
        F_LF1 = 58.0  #From NERC PRC - 024
        F_LF2 = 58.4  #From NERC PRC - 024
        F_LF3 = 59.4  #From NERC PRC - 024
        t_LF1_limit = 2.0      #Time limit for LF1 zone
        t_LF2_limit = 30.0      #Time limit for LF2 zone
        t_LF3_limit = 540.0      #Time limit for LF2 zone
        """
        #IEEE 1547-2018 standards
        F_LF1 = 57.0  #From IEEE 1557-2018 Category III    
        F_LF2 = 58.8  #From IEEE 1557-2018 Category III    
        t_LF1_limit = 0.16      #Time limit for LF1 zone From IEEE 1557-2018 Category III (Table 19 and figure H10)   
        t_LF2_limit = 299.0      #Time limit for LF2 zone From IEEE 1557-2018 Category III (Table 19 and figure H10)
        
        #Dummy variables to preserve logic
        F_LF3 = 0.0
        t_LF3_limit = 0.0
               
        del_f =0.02
        assert (F_LF1 < F_LF2) and (t_LF1_limit < t_LF2_limit) == True, "Frequency level 2 should be greater than frequency level 1"

        #Use grid frequency as estimated by PLL
        fgrid = self.we/(2.0*math.pi)
        
        if self.LFRT_ENABLE == True and t >  self.t_stable and self.LFRT_TRIP == False:
            
            if self.LFRT_DEBUG == True:
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

        elif self.LFRT_ENABLE == True and t > self.t_stable and self.LFRT_TRIP == True:

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