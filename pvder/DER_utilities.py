from __future__ import division
import operator
import math
import cmath
import numpy as np
    
from pvder.grid_components import BaseValues
from pvder import utility_functions

class PVDER_ModelUtilities(BaseValues):
    """
       Utility class for single phase and three phase PV-DER model.
    """
   
    Vdcbase = BaseValues.Vbase #DC side base value is same as AC side base value
    
    #Ramp control
    RAMP_ENABLE = False
    RAMP_FLAG = False
    ramp_list = []
    n_steps = 0
    ramp_del_t = 0.5 #0.025 
            
    #Limits
    m_limit = 1.0 #Maximum duty cycle
    
     #Flags
    PRINT_INLINE = False
    VERBOSE = False
    MPPT_ENABLE = False
    
    DO_EXTRA_CALCULATIONS = False #Do calculations not essential to ODE model (useful for debugging)
    Qref_EXTERNAL = False #Allow direct manipulation of Qref from outside ODE model.
    
    Vdc_ref_list = []
    Vdc_ref_total = len(Vdc_ref_list) #Get total events
    Vdc_ref_counter = 0
    del_Vdc_ref = 1.0
        
    @property                         #Decorator used for auto updating
    def Vdc_actual(self):
        """Actual DC link voltage"""
        return min(self.Vdcmpp_max,self.Vdc*self.Vdcbase)  #Calculate actual voltage
    
    #Average duty cycle - Phase A
    @property                         #Decorator used for auto updating
    def ma(self):
        """Phase A duty cycle"""
        
        return self.Kp_GCC*self.ua + self.xa #PI controller equation
        #return utility_functions.m_calc(self.Kp_GCC,self.ua,self.xa)
    
    #Average duty cycle - Phase B
    @property                         #Decorator used for auto updating
    def mb(self):
        """Phase B duty cycle"""
        
        return self.Kp_GCC*self.ub + self.xb #PI controller equation
    
    #Average duty cycle - Phase C
    @property                         #Decorator used for auto updating
    def mc(self):
        """Phase C duty cycle"""
        
        return self.Kp_GCC*self.uc + self.xc #PI controller equation
    
    #Controller outer loop equations (Current set-point)    
    def ia_ref_calc(self):
        """Phase A current reference"""
        return self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ  - self.Kp_Q*(self.Q_ref - self.S_PCC.imag)) #PI controller equation
           
    def ib_ref_calc(self):
        """Phase B current reference"""
        
        return utility_functions.Ub_calc(self.ia_ref)
    
    def ic_ref_calc(self):
        """Phase C current reference"""
        
        return utility_functions.Uc_calc(self.ia_ref)    
    
    def iphload1_calc(self,vph):
        """Current counsumed by load connected at PCC LV side -  Phase A/B/C."""
        
        return vph/self.Zload1
        #return self.va/self.Zload1
    
    def vta_calc(self):
        """Inverter terminal voltage -  Phase A"""
        
        return self.ma*(self.Vdc/2)
        
    def vtb_calc(self):
        """Inverter terminal voltage -  Phase B"""
        
        return self.mb*(self.Vdc/2)
    
    def vtc_calc(self):
        """Inverter terminal voltage -  Phase C"""
        
        return self.mc*(self.Vdc/2)
    
    
    def va_calc(self):
        """PCC - LV side - Phase A"""
        
        if self.standAlone:
            val=((self.grid_model.vag+(self.ia/self.a)*self.grid_model.Z2)/(self.a) +self.ia*self.Z1)*((self.Zload1*self.a*self.a)/((self.a*self.a*(self.Z1+self.Zload1))+self.grid_model.Z2))
        
        else:
            val=self.gridVoltagePhaseA
        
        return val    
    
    #@property
    def vb_calc(self):
        """PCC - LV side - Phase B"""
        
        if self.standAlone:
            if type(self).__name__ == 'SolarPV_DER_SinglePhase':
                val=((self.grid_model.vbg)/(self.a))*((self.Zload1*self.a*self.a)/((self.a*self.a*(self.Z1+self.Zload1))+self.grid_model.Z2))
            
            elif type(self).__name__ == 'SolarPV_DER_ThreePhase':
                val = ((self.grid_model.vbg+(self.ib/self.a)*self.grid_model.Z2)/(self.a) +self.ib*self.Z1)*((self.Zload1*self.a*self.a)/((self.a*self.a*(self.Z1+self.Zload1))+self.grid_model.Z2))
        
        else:
            val=self.gridVoltagePhaseB
        
        return val
    
    #@property
    def vc_calc(self):
        """PCC - LV side - Phase C"""
        
        if self.standAlone:
            if type(self).__name__ == 'SolarPV_DER_SinglePhase':
                val=((self.grid_model.vcg)/(self.a))*((self.Zload1*self.a*self.a)/((self.a*self.a*(self.Z1+self.Zload1))+self.grid_model.Z2))
            
            elif type(self).__name__ == 'SolarPV_DER_ThreePhase':
                val = ((self.grid_model.vcg+(self.ic/self.a)*self.grid_model.Z2)/(self.a) +self.ic*self.Z1)*((self.Zload1*self.a*self.a)/((self.a*self.a*(self.Z1+self.Zload1))+self.grid_model.Z2))
        
        else:
            val=self.gridVoltagePhaseC
        
        return val
    
    def wgrid_calc(self):
        """Frequency of grid voltage source."""
        if self.standAlone:
            val = self.grid_model.wgrid
        else:
            val = self.gridFrequency
        return val
    
    #PLL equation (inverter frequency)
    def we_calc(self):
        """Calculate inverter frequency from PLL."""
        
        return  self.Kp_PLL*(self.vd) + self.xPLL + 2*math.pi*60.0
        
    def update_Ppv(self,t):
        """Update load impedance at PCC-LV side."""
       
        Sinsol_new,Tactual_new = self.events.solar_events(t)
        if abs(self.Sinsol- Sinsol_new) or abs(self.Tactual- Tactual_new) > 0.0:    #Update Iph only if solar insolation changes
             
            self.Sinsol = Sinsol_new
            self.Tactual = Tactual_new
            utility_functions.print_to_terminal("{}:PV module current output changed from {:.3f} A to {:.3f} A at {:.3f}".format(self.name,self.Iph,self.Iph_calc(),t))
            self.Iph = self.Iph_calc()
        
        self.Ppv = self.Ppv_calc(self.Vdc_actual)
    
    def update_Q_Vdc_ref(self,t):
        """Update DC link voltage reference, and reactive power reference."""
        
        self.Q_ref = self.get_Qref(t)
        self.Vdc_ref = self.get_Vdc_ref(t)
    
    def update_Zload1(self,t):
        """Update load impedance at PCC-LV side."""
        
        if self.standAlone:   #Update load at PCC LV side only in stand alone mode
            Zload1_actual_new =  self.events.load_events(t)
            Zload1_new = Zload1_actual_new/BaseValues.Zbase
        
            if abs(self.Zload1- Zload1_new)> 0.0:
                self.Zload1 = Zload1_new
                utility_functions.print_to_terminal("Load at PCC LV side changed from {:.3f} VA to {:.3f} VA at {:.3f}".format(self.S_load1,self.S_load1_calc(),t))
    
    def S_PCCph_calc(self,vph,iph):
        """Inverter apparent power output - phase a/b/c"""
        
        return (1/2)*(vph*iph.conjugate())
        #return (1/2)*(self.va*self.ia.conjugate())*1.0
    
    def show_PV_DER_states(self,quantity='voltage'):
        """Display values of states in the DER model quantities.
        Args:
          quantity: A string ('voltage','current','power','duty cycle') specifying the electrical quantity to be displayed.
        """
        
        if quantity not in {'voltage','current','power','duty cycle'}:
            raise ValueError('Unknown quantity: ' + str(quantity))
        print('\n______{} - {}_____'.format(self.name,quantity.capitalize()))
        
        if quantity ==  'voltage':
            print('Vdc:{:.2f}\nVta:{:.2f} V'.format(self.Vdc*self.Vbase,self.vta*self.Vbase))
            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                print('Vtb:{:.2f} V,Vtb:{:.2f} V\nVtn:{:.2f} V'.format(self.vtb*self.Vbase,self.vtc*self.Vbase,(self.vta+self.vtb+self.vtc)*self.Vbase))
            
            print('Va:{:.2f} V'.format(self.va*self.Vbase))
            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                print('Vb:{:.2f} V,Vc:{:.2f} V\nVn:{:.2f} V'.format(self.vb*self.Vbase,self.vc*self.Vbase,(self.vta+self.vtb+self.vtc)*self.Vbase))
            
            print('Vtrms:{:.2f} V\nVpccrms:{:.2f} V'.format(self.Vtrms*self.Vbase,self.Vrms*self.Vbase))
        
        elif quantity ==  'current':
            print('ia:{:.2f} A'.format(self.ia*self.Ibase))
            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                print('ib:{:.2f} A,ic:{:.2f} A\nIn:{:.2f} A'.format(self.ib*self.Ibase,self.ic*self.Ibase,(self.ia+self.ib+self.ic)*self.Ibase))
            print('Irms:{:.2f} V'.format(self.Irms*self.Ibase))
             
        elif quantity ==  'power':
            print('Ppv:{:.1f} W\nS:{:.1f} VA\nS_PCC:{:.1f} VA'.format(self.Ppv*self.Sbase,self.S*self.Sbase,self.S_PCC*self.Sbase)) 
        
        elif quantity ==  'duty cycle':
            print('ma:{:.2f}'.format(self.ma))
            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                print('mb:{:.2f},mc:{:.2f}\nm0:{:.2f}'.format(self.mb,self.mc,(self.ma+self.mb+self.mc)))
    
    def show_PV_DER_parameters(self,quantity='inverter_ratings'):
        """Display rated values.
        Args:
          quantity: A string ('inverter_ratings','controller_gains','circuit_parameters') specifying the parameter to be displayed.
        """
        
        if quantity not in {'inverter_ratings','controller_gains','circuit_parameters'}:
            raise ValueError('Unknown quantity: ' + str(quantity))
        
        if quantity ==  'inverter_ratings':
            print('Vdcrated:{:.3f} V'.format(self.Vdcrated))
            print('Vtrated (L-G peak):{:.3f} V\nVrated (L-G peak):{:.3f} V'.format((self.Vdcrated/2)*self.m_steady_state,self.Varated))
        
        elif quantity == 'circuit_parameters':
            print('Cdc:{:.9f} F\nLf:{:.6f} H\nRf:{:.3f} Ohm'.format(self.C*self.Cbase,self.Lf*self.Lbase,self.Rf*self.Zbase))
        
        elif quantity == 'controller_gains':
            print('Current controller:\nKp_GCC:{:.3f}, Ki_GCC:{:.3f}, wp:{:.3f}'.format(self.Kp_GCC,self.Ki_GCC,self.wp))
            print('DC link voltage controller:\nKp_DC:{:.3f}, Ki_DC:{:.3f}'.format(self.Kp_DC,self.Ki_DC))
            print('Reactive power controller:\nKp_Q:{:.3f}, Ki_Q:{:.3f}'.format(self.Kp_Q,self.Ki_Q))
            print('PLL controller:\nKp_PLL:{:.3f}, Ki_PLL:{:.3f}'.format(self.Kp_PLL,self.Ki_PLL))     

    def validate_model(self,PRINT_ERROR = True):
        """Compare error between RMS quantities and Phasor quantities."""
        
        #Calculation with phasor quantities
        self.Pf_phasor = self.S_calc().real-self.S_PCC_calc().real  #Active power consumed by filter resistor
        self.Qf_phasor = self.S_calc().imag-self.S_PCC_calc().imag  #Reactive power consumed by filter inductor 
        
        #Caculation with RMS quantities        
        if type(self).__name__ == 'SolarPV_DER_SinglePhase':
            _phases = 1
        elif type(self).__name__ == 'SolarPV_DER_ThreePhase':
            _phases = 3
        
        self.Pf_RMS = _phases*((self.Irms)**2)*self.Rf   #Active power consumed by filter resistor
        self.Qf_RMS = _phases*((self.Irms)**2)*self.Xf   #Reactive power consumed by filter inductor 
        
        #Calculation with phasor quantities
        self.Pt_phasor = self.S_calc().real   #Active power output at inverter terminal
        self.Qt_phasor = self.S_calc().imag   #Reactive power output at inverter terminal
                
        #Caculation with RMS quantities 
        ra1,pha1 = cmath.polar(self.vta)
        ra2,pha2 = cmath.polar(self.ia)        
        
        self.Pt_RMS = (abs(self.vta)/math.sqrt(2))*(abs(self.ia)/math.sqrt(2))*math.cos(pha1-pha2) #Active power at inverter terminal
        
        self.Qt_RMS = (abs(self.vta)/math.sqrt(2))*(abs(self.ia)/math.sqrt(2))*math.sin(pha1-pha2) #Reactive power output        
        
        if type(self).__name__ == 'SolarPV_DER_ThreePhase':
            rb1,phb1 = cmath.polar(self.vtb)
            rb2,phb2 = cmath.polar(self.ib) 
            rc1,phc1 = cmath.polar(self.vtc)
            rc2,phc2 = cmath.polar(self.ic)
        
            self.Pt_RMS = self.Pt_RMS +\
                         (abs(self.vtb)/math.sqrt(2))*(abs(self.ib)/math.sqrt(2))*math.cos(phb1-phb2) +\
                         (abs(self.vtc)/math.sqrt(2))*(abs(self.ic)/math.sqrt(2))*math.cos(phc1-phc2) #Active power at inverter terminal
        
            self.Qt_RMS = self.Qt_RMS +\
                          (abs(self.vtb)/math.sqrt(2))*(abs(self.ib)/math.sqrt(2))*math.sin(phb1-phb2) +\
                          (abs(self.vtc)/math.sqrt(2))*(abs(self.ic)/math.sqrt(2))*math.sin(phc1-phc2) #Reactive power output
        
        #self.Pt_RMS = 3*(self.Vtrms)*(self.Irms)*math.cos(ph1-ph2) #Active power output at inverter terminal
        #self.Qt_RMS = 3*(self.Vtrms)*(self.Irms)*math.sin(ph1-ph2) #Reactive power output at inverter terminal
        
        if PRINT_ERROR:
            print('Active power output error:{:.4f}\nReactive power output error:{:.4f}'.format(abs(self.Pt_phasor-self.Pt_RMS),abs(self.Qt_phasor-self.Qt_RMS)))    
            print('Inverter filter active power loss error:{:.4f}\nInverter filter reactive power loss error:{:.4f}'.format(abs(self.Pf_phasor-self.Pf_RMS),abs(self.Qf_phasor-self.Qf_RMS)))
    
    
    def get_Qref(self,t):
        """Output reactive power set-point."""
        if self.VOLT_VAR_ENABLE == True:
            _Qref= self.Volt_VAR_logic(t)
        elif self.Qref_EXTERNAL:
            _Qref = self.Q_ref
        else:
            _Qref = 0.0/self.Sbase
        return _Qref
    
    def MPP_table(self):
        """Method to output Vdc reference corresponding to MPP at different insolation levels values."""
        
        if self.USE_POLYNOMIAL_MPP == True:
             _Vdcmpp = np.polyval(self.z , self.Sinsol)
        else:
            _Vdcmpp = self.Vdcrated

        _Vdcmpp=  max(min(_Vdcmpp,self.Vdcmpp_max),self.Vdcmpp_min)
        
        return _Vdcmpp/self.Vdcbase
  
    def add_Vdc_ref(self,t,Vdc_ref):
        """Add new solar event."""
        
        t = float(t)
        Vdc_ref = float(Vdc_ref)
        
        if Vdc_ref <self.Vdcmpp_min or Vdc_ref > self.Vdcmpp_max:
            raise ValueError('{} V is not a valid value for DC link voltage!'.format(Vdc_ref))       
        
        for ref in self.Vdc_ref_list:
            if t==ref['t']:
                print('Removing existing Vdc_ref at {:.2f}!'.format(event['t']))
                self.Vdc_ref_list.remove(ref)   # in {}Remove exi,self.events_IDsting event at same time stamp
        
        print('Adding new Vdc reference at {:.2f} s'.format(t))
        self.Vdc_ref_list.append({'t':t,'Vdc_ref':Vdc_ref/self.Vdcbase})  #Append new event to existing event list
        self.Vdc_ref_list.sort(key=operator.itemgetter('t'))  #Sort new events list
        self.Vdc_ref_total = len(self.Vdc_ref_list) #Get total events
    
    def get_Vdc_ref(self,t):
        """Output Vdc reference."""
        if self.Vdc_ref_list: #Check whether list is empty
            
            if t<self.Vdc_ref_list[0]['t']: 
                Vdc_ref = self.Vdc_ref
                                
            elif t<self.Vdc_ref_list[self.Vdc_ref_counter]['t']  and self.Vdc_ref_counter >=1:
                Vdc_ref = self.Vdc_ref_list[self.Vdc_ref_counter-1]['Vdc_ref']
               
            elif t>=self.Vdc_ref_list[self.Vdc_ref_counter]['t']:
                Vdc_ref = self.Vdc_ref_list[self.Vdc_ref_counter]['Vdc_ref']
                
                self.Vdc_ref_counter = min(self.Vdc_ref_total-1,self.Vdc_ref_counter+1)
        else:
            Vdc_ref = self.Vdc_ref
            
        return Vdc_ref
    
    def Vdc_ref_ramp(self,tstart,Vdc_ref_target):
        """Create a ramp signal for voltage reference that ramps at 1 V/s.
        Args:
           tstart: A scalr specifying start time of ramp in seconds.
           Vdc_ref_target: A scalar specifying target Vdc reference seconds in volts.
        """
        
        Vdc_ref_start = self.get_Vdc_ref(t=tstart)*self.Vdcbase
        if abs(Vdc_ref_start-Vdc_ref_target) <= self.del_Vdc_ref:
            self.add_Vdc_ref(t=tstart,Vdc_ref=Vdc_ref_target)
            
        else:
            Vdc_ref_range = np.arange(Vdc_ref_start+self.del_Vdc_ref,Vdc_ref_target+self.del_Vdc_ref,self.del_Vdc_ref)
            trange = np.arange(tstart,tstart+len(Vdc_ref_range),1.0)
            for i,Vdc_ref in enumerate(Vdc_ref_range):
                self.add_Vdc_ref(t=trange[i],Vdc_ref=Vdc_ref)
            
    def show_references(self):
        """Print references."""
        
        print('Showing all references in {}!'.format(self.name))
        print('Total references:{}'.format(len(self.Vdc_ref_list)))
        if self.Vdc_ref_list:
            for ref in self.Vdc_ref_list:
                print('t:{:.3f},Vdc_ref:{:.3f} V'.format(ref['t'],ref['Vdc_ref']*self.Vdcbase))
        else:
            print("No Vdc references!!!")
    
    def reset_reference_counters(self):
        self.Vdc_ref_counter = 0
        
        print('Reference event counters reset!')
        
    def initialize_jacobian(self):
        """Create a Jacobian matrix with zero values."""
        
        self.J = np.zeros((self.n_total_ODE,self.n_total_ODE))
        self.varInd={}; n=0
        
        if type(self).__name__ == 'SolarPV_DER_SinglePhase':
            state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
                          'Vdc','xDC','xQ','xPLL','wte']
        
        elif type(self).__name__ == 'SolarPV_DER_ThreePhase':
            state_list = ['iaR','iaI','xaR','xaI','uaR','uaI',
                          'ibR','ibI','xbR','xbI','ubR','ubI',
                          'icR','icI','xcR','xcI','ucR','ucI',
                          'Vdc','xDC','xQ','xPLL','wte']            
            
        for entry in state_list:
            self.varInd[entry]=n
            n+=1