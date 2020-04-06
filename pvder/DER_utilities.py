"""Code containing utilities used by PV-DER model instances."""

from __future__ import division
import operator
import logging
import pprint
import json
import pickle

import math
import cmath
import numpy as np

from pvder.utility_classes import Logging
from pvder.grid_components import BaseValues
from pvder import utility_functions
from pvder import config, templates

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
    del_Vdc_ref = config.DEFAULT_del_Vdc_ref
    del_t_Vdc_ref = config.DEFAULT_del_t_Vdc_ref
    
    #Grid frequency estimate variables
    use_frequency_estimate = config.DEFAULT_USE_FREQUENCY_ESTIMATE
    _del_t_frequency_estimate =  config.DEFAULT_DELTA_T
    _t_estimate_frequency_previous = 0.0  
    _westimate = 2*math.pi*60.0
        
    @property                         #Decorator used for auto updating
    def Vdc_actual(self):
        """Actual DC link voltage.
                
        Returns:
            float: DC link voltage in Volts.
        """
        
        return min(self.Vdcmpp_max,self.Vdc*self.Vdcbase)  #Calculate actual voltage
    
    #Average duty cycle - Phase A
    @property                         #Decorator used for auto updating
    def ma(self):
        """Phase A duty cycle.
        
        Returns:
            complex: Duty cycle.
        
        """
        
        return self.Kp_GCC*self.ua + self.xa #PI controller equation
        #return utility_functions.m_calc(self.Kp_GCC,self.ua,self.xa)
    
    #Average duty cycle - Phase B
    @property                         #Decorator used for auto updating
    def mb(self):
        """Phase B duty cycle.
        
        Returns:
            complex: Duty cycle.
        
        """        
        
        return self.Kp_GCC*self.ub + self.xb #PI controller equation
    
    #Average duty cycle - Phase C
    @property                         #Decorator used for auto updating
    def mc(self):
        """Phase C duty cycle.
        
        Returns:
            complex: Duty cycle.
        
        """
        
        return self.Kp_GCC*self.uc + self.xc #PI controller equation
    
    def get_DER_config(self,config_file,DER_id):
        """Check DER ID in config file."""
        
        config_dict = self.read_config(config_file) #Get bounds dict
        
        available_ids = list(config_dict.keys())
        if DER_id in available_ids:
            print('DER configuration with ID:{} was found in {}'.format(DER_id,config_file))
        else:
            raise KeyError('DER configuration with ID:{} could not be found in {}! - Available IDs are:{}'.format(DER_id,config_file,available_ids))
        
        return config_dict[DER_id]
    
    def get_pvderconfig(self,pvderConfig):
        """Check DER ID in config file."""
        
        if pvderConfig is not None:
            if not isinstance(pvderConfig,dict):
                raise ValueError('Expected pvderConfig to be a dictionary but got {}'.format(type(pvderConfig)))
        else:
            pvderConfig = {}    
        
        return pvderConfig
    
    #Controller outer loop equations (Current set-point)    
    def ia_ref_calc(self):
        """Phase A current reference"""
        
        return self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ  - self.Kp_Q*(self.Q_ref - self.S_PCC.imag)) #PI controller equation
    
    def ia_ref_constantVdc_calc(self):
        """Phase A current reference for constant Vdc"""
        
        return self.xDC + self.Kp_DC*(self.Ppv -  self.S.real) + 1j*(self.xQ  - self.Kp_Q*(self.Q_ref - self.S_PCC.imag)) #PI controller equation
      
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
    
    def wgrid_calc(self,t):
        """Frequency of grid voltage source."""
        
        if self.use_frequency_estimate:
            val = self.wgrid_estimate(t)
        elif self.standAlone:
            val = self.grid_model.wgrid
        else:
            val = self.gridFrequency
        return val
    
    def wgrid_estimate(self,t):
        """Estimate frequency from phasor angle."""
                
        if t > self._t_estimate_frequency_previous: #Prevent going back in time 
            if self.standAlone:
                _,phia_new = cmath.polar(self.grid_model.vag)
                _,phia_previous = cmath.polar(self._vag_previous)
                if abs(phia_previous - phia_new)>0:
                    self._vag_previous = self.grid_model.vag                   
                    self._t_estimate_frequency_previous = t
                    self._westimate = self.wgrid_estimate_calc(t,phia_new,phia_previous)                   
        
            else:
                _,phia_new = cmath.polar(self.va)
                _,phia_previous = cmath.polar(self._va_previous)                    
                if abs(phia_previous - phia_new)>0:
                    self._va_previous = self.va
                    self._t_estimate_frequency_previous = t
                    self._westimate = self.wgrid_estimate_calc(t,phia_new,phia_previous)
        
        return self._westimate   
    
    def wgrid_estimate_calc(self,t,phi_new,phi_previous):
        """Estimate frequency."""
        
        del_f = (phi_new-phi_previous)/((2*math.pi)*self._del_t_frequency_estimate) #(theta_N - theta_N-1)/(2*pi*dt)\
        festimate = self.wgrid_measured/(2*math.pi) + del_f
        self.logger.debug('t:{}:{}:Phase angle changed from {:.4f} rad to {:.4f} rad -> Estimated frequency from phase angle change:{:.3f} Hz'.format(t,self.name,phi_previous,phi_new,festimate))
        
        return 2*math.pi*(festimate)
    
    #PLL equation (inverter frequency)
    def we_calc(self):
        """Calculate inverter frequency from PLL."""
        
        return  self.Kp_PLL*(self.vd) + self.xPLL + 2*math.pi*60.0
        
    def update_Ppv(self,t):
        """Update PV module power output based on solar events and DC link voltage."""
       
        Sinsol_new,Tactual_new = self.events.solar_events(t)
        if abs(self.Sinsol- Sinsol_new) or abs(self.Tactual- Tactual_new) > 0.0:    #Update Iph only if solar insolation changes
             
            self.Sinsol = Sinsol_new
            self.Tactual = Tactual_new
            utility_functions.print_to_terminal("{}:PV module current output changed from {:.3f} A to {:.3f} A at {:.3f} s".format(self.name,self.Iph,self.Iph_calc(),t))
            self.Iph = self.Iph_calc()
        
        self.Ppv = self.Ppv_calc(self.Vdc_actual)
    
    def update_Qref(self,t):
        """Update reactive power set-point."""
        
        if self.VOLT_VAR_ENABLE:
            Qref= self.Volt_VAR_logic(t)
        elif self.Qref_EXTERNAL:
            Qref = self.Q_ref
        else:
            Qref = 0.0/self.Sbase
        
        self.Q_ref = Qref
        
        return Qref
    
    def update_Vdc_ref(self,t):
        """Update DC link voltage reference."""
        
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
    
    def show_PV_DER_states(self,quantity='voltage'):
        """Display values of states in the DER model quantities.
        
        Arguments
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
    
    def show_PV_DER_parameters(self,parameter_type='inverter_ratings'):
        """Display rated values.
        
        Args:
          parameter_type (str): A string ('inverter_ratings','controller_gains','circuit_parameters') specifying the parameter to be displayed.
        
        """
        
        if parameter_type not in {'module_parameters','inverter_ratings','controller_gains','circuit_parameters','all'}:
            raise ValueError('Unknown quantity: ' + str(parameter_type))
        
        print('Parameters for DER with ID:{}'.format(self.parameter_ID))
        if parameter_type ==  'module_parameters' or parameter_type ==  'all':
            print('Np:{},Ns:{}'.format(self.Np,self.Ns))
            print('Vdcmpp0:{:.3f} V\nVdcmpp_min:{:.3f} V\nVdcmpp_max:{:.3f} V'.format(self.Vdcmpp0,self.Vdcmpp_min,self.Vdcmpp_max))
        
        if parameter_type ==  'inverter_ratings' or parameter_type ==  'all':
            print('Srated:{:.3f} VA'.format(self.Sinverter_rated))
            print('Vdcrated:{:.3f} V'.format(self.Vdcrated))
            print('Vtrated (L-G peak):{:.3f} V\nVrated (L-G peak):{:.3f} V'.format((self.Vdcrated/2)*self.m_steady_state,self.Varated))
        
        if parameter_type == 'circuit_parameters' or parameter_type ==  'all':
            print('Cdc:{:.9f} F\nLf:{:.6f} H\nRf:{:.3f} Ohm'.format(self.C*self.Cbase,self.Lf*self.Lbase,self.Rf*self.Zbase))
        
        if parameter_type == 'controller_gains' or parameter_type ==  'all':
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
        
    def set_Vdc_ref(self):
        """Return the correct Vdc reference voltage."""
        
        if self.MPPT_ENABLE:
            Vdc_ref = self.Vdcmpp/self.Vdcbase
        else:
            Vdc_ref = self.Vdcnominal
        
        return Vdc_ref        
    
    def MPP_table(self):
        """Method to output Vdc reference corresponding to MPP at different insolation levels values."""
        
        if self.USE_POLYNOMIAL_MPP:
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
                print('Removing existing Vdc_ref at {:.2f}!'.format(t))
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
        
        Arguments:
          tstart (float): A scalar specifying start time of ramp in seconds.
          Vdc_ref_target (float): A scalar specifying target Vdc reference seconds in volts.
        """
        
        Vdc_ref_start = self.get_Vdc_ref(t=tstart)*self.Vdcbase
        if abs(Vdc_ref_start-Vdc_ref_target) <= self.del_Vdc_ref:
            self.add_Vdc_ref(t=tstart,Vdc_ref=Vdc_ref_target)
            
        else:
            if Vdc_ref_start>Vdc_ref_target:
                self.del_Vdc_ref = -self.del_Vdc_ref
            Vdc_ref_range = np.arange(Vdc_ref_start+self.del_Vdc_ref,Vdc_ref_target+self.del_Vdc_ref,self.del_Vdc_ref)
            Vdc_ref_range[-1]=Vdc_ref_target
            trange = np.arange(tstart,tstart+len(Vdc_ref_range),self.del_t_Vdc_ref)
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
            print("{}:No Vdc references!!!".format(self.name))
    
    def reset_reference_counters(self):
        """Reset counter for reference change events."""
        
        self.Vdc_ref_counter = 0
        
        self.logger.debug('{}:Reference event counters reset!'.format(self.name))        
    
    def create_parameter_dict(self,parameter_ID):
        """Create a DER mode."""
        
        assert isinstance(parameter_ID, str), 'Expected parameter_ID to be a string, but got {}!'.format(type(parameter_ID))
                
        default_ID = self.get_default_parameter_ID()
        
        self.module_parameters[parameter_ID] = dict.fromkeys(list(self.module_parameters[default_ID].keys()), None)
        self.inverter_ratings[parameter_ID] = dict.fromkeys(list(self.inverter_ratings[default_ID].keys()), None)
        self.circuit_parameters[parameter_ID] = dict.fromkeys(list(self.circuit_parameters[default_ID].keys()), None)
        self.controller_gains[parameter_ID] = dict.fromkeys(list(self.controller_gains[default_ID].keys()), None)
        self.steadystate_values[parameter_ID] = dict.fromkeys(list(self.steadystate_values[default_ID].keys()), None)
        
        self.logger.debug('{}:Creating parameter dicitonary with ID {}!'.format(self.name,parameter_ID))            
    
    def get_default_parameter_ID(self):
        """Return default parameter ID."""
        
        if type(self).__name__ == 'SolarPV_DER_SinglePhase':
            default_ID = config.DEFAULT_PARAMETER_ID_1PH  
        elif type(self).__name__ == 'SolarPV_DER_ThreePhase':
            default_ID = config.DEFAULT_PARAMETER_ID_3PH 
        
        return default_ID
    
    def initialize_parameter_dict(self,parameter_ID,source_parameter_ID):
        """Initialize a new parameter dictinary with values from an existing parameter dictionary."""
        
        self.create_parameter_dict(parameter_ID)
        
        self.update_parameter_dict(parameter_ID,'module_parameters',self.module_parameters[source_parameter_ID])
        self.update_parameter_dict(parameter_ID,'inverter_ratings',self.inverter_ratings[source_parameter_ID])
        self.update_parameter_dict(parameter_ID,'circuit_parameters',self.circuit_parameters[source_parameter_ID])  
        self.update_parameter_dict(parameter_ID,'controller_gains',self.controller_gains[source_parameter_ID])
        self.update_parameter_dict(parameter_ID,'steadystate_values',self.steadystate_values[source_parameter_ID])
        
        self.logger.info('{}:Created and initialized new parameter dicitonary {} with source dictionary {}.'.format(self.name,parameter_ID,source_parameter_ID))
    
    def update_parameter_dict(self,parameter_ID,parameter_type,parameter_dict):
        """Update parameters."""
                
        if parameter_type not in templates.DER_design_template:
            raise ValueError('Unknown parameter type: ' + str(parameter_type))
        
        for parameter in parameter_dict.keys():
            
            if parameter_type == 'module_parameters':            
                self.module_parameters[parameter_ID][parameter] = parameter_dict[parameter]        
                      
            elif parameter_type == 'inverter_ratings':        
                self.inverter_ratings[parameter_ID][parameter] = parameter_dict[parameter]
            
            elif parameter_type == 'circuit_parameters':        
                self.circuit_parameters[parameter_ID][parameter] = parameter_dict[parameter]
                
            elif parameter_type == 'controller_parameters':        
                self.controller_parameters[parameter_ID][parameter] = parameter_dict[parameter]
            
            elif parameter_type == 'steadystate_values':        
                self.steadystate_values[parameter_ID][parameter] = parameter_dict[parameter]
            
            else:
                print('{} is invalid parameter!'.format(parameter_type))
                
            self.logger.debug('{}:Updating {} in parameter dicitonary {} with {}!'.format(self.name,parameter,parameter_ID,parameter_dict[parameter]))
        
        if self.parameter_ID == parameter_ID:
            self.initialize_DER()
    
    def show_parameter_dictionaries(self):
        """Show all parameter dictionary types and their ID's."""
        
        print('-----Parameter dictionary: Parameter IDs-----')
        
        utility_functions.print_dictionary_keys(self.module_parameters,'module_parameters')
        utility_functions.print_dictionary_keys(self.inverter_ratings,'inverter_ratings')
        utility_functions.print_dictionary_keys(self.circuit_parameters,'circuit_parameters')
        utility_functions.print_dictionary_keys(self.controller_parameters,'controller_parameters')
        utility_functions.print_dictionary_keys(self.steadystate_values,'steadystate_values')
    
    def show_parameter_types(self):
        """Show all parameters within all parameter dictionary types."""
        
        print('-----Parameter dictionary: Parameter types-----')        
        
        key1 = list(self.module_parameters.keys())[0]
        key2 = list(self.inverter_ratings.keys())[0]
        
        utility_functions.print_dictionary_keys(self.module_parameters[key1],'module_parameters')
        utility_functions.print_dictionary_keys(self.inverter_ratings[key2],'inverter_ratings')
        utility_functions.print_dictionary_keys(self.circuit_parameters[key2],'circuit_parameters')
        utility_functions.print_dictionary_keys(self.controller_parameters[key2],'controller_parameters')
        utility_functions.print_dictionary_keys(self.steadystate_values[key2],'steadystate_values')
    
    def get_parameter_dictionary(self,parameter_type,parameter_ID,SHOW_DICTIONARY=True):
        """Return parameter dictionary for specified parameter type and parameter ID.
        
        Args:
          parameter_type (str): Specify type of parameter.
          parameter_ID (str): Specify parameter ID or 'all'.
          SHOW_DICTIONARY (bool): Print the dictionary.
          
        Returns:
             dict: Parameters and their values
        """
        
        if parameter_type not in {'module_parameters','inverter_ratings','controller_gains','circuit_parameters','all'}:
            raise ValueError('Unknown parameter type: ' + str(parameter_type))
                
        parameter_dict = {}
        
        if parameter_type == 'module_parameters' or parameter_type == 'all':
            parameter_dict.update(self.module_parameters[parameter_ID]) 
                      
        if parameter_type == 'inverter_ratings' or parameter_type == 'all':
            parameter_dict.update(self.inverter_ratings[parameter_ID])            
            
        if parameter_type == 'circuit_parameters' or parameter_type == 'all':
            parameter_dict.update(self.circuit_parameters[parameter_ID])                
                
        if parameter_type == 'controller_parameters' or parameter_type == 'all':
            parameter_dict.update(self.controller_parameters[parameter_ID])
            
        if parameter_type == 'steadystate_values' or parameter_type == 'all':
            parameter_dict.update(self.steadystate_values[parameter_ID])
        
        if SHOW_DICTIONARY:
            self.pp.pprint(parameter_dict)
        
        return parameter_dict     
    
    def save_parameter_dictionary(self,parameter_ID,save_format='pickle',SHOW_DICTIONARY=False):
        """Save parameter dictionary."""
    
        parameter_dict = self.get_parameter_dictionary(parameter_type='all',parameter_ID=parameter_ID,SHOW_DICTIONARY=SHOW_DICTIONARY)
        
        if save_format == 'pickle':
            file_name = parameter_ID + '.pkl'
            
            pickle_out = open(file_name,"wb")
            pickle.dump(parameter_dict, pickle_out)
            pickle_out.close()
            
            self.logger.info('{}:Saved all the parameter dicitonaries as a {} file in {}.'.format(self.name,save_format,file_name))
        
        #elif save_format == 'json':
        #    file_name = parameter_ID + '.json'            
        #    json = json.dumps(parameter_dict)
        #    f = open(file_name,"w")
        #    f.write(json)
        #    f.close()    
        
        else:
            print('Unknown file format!')
            
        return file_name
            
    def load_parameter_dictionary(self,file_name): 
        """Load parameter dictionary from saved file."""
        
        pickle_in = open(file_name,"rb")
        
        parameter_dict = pickle.load(pickle_in)
        
        if isinstance(parameter_dict,dict):
            print('Read following dictionary from {}:'.format(file_name))
            self.pp.pprint(parameter_dict)
            dict_name = file_name.split('.')[0] 
            
            self.logger.debug('{}:Loading parameters into DER parameter dictionary...'.format(self.name))
            
            self.initialize_parameter_dict(parameter_ID = dict_name,source_parameter_ID=self.get_default_parameter_ID())
            self.update_parameter_dict(parameter_ID = dict_name,parameter_type='module_parameters',parameter_dict = parameter_dict)
            self.update_parameter_dict(parameter_ID = dict_name,parameter_type='inverter_ratings',parameter_dict = parameter_dict)
            self.update_parameter_dict(parameter_ID = dict_name,parameter_type='circuit_parameters',parameter_dict = parameter_dict)
            self.update_parameter_dict(parameter_ID = dict_name,parameter_type='controller_parameters',parameter_dict = parameter_dict)
            self.update_parameter_dict(parameter_ID = dict_name,parameter_type='steadystate_values',parameter_dict = parameter_dict)
            
            self.logger.info('{}:Succesfully loaded parameters from {} into DER parameter dictionary with parameter ID {}.'.format(self.name,file_name,dict_name))
            
        else:
            raise ValueError('Expected to read dictionary but found {}!'.format(type(parameter_dict)))
           
        return parameter_dict            
            