from __future__ import division
import math
import cmath

from pvder.grid_components import Grid
from pvder import utility_functions

class PVDER_ModelUtilities(Grid):
    """
       Utility class for single phase and three phase PV-DER model.
    """
   
    Vdcbase = Grid.Vbase #DC side base value is same as AC side base value
    
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
        
    @property                         #Decorator used for auto updating
    def Vdc_actual(self):
        """Actual DC link voltage"""
        return min(self.Vdcrated,self.Vdc*self.Vdcbase)  #Calculate actual voltage
    
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
    
    def S_PCCph_calc(self,vph,iph):
        """Inverter apparent power output - phase a/b/c"""
        
        return (1/2)*(vph*iph.conjugate())
        #return (1/2)*(self.va*self.ia.conjugate())*1.0
    
    def show_PV_DER_states(self,quantity='voltage'):
        """Display values of voltage and current quantities."""
        
        if quantity not in {'voltage','current','power','duty_cycle'}:
            raise ValueError('Unknown quantity: ' + str(quantity))
        print('\n______{} - {}_____'.format(self.name,quantity.capitalize()))
        
        if quantity ==  'voltage':
            print('Vdc:{:.2f}\nVta:{:.2f} V'.format(self.Vdc*self.Vbase,self.vta*self.Vbase))
            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                print('Vtb:{:.2f} V,Vtb:{:.2f} V\nVtn:{:.2f} V'.format(self.vtb*self.Vbase,self.vtc*self.Vbase,(self.vta+self.vtb+self.vtc)*self.Vbase))
            print('Vtrms:{:.2f} V\nVpccrms:{:.2f} V'.format(self.Vtrms*self.Vbase,self.Vrms*self.Vbase))
        
        elif quantity ==  'current':
            print('ia:{:.2f} A'.format(self.ia*self.Ibase))
            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                print('ib:{:.2f} A,ic:{:.2f} A\nIn:{:.2f} A'.format(self.ib*self.Ibase,self.ic*self.Ibase,(self.ia+self.ib+self.ic)*self.Ibase))
            print('Irms:{:.2f} V'.format(self.Irms*self.Ibase))
             
        elif quantity ==  'power':
            print('Ppv:{:.1f} W\nS:{:.1f} VA\nS_PCC:{:.1f} VA'.format(self.Ppv*self.Sbase,self.S*self.Sbase,self.S_PCC*self.Sbase)) 
        
        elif quantity ==  'duty_cycle':
            print('ma:{:.2f}'.format(self.ma))
            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                print('mb:{:.2f},mc:{:.2f}\nm0:{:.2f}'.format(self.ib,self.ic,(self.ma+self.mb+self.mc)))
    
    def show_PV_DER_parameters(self,quantity='inverter_ratings'):
        """Display rated values."""
        
        if quantity not in {'inverter_ratings','controller_gains','circuit_parameters'}:
            raise ValueError('Unknown quantity: ' + str(quantity))
        
        if quantity ==  'inverter_ratings':
            print('Vdcrated:{:.3f} V'.format(self.Vdcrated))
            print('Vtrated (L-G peak):{:.3f} V\nVrated (L-G peak):{:.3f} V'.format((self.Vdcrated/2)*self.m_steady_state,self.Varated))
        
        elif quantity == 'circuit_parameters':
            print('Cdc:{:.9f} F\nLf:{:.6f} H\nRf:{:.3f} Ohm'.format(self.C*self.Cbase,self.Lf*self.Lbase,self.Rf*self.Zbase))
        
        elif quantity == 'controller_parameters':
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
    
    def MPP_table(self):
        """Method to output Vdc reference corresponding to MPP at different insolation levels values."""
        
        if self.USE_POLYNOMIAL_MPP == True:
             _Vdcmpp = np.polyval(self.z , self.Sinsol)
        else:
            _Vdcmpp = self.Vdcrated

        _Vdcmpp=  max(min(_Vdcmpp,self.Vdcmpp_max),self.Vdcmpp_min)
        
        return _Vdcmpp/self.Vdcbase
    
    def Vdc_ramp(self,t):
        """Method to ramp up Vdc setpoint changes."""
        _del_Vdc = 0.5
        _min_Vdc_ramp = 5.0/self.Vdcbase        
        
        if self.RAMP_ENABLE == True and t > self.t_stable:
             
            if len(self.ramp_list) == 0 and self.Vdc_ref != self.Vdc_ref_new and abs(self.Vdc_ref - self.Vdc_ref_new) > _min_Vdc_ramp:
                self.n_steps = 0
                self.t_ramp_start = t
                print('Vdc_ref will ramp from {:.2f} to {:.2f} starting at {:.4f} s'.format(self.Vdc_ref*self.Vdcbase,self.Vdc_ref_new*self.Vdcbase,self.t_ramp_start))

                if self.Vdc_ref_new  < self.Vdc_ref:
                    self.ramp_list = np.append(np.arange(self.Vdc_ref*self.Vdcbase,self.Vdc_ref_new*self.Vdcbase,-_del_Vdc),self.Vdc_ref_new*self.Vdcbase)
                    #utility_functions.print_to_terminal('Creating smaller to larger ramp_list:{}'.format(self.ramp_list))
                else:
                    self.ramp_list = np.append(np.arange(self.Vdc_ref*self.Vdcbase,self.Vdc_ref_new*self.Vdcbase,_del_Vdc),self.Vdc_ref_new*self.Vdcbase)
                    #utility_functions.print_to_terminal('Creating larger to smaller ramp_list:{}'.format(self.ramp_list))
                
                self.ramp_list = np.delete(self.ramp_list, 0)
                self.ramp_list = self.ramp_list/self.Vdcbase
                utility_functions.print_to_terminal('Creating new ramp_list:{} of length {}'.format(self.ramp_list*self.Vdcbase,len(self.ramp_list)))                
               
                self.Vdc_ref =  self.ramp_list[0]
                self.n_steps = self.n_steps+1
                self.RAMP_FLAG = True
                utility_functions.print_to_terminal('Vdc ramp started at {:.3f} s and is expected to end at {:.3f} s'.format(t,t+len(self.ramp_list)*self.ramp_del_t))
                #utility_functions.print_to_terminal('Ramp delta t is {} s'.format(self.ramp_del_t))
                utility_functions.print_to_terminal('Ramp step:{},Vdc_ref is updated to {:.4f} at {:.4f}!'.format(self.n_steps,self.Vdc_ref*self.Vdcbase,t))
                
            elif len(self.ramp_list) >= 1 and self.Vdc_ref != self.Vdc_ref_new:
                  
                if t> self.t_ramp_start + self.ramp_del_t*self.n_steps:
                    self.ramp_list = np.delete(self.ramp_list, 0)
                    self.Vdc_ref =  self.ramp_list[0]
                   
                    self.n_steps = self.n_steps+1
                    #utility_functions.print_to_terminal('Ramp step:{},Vdc_ref is updated to {:.4f} at {:.4f}!'.format(self.n_steps,self.Vdc_ref*self.Vdcbase,t))
                   
                    if self.Vdc_ref == self.Vdc_ref_new:
                        self.ramp_list = np.delete(self.ramp_list, 0)
                else:
                    self.Vdc_ref = self.Vdc_ref 
                
            elif len(self.ramp_list) == 0 and self.n_steps !=0 :
                self.Vdc_ref =   self.Vdc_ref_new
                print('Ramp list has been emptied at {:.4f}! Vdc_ref is {:.4f} and Vdc is {:.4f} !'.format(t,self.Vdc_ref*self.Vdcbase,self.Vdc*self.Vdcbase))
                self.n_steps =0 
                self.RAMP_FLAG = False
            else:
                if abs(self.Vdc_ref - self.Vdc_ref_new) < 2.0/self.Vdcbase:
                   self.Vdc_ref =   self.Vdc_ref_new
                else:
                   self.Vdc_ref =   self.Vdc_ref
        else:
            self.Vdc_ref = self.Vdc_ref     
    