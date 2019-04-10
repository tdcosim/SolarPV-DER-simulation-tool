from __future__ import division
import math
import numpy as np
from scipy.optimize import fsolve, minimize

from pvder.grid_components import Grid
from pvder import utility_functions

class PVDER_SetupUtilities():
    """
       Utility class for error checking during model initialization.
    """
   
    def initialize_DER(self,Sinverter_rated):
        
        if str(int(Sinverter_rated/1e3)) in self.Sinverter_list:
            print('Creating PV inverter instance for DER with rating:' + str(int(Sinverter_rated/1e3)) + ' kVA')
            self.Sinverter_rated = Sinverter_rated #Inverter rating in kVA
            self.Sinverter_nominal = (self.Sinverter_rated/Grid.Sbase) #Converting to p.u. value
           
            #Initialize inverter parameters according to DER rating
            self.initialize_inverter_parameters()
            #Initialize control loop gains according to DER rating
            self.initialize_controller_gains()
           
        else:
            raise ValueError('PV inverter parameters not available for DER with rating: ' + Sinverter_rated +' kVA')
    
    def initialize_inverter_parameters(self):
        """Initialize ratings and C, Lf, and Rf parameters."""
        _DER_rating = str(int(self.Sinverter_rated/1e3))
        
        self.Vdcrated = self.inverter_ratings[_DER_rating]['Vdcrated'] #Rated DC voltage
        if self.Vrms_rated is None:
            self.Varated = self.inverter_ratings[_DER_rating]['Varated'] #L-G peak to peak equivalent to 300 V L-L RMS
            self.Vrms_rated = self.Varated/math.sqrt(2)
        else:
            self.Varated = self.Vrms_rated*math.sqrt(2)
                
        self.a = Grid.Vgridrated/self.Varated  #Transformer turns ratio
        
        self.Vdcnominal = (self.Vdcrated/self.Vdcbase)   #Converting to p.u. value
        self.Vanominal = self.Varated/Grid.Vbase #Converting to p.u. value
        self.Vrms_ref =  self.Vanominal/math.sqrt(2) 
        
        if type(self).__name__ == 'SolarPV_DER_SinglePhase':
            self.Iarated = (self.Sinverter_rated/(self.Varated/math.sqrt(2)))*math.sqrt(2)
        elif type(self).__name__ == 'SolarPV_DER_ThreePhase':
            self.Iarated = (self.Sinverter_rated/(3*(self.Varated/math.sqrt(2))))*math.sqrt(2)
        
        self.iref_limit = (self.Iarated/self.Ibase)*self.Ioverload #Maximum current reference
        
        self.Rf_actual =  self.circuit_parameters[_DER_rating]['Rf_actual'] #Filter resistance
        self.Lf_actual = self.circuit_parameters[_DER_rating]['Lf_actual']  #Filter inductance
        self.C_actual =  self.circuit_parameters[_DER_rating]['C_actual']   #DC link capacitance
        self.Zf_actual =  self.Rf_actual + 1j*self.Lf_actual*Grid.wbase
        
        self.Rf = self.Rf_actual/Grid.Zbase        # VSC Filter resistance 
        self.Lf = self.Lf_actual/Grid.Lbase        # VSC Filter inductance
        self.Xf = (self.Lf_actual*Grid.wbase)/Grid.Zbase # VSC Filter reactance
        self.C = self.C_actual/Grid.Cbase          #DC link capacitor capacitance
        self.Zf = self.Zf_actual/Grid.Zbase
        
        #Interconnection to PCC - HV side
        #Actual values
        self.Z1_actual = self.circuit_parameters[_DER_rating]['Z1_actual']
        self.R1_actual = self.Z1_actual.real
        self.L1_actual = self.Z1_actual.imag/(2*math.pi*60.0)
        
        #Per-unit values
        self.R1 = self.R1_actual/Grid.Zbase  #Line/transformer resistance
        self.L1 = self.L1_actual/Grid.Lbase  #Line/transformer inductance
        self.Z1 =self.Z1_actual/Grid.Zbase    #Line/transformer impedance
        
        self.transformer_name = 'transformer_'+str(self.PV_DER_ID)
        
        self.check_PV_DER_parameters()  #Check PV-DER parameters
        
    def initialize_controller_gains(self):
        """Initialize controller settings."""
        #Current controller parameters
        _DER_rating = str(int(self.Sinverter_rated/1e3))
        self. Kp_GCC = 300/self.controller_parameters[_DER_rating]['scale_Kp_GCC'] #Current controller Proportional constant
        self.Ki_GCC = 70/self.controller_parameters[_DER_rating]['scale_Ki_GCC']  #Current controller Integral constant
        self.wp =  self.controller_parameters[_DER_rating]['wp']      #First order filter gain for GCC
        
        #Power (active and reactive) controller parameters
        self.Kp_DC = -1.0/self.controller_parameters[_DER_rating]['scale_Kp_DC']   #Active power controller Proportional constant
        self.Ki_DC = -0.5/self.controller_parameters[_DER_rating]['scale_Ki_DC'] #Active power controller Integral constant
        self.Kp_Q = 0.01/self.controller_parameters[_DER_rating]['scale_Kp_Q']  #Reactive power controller Proportional constant
        self.Ki_Q = 0.5/self.controller_parameters[_DER_rating]['scale_Ki_Q']   #Reactive power controller Integral constant    
    
    def attach_grid_model(self,grid_model):
        
        #Connect grid instance only if working in stand alone mode
        if self.standAlone and grid_model is not None:
            self.grid_model = grid_model
        elif self.standAlone and grid_model is None:
            raise ValueError('`Grid` instance need to be provided in stand alone mode for creating `SolarPV_DER` instance`!')
        else: #Grid model is not connected
            print('{self.name}:No grid model attached since PV-DER instance is not stand alone!')
                             
    def check_voltage(self):
                             
        if not self.standAlone and abs(abs(self.gridVoltagePhaseA) - self.Vanominal)/self.Vanominal > 0.1:
            raise ValueError('The rated PCC-LV voltage {} V has more than 10% deviation from the voltage input from external program {} V!'.format(self.Vanominal, abs(self.gridVoltagePhaseA)))
                             
        if self.m_steady_state*(self.Vdcrated/2) < self.Varated:
            raise ValueError('The nominal DC link voltage {:.1f} V is not sufficient for the inverter to generate the nominal voltage at PCC - LV side {:.1f} (L-G peak). Increase nominal DC link voltage to {:.1f} V.'.format(self.Vdcrated,self.Varated,math.ceil((self.Varated/self.m_steady_state)*2)))
            
    def check_PV_DER_parameters(self):
        """Method to check whether DER parameter's are feasible."""
        
        _del_I1max = 0.1*self.Iarated
        _Lf_min = self.Vdcrated/(16*self.fswitching*_del_I1max)
        _del_I1max_actual = self.Vdcrated/(16*self.fswitching*self.Lf_actual)
        if _del_I1max_actual > _del_I1max:   #Check if ripple current is less than 10 %
            print('Filter inductance {:.4} H is acceptable since AC side current ripple is {:.2}% (< 10%)'.format(_Lf_min,_del_I1max_actual/self.Iarated))
        else:
            print('Warning:Filter inductance {:.4} H results in AC side current ripple of {:.2}% (> 10%)'.format(_Lf_min,_del_I1max_actual/self.Iarated))
        
        _I_ripple = (0.25*self.Vdcrated)/(self.Lf_actual*self.fswitching)  #Maximum ripple voltage (p-p) at DC link
        _V_ripple = self.Vdcrated/(32*self.Lf_actual*self.C_actual*(self.fswitching**2))  #Maximum ripple voltage (p-p) at DC link
        _V_ripple_percentage = (_V_ripple/self.Vdcrated)*100
        if _V_ripple_percentage <= 1.0:   #Check if voltage ripple on DC link is less than 1%
            print('DC link capacitance of {:.4} F is acceptable since voltage ripple is only {:.2}% (< 1%)'.format(self.C_actual,_V_ripple_percentage))
        
        else:
            _V_ripple_ideal = self.Vdcrated*0.01  #1% ripple is acceptable
            _C = self.Vdcrated/(32*self.Lf_actual*_V_ripple_ideal*(self.fswitching**2))
            print('Warning:DC link capacitance of {:.4} F results in DC link voltage ripple of {:.2}% (> 1%)!Please use at least {} F.'.format(self.C_actual,_V_ripple_percentage,_C))
            #warnings.warn('Warning:DC link capacitance of {} F results in DC link voltage ripple of {:.3}% (> 1%)!Please use at least {} F.'.format(self.C_actual,_V_ripple_percentage,_C))               
     
    def power_error_calc(self,x):
        """Function for power."""
        maR = x[0]
        maI = x[1]
        iaR = x[2]
        iaI = x[3]
        
        ma = maR + 1j*maI
        
        ia = self.ia = iaR + 1j*iaI
        
        vta = ma*(self.Vdc/2)
        
        va = self.va_calc()
        
        St = (vta*ia.conjugate())/2
        S_PCC = (va*ia.conjugate())/2

        Ploss_filter = ((abs(ia)/math.sqrt(2))**2)*self.Rf 
        Qloss_filter = ((abs(ia)/math.sqrt(2))**2)*self.Xf
    
        if type(self).__name__ == 'SolarPV_DER_ThreePhase':
            mb = utility_functions.Ub_calc(ma)
            mc = utility_functions.Uc_calc(ma)
            
            ib = self.ib = utility_functions.Ub_calc(ia)
            ic = self.ic = utility_functions.Uc_calc(ia)
            
            vtb = mb*(self.Vdc/2)
            vtc = mc*(self.Vdc/2)
            
            vb = self.vb_calc()
            vc = self.vc_calc()
            
            St = St + (vtb*ib.conjugate() + vtc*ic.conjugate())/2
            S_PCC = S_PCC + (vb*ib.conjugate() + vc*ic.conjugate())/2
            
            Ploss_filter = Ploss_filter + ((abs(ib)/math.sqrt(2))**2)*self.Rf + ((abs(ic)/math.sqrt(2))**2)*self.Rf
            Qloss_filter = Qloss_filter + ((abs(ib)/math.sqrt(2))**2)*self.Xf + ((abs(ic)/math.sqrt(2))**2)*self.Xf    
    
        P_PCC_error = ((S_PCC.real + Ploss_filter)   - self.Ppv)**2 
        Q_PCC_error = (S_PCC.imag - self.Q_ref)**2   
        P_error = (St.real - self.Ppv)**2
        Q_error = (St.imag - Qloss_filter - self.Q_ref)**2
        return P_PCC_error  + Q_PCC_error + P_error# + Q_error
    
    def steady_state_calc(self):
        """Return steady state values."""
        
        #Find duty cycle that minimize steady state error
        print('Solving for steady state at current operating point.')
        x0 = np.array([0.89,0.0,124.0,3.59])
        result = minimize(self.power_error_calc, x0, method='nelder-mead',options={'xtol': 1e-8, 'disp': True})
        
        if result.success == False:
            raise ValueError('Steady state solution did not converge! Change operating point or disable steady state flag and try again.')
        
        ma0 = result.x[0] + 1j*result.x[1]
        self.ia = result.x[2] + 1j*result.x[3]
        self.ua = 0.0+0.0j
        self.xa = ma0
        self.vta = self.vta_calc()
        self.va = self.va_calc()
        
        self.xDC = self.ia.real
        self.xQ = self.ia.imag
        
        self.xPLL = 0.0
        self.wte = 2*math.pi
        
        if type(self).__name__ == 'SolarPV_DER_ThreePhase':
            mb0 = utility_functions.Ub_calc(ma0)
            mc0 = utility_functions.Uc_calc(ma0)
            
            self.xb = mb0
            self.xc = mc0
            
            self.ib = utility_functions.Ub_calc(self.ia)
            self.ic = utility_functions.Uc_calc(self.ia)
            
            self.ub = 0.0+0.0j
            self.uc = 0.0+0.0j            
            
            self.vtb = self.vtb_calc()
            self.vtc = self.vtc_calc()
            
            self.vb = self.vb_calc()
            self.vc = self.vc_calc()
        
        self.S =  self.S_calc()
        self.S_PCC = self.S_PCC_calc()
        self.Vtrms = self.Vtrms_calc()
        self.Vrms = self.Vrms_calc()
        self.Irms = self.Irms_calc()
        
        print('Steady state values for operating point defined by Ppv:{:.2f} W, Vdc:{:.2f} V, va:{:.2f} V found at:'.format(self.Ppv*self.Sbase,self.Vdc*self.Vdcbase,self.va*self.Vbase))
            
        self.show_PV_DER_states(quantity='power')
        self.show_PV_DER_states(quantity='duty_cycle')
        self.show_PV_DER_states(quantity='voltage')
        self.show_PV_DER_states(quantity='current')
        #print('Pt:{:.2f} W,Qt:{:.2f} VAR'.format(St0.real*self.Sbase,St0.imag*self.Sbase))
        #print('P_PCC:{:.2f} W,Q_PCC:{:.2f} VAR'.format(S_PCC0.real*self.Sbase,S_PCC0.imag*self.Sbase))
        #print('Voltage sum:{:.4f}'.format(vta0+vtb0+vtc0))
        #print('Current sum:{:.4f}'.format(ia0+ib0+ic0))
        """
        return [ia0,xa0,ua0,\
                ib0,xb0,ub0,\
                ic0,xc0,uc0,\
                self.Vdc,xDC0,xQ0,xPLL0,wte0]
        """