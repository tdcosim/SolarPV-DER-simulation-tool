from __future__ import division
import numpy as np
import math
import cmath
import scipy
import six
import pdb
import warnings

from scipy.optimize import fsolve, minimize
from pvder.grid_components import Grid
from pvder.utilities import SimulationUtilities
from pvder import utility_functions

class PV_module(object):
    """
    Class for describing PV module."
    """
    
    #Select either NN or polyfit model for MPP
    USE_POLYNOMIAL_MPP = True
    _MPP_fit_points = 50
    _Tactual_min = 273.15 #0 degree celsius
    _S_min = 10.0
    _Tactual_max = _Tactual_min + 35.0 #35 degree celsius
    _S_max = 100.0
    
    Iscr = 8.03 #Cell short-circuit current at reference temperature and radiation
    Kv = 0.0017  #Short-circuit current temperature co-efficient 
    T0 = 273.15 + 25.0 #Cell reference temperature in Kelvin
    Irs = 1.2e-7         #Cell reverse saturation current
    q = 1.602e-19        #Charge of an electron
    k = 1.38e-23        #Boltzmann's constant
    A = 1.92      #p-n junction ideality factor
    
    module_parameters = {'50':{'Np':11,'Ns':735,'Vdcmpp0':550.0,'Vdcmpp_min': 545.0,'Vdcmpp_max': 650.0},
                         '250':{'Np':45,'Ns':1000,'Vdcmpp0':750.0,'Vdcmpp_min': 650.0,'Vdcmpp_max': 1000.0}}
    
    def __init__(self,events,Sinverter_rated):
        #Events object
        self.events = events
        if Sinverter_rated in {50e3,100e3,250e3}:
           _DER_rating = str(int(Sinverter_rated/1e3))
           print('Creating PV module instance for DER with rating: ' + str(Sinverter_rated/1e3) + ' kVA')
           
           self.Np = self.module_parameters[str(_DER_rating)]['Np']
           self.Ns = self.module_parameters[str(_DER_rating)]['Ns']
           self.Vdcmpp0 = self.module_parameters[str(_DER_rating)]['Vdcmpp0']
           self.Vdcmpp_min = self.module_parameters[str(_DER_rating)]['Vdcmpp_min']
           self.Vdcmpp_max = self.module_parameters[str(_DER_rating)]['Vdcmpp_max']
           
        else:
           raise ValueError('PV module parameters not available for DER with rating: ' + str(Sinverter_rated/1e3)+' kVA')
        
        #Fit polynomial
        if self.USE_POLYNOMIAL_MPP:
           self.fit_MPP_poly()
        
        #PV conditions
        self.Sinsol,self.Tactual= events.solar_events(t=0.0)
        
        self.Iph = self.Iph_calc()       
        utility_functions.print_to_terminal('after PVCondition')
    
    @property    
    def Vdcmpp(self):
        """Voltage at maximum power point for given insolation and temperature"""
        
        #sol, = fsolve(lambda x: 2.5 - np.sqrt(x), 8)
        self.Iph = self.Iph_calc() #This function uses solar insolation
        
        return fsolve(lambda Vdc0:-((self.Np*self.Irs*(scipy.exp((self.q*Vdc0)/(self.k*self.Tactual*self.A*self.Ns))))*(self.q/(self.k*self.Tactual*self.A*self.Ns))*Vdc0)-((self.Np*self.Irs*(scipy.exp((self.q*Vdc0)/(self.k*self.Tactual*self.A*self.Ns))-1)))\
                       +(self.Np*self.Iph),self.Vdcmpp0)[0] #This is a time consuming operation 
    
    def Iph_calc(self):
        """Panel current for given insolation and temperature"""
        return (self.Iscr+(self.Kv*(self.Tactual-self.T0)))*(self.Sinsol/100.0)
    
    def Ppv_calc(self,Vdc_actual):
       """PV panel power output from  solar insolation."""
       self.Iph = self.Iph_calc()
       self.Ipv = (self.Np*self.Iph)-(self.Np*self.Irs*(math.exp((self.q*Vdc_actual)/(self.k*self.Tactual*self.A*self.Ns))-1))   #Faster  with Pure Python functions
       
       return max(0,(self.Ipv*Vdc_actual))/Grid.Sbase
       
       #return utility_functions.Ppv_calc(self.Iph,self.Np,self.Ns,Vdc_actual,self.Tactual,Grid.Sbase)
    
    def fit_MPP_poly(self):
        """Method to fit MPP to a polynomial function."""
        
        self.Tactual =  298.15  #Use constant temperature
        self._MPP_fit_points = 10
        _Srange = np.linspace(self._S_min,self._S_max,self._MPP_fit_points+1)
        _Vdcmpp_list = []
        _Ppvmpp_list =[]
        _Sinsol_list= []
        utility_functions.print_to_terminal('Calculating {} values for MPP polynomial fit!'.format(len(_Srange)))
        for S in _Srange:
            self.Sinsol = S
            _Vdcmpp = self.Vdcmpp
            _Ppvmpp = self.Ppv_calc(self.Vdcmpp)*Grid.Sbase
            
            _Sinsol_list.append(S)
            _Vdcmpp_list.append(_Vdcmpp)
            _Ppvmpp_list.append(_Ppvmpp)
            #print(f'S:{self.Sinsol},Tactual:{self.Tactual},Vdcmpp:{_Vdcmpp},Ppvmpp:{_Ppvmpp/1e3}')
        
        _x = np.array(_Sinsol_list)
        _y = np.array(_Vdcmpp_list)
        self.z = np.polyfit(_x, _y, 3)
        utility_functions.print_to_terminal('Found polynomial for MPP :{:.4f}x^3 + {:.4f}x^2 +{:.4f}x^1 + {:.4f}!'.format(self.z[0],self.z[1],self.z[2], self.z[3]))
    
class SolarPV_DER(PV_module,Grid):
    """
       Class for describing a Solar Photo-voltaic Distributed Energy Resource consisting of panel, converters, and
       control systems.
    """
    DER_count = 0
    #Number of ODE's
    n_ODE = 23
    
    Vdcbase = Grid.Vbase #DC side base value is same as AC side base value
    #PLL controller parameters
    Kp_PLL = 180 #1800
    Ki_PLL = 320 #32000
    
    #Inverter current overload rating (Max 10s)
    Ioverload = 1.5
    inverter_ratings = {'50':{'Varated':245.0,'Vdcrated':550.0},
                        '250':{'Varated':360.0,'Vdcrated':750.0}}
    
    circuit_parameters = {'50':{'Rf_actual':0.002,'Lf_actual' :25.0e-6,'C_actual':300.0e-6,'Z1_actual':0.0019 + 1j*0.0561},
                          '250':{'Rf_actual':0.002,'Lf_actual':300.0e-6,'C_actual':300.0e-6,'Z1_actual':0.0019 + 1j*0.0561}}
    
    controller_parameters = {'50':{'scale_Kp_GCC':0.05,'scale_Ki_GCC':0.05,\
                                   'scale_Kp_DC':0.05,'scale_Ki_DC' : 0.05,\
                                   'scale_Kp_Q' : 0.05,'scale_Ki_Q' : 0.05,'wp' : 20e4},
                             '250':{'scale_Kp_GCC':0.1,'scale_Ki_GCC':0.1,\
                                    'scale_Kp_DC':0.01,'scale_Ki_DC' : 0.01,\
                                    'scale_Kp_Q' : 0.01,'scale_Ki_Q' : 0.01,'wp' : 20e4}}
    #Frequency
    winv = we = 2.0*math.pi*60.0
    fswitching  = 10e3
    
     #Time delay before activating logic for MPP, Volt-VAR control,  LVRT/LFRT 
    t_stable = 1.0
    
    #Limits
    
    m_limit = 1.0 #Maximum duty cycle
    
    #Simulation time steps
    tStart = 0.0
    tStop = 0.5
    tInc = 0.001
    t = np.arange(tStart, tStop, tInc)
    
    #Duty cycle
    m_steady_state = 0.96 #Expected duty cycle at steady state    
    
    #Flags
    VOLT_VAR_ENABLE = False
    VOLT_VAR_FLAG = False
    
    VOLT_WATT_ENABLE = False
    VOLT_WATT_FLAG = False
    
    LVRT_ENABLE = False
    LFRT_ENABLE = False
    
    PRINT_INLINE = False
    VERBOSE = False
    MPPT_ENABLE = False
    
    DO_EXTRA_CALCULATIONS = False #Do calculations not essential to ODE model (useful for debugging)

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
    
    #Ramp control
    RAMP_ENABLE = False
    RAMP_FLAG = False
    ramp_list = []
    n_steps = 0
    ramp_del_t = 0.5 #0.025 
    
    def __init__(self,grid_model,events,\
                                       Sinverter_rated = 50.0e3,\
                                       ia0 = 0+0j,xa0 =0+0j , ua0 = 0+0j,\
                                       xDC0 = 0,xQ0 = 0,xPLL0 = 0.0,wte0 = 2*math.pi,\
                                       standAlone=True,gridVoltagePhaseA=.50+0j,\
                                       gridVoltagePhaseB=-.25-.43301270j,\
                                       gridVoltagePhaseC=-.25+.43301270j,STEADY_STATE_INITIALIZATION=False,\
                                       pvderConfig=None):   # Rl = 1.0,Ll = 0.0, ia0 = 0.26 -  1j*0.04,xa0 = 1.0 + 1j*0.15,xDC0 = 0.25,,xQ0 = -0.03
        
        ####
        self.standAlone,self.gridVoltagePhaseA,self.gridVoltagePhaseB,self.gridVoltagePhaseC=\
                                                                                             standAlone,\
                                                                                             gridVoltagePhaseA,gridVoltagePhaseB,gridVoltagePhaseC
        #Increment count to keep track of number of PV-DER model instances
        SolarPV_DER.DER_count = SolarPV_DER.DER_count+1
        self.PV_DER_ID = SolarPV_DER.DER_count
        
        if six.PY3:
            super().__init__(events,Sinverter_rated)  #Initialize PV module class (base class)
        elif six.PY2:
            super(SolarPV_DER,self).__init__(events,Sinverter_rated)
        
        self.STEADY_STATE_INITIALIZATION = STEADY_STATE_INITIALIZATION
        
        #Grid object
        self.grid_model = grid_model
        #Events object
        self.events = events
        #Object name
        self.name = 'PV_DER_'+str(self.PV_DER_ID)
        utility_functions.print_to_terminal('before Sinverter_rated')
        if Sinverter_rated in {50e3,100e3,250e3}:
           print('Creating PV inverter instance for DER with rating:' + str(Sinverter_rated/1e3) + ' kVA')
           self.Sinverter_rated = Sinverter_rated #Inverter rating in kVA
           self.Sinverter_nominal = (Sinverter_rated/Grid.Sbase) #Converting to p.u. value
           
           #Initialize inverter parameters according to DER rating
           self.initialize_inverter_parameters()
           #Initialize control loop gains according to DER rating
           self.initialize_controller_gains()
           
        else:
           raise ValueError('PV inverter parameters not available for DER with rating: ' + str(Sinverter_rated/1e3)+' kVA')
        
        if self.m_steady_state*(self.Vdcrated/2) < self.Varated:
            raise ValueError('The nominal DC link voltage {:.1f} V is not sufficient for the inverter to generate the nominal voltage at PCC - LV side {:.1f} (L-G peak). Increase nominal DC link voltage to {:.1f} V.'.format(self.Vdcrated,self.Varated,math.ceil((self.Varated/self.m_steady_state)*2)))    
        
        #Exhibit different behavior when used in standalone mode
        if self.standAlone:
            self.n_total_ODE = self.n_ODE + self.grid_model.n_grid_ODE
            self.Zload1_actual =  self.events.load_events(t=0.0)          #Load at PCC   
        
        else:
            self.n_total_ODE = self.n_ODE
            self.Zload1_actual =  10e6+1j*0.0     #Using a large value of impedance to represent a no load condition
        utility_functions.print_to_terminal('after Exhibit different behavior')
        
        self.Rload1_actual = self.Zload1_actual.real
        self.Lload1_actual = self.Z1_actual.imag/(2*math.pi*60.0)
        
        #Per-unit values
        self.Rload1 = self.Rload1_actual/self.grid_model.Zbase  #resistance in per unit value
        self.Lload1 = self.Lload1_actual/self.grid_model.Lbase  #inductance in per unit value
        self.Zload1 =self.Zload1_actual/self.grid_model.Zbase
        
        #Reference
        self.Q_ref = 0.0
        
        if self.MPPT_ENABLE:
            self.Vdc_ref = self.Vdcmpp/self.Vdcbase
            self.Vdc_ref_new = self.Vdcmpp/self.Vdcbase
        else:
            self.Vdc_ref = self.Vdcnominal
            self.Vdc_ref_new = self.Vdcnominal
        
        #DC link voltage
        self.Vdc = self.Vdc_ref
        #PV module power output
        self.Ppv = self.Ppv_calc(self.Vdc_actual)
        utility_functions.print_to_terminal('before Initialize all states with steady state')
        #Initialize all states with steady state values at current operating point
        if self.STEADY_STATE_INITIALIZATION == True:
           ia0,xa0,ua0,ib0,xb0,ub0,ic0,xc0,uc0,Vdc0,xDC0,xQ0,xPLL0,wte0 =  self.steady_state_calc()
        else:
           ib0 = utility_functions.Ub_calc(ia0)
           xb0 = utility_functions.Ub_calc(xa0)
           ub0 = utility_functions.Ub_calc(ua0)
           ic0 = utility_functions.Uc_calc(ia0)
           xc0 = utility_functions.Uc_calc(xa0)
           uc0 = utility_functions.Uc_calc(ua0)
        
        #Phase a
        self.ia = ia0
        self.xa = xa0
        self.ua = ua0
        
        #Phase b
        self.ib = ib0
        self.xb = xb0  #Shift by -120 degrees
        self.ub = ub0
        
        #Phase c
        self.ic = ic0
        self.xc = xc0   #Shift by +120 degrees
        self.uc = uc0
        
        #DC link voltage and reactive power controller
        self.xDC = xDC0
        self.xQ = xQ0
        
        #PLL
        self.xPLL = xPLL0
        self.wte = wte0
        utility_functions.print_to_terminal('before Derived voltages')
        #Derived voltages
        self.vta = self.vta_calc()
        self.vtb = self.vtb_calc()
        self.vtc = self.vtc_calc()
        self.va = self.va_calc()
        self.vb = self.vb_calc()
        self.vc = self.vc_calc()
        utility_functions.print_to_terminal('before Load current')
        #Load current
        self.iaload1 = self.iaload1_calc()
        self.ibload1 = self.ibload1_calc()
        self.icload1 = self.icload1_calc()
                
        #Derived powers
        utility_functions.print_to_terminal('before Derived powers')
        self.S = self.S_calc()
        self.S_PCC = self.S_PCC_calc()
        self.S_G = self.S_G_calc()
        self.S_load1 = self.S_load1_calc()
        self.S_PCCa = self.S_PCCa_calc()
        self.S_PCCb = self.S_PCCb_calc()
        self.S_PCCc = self.S_PCCc_calc()
        utility_functions.print_to_terminal('before RMS voltages')
        #RMS voltages
        self.Vrms = self.Vrms_calc()
        self.Vrms_ref =  self.Vanominal/math.sqrt(2) 
        if self.DO_EXTRA_CALCULATIONS:
            self.Vtrms = self.Vtrms_calc()
            self.Vtabrms = self.Vtabrms_calc()
            self.Vabrms = self.Vabrms_calc()
            #Inverter RMS current
            self.Irms = self.Irms_calc()
        utility_functions.print_to_terminal('before Reference currents')
        #Reference currents
        self.ia_ref = self.ia_ref_calc()
        self.ib_ref = self.ib_ref_calc()
        self.ic_ref = self.ic_ref_calc()
        utility_functions.print_to_terminal('before Convert to time domain')
        #Convert to time domain
        self.vat,self.vbt,self.vct = utility_functions.phasor_to_time(upha = self.va,uphb = self.vb,uphc = self.vc,w=self.grid_model.wgrid,t=0.0)
        utility_functions.print_to_terminal('before Convert from 3ph time domain')
        #Convert from 3ph time domain to d-q using Parks transformation
        self.vd,self.vq,self.v0 = utility_functions.abc_to_dq0(self.vat,self.vbt,self.vct,self.wte) #Grid voltage
        utility_functions.print_to_terminal('before LVRT settings')
        #LVRT settings
        self.LVRT_initialize(pvderConfig)
        
        if self.standAlone:
            self.vagt,self.vbgt,self.vcgt = utility_functions.phasor_to_time(upha = self.grid_model.vag,uphb = self.grid_model.vbg,uphc = self.grid_model.vcg,w=self.grid_model.wgrid,t=0.0)
            self.vdg,self.vqg,self.v0g = utility_functions.abc_to_dq0(self.vagt,self.vbgt,self.vcgt,self.wte)
        utility_functions.print_to_terminal('after standAlone')
        #PLL frequency
        self.we = self.we_calc()
        utility_functions.print_to_terminal('after standAlone')
    
    @property                         #Decorator used for auto updating
    def y0(self):
        """List of initial states"""
        return  [self.ia.real, self.ia.imag, self.xa.real, self.xa.imag, self.ua.real,self.ua.imag,\
                 self.ib.real, self.ib.imag, self.xb.real, self.xb.imag, self.ub.real,self.ub.imag,\
                 self.ic.real, self.ic.imag, self.xc.real, self.xc.imag, self.uc.real,self.uc.imag,\
                 self.Vdc,self.xDC,self.xQ,self.xPLL,self.wte]

    def initialize_inverter_parameters(self):
        """Initialize ratings and C, Lf, and Rf parameters."""
        _DER_rating = str(int(self.Sinverter_rated/1e3))
        
        self.Vdcrated = self.inverter_ratings[_DER_rating]['Vdcrated'] #Rated DC voltage
        self.Vdcnominal = (self.Vdcrated/self.Vdcbase)*1.0             #Converting to p.u. value
        
        self.Varated = self.inverter_ratings[_DER_rating]['Varated'] #L-G peak to peak equivalent to 300 V L-L RMS
        self.Vanominal = self.Varated/Grid.Vbase #Converting to p.u. value
        self.a = Grid.Vgridrated/self.Varated  #Transformer turns ratio
        
        self.Iarated = (self.Sinverter_rated/(3*(self.Varated/math.sqrt(2))))*math.sqrt(2)
        
        self.iref_limit = ((self.Sinverter_nominal/(3*(self.Vanominal/math.sqrt(2))))*math.sqrt(2))*self.Ioverload #Maximum current reference
        
        self.Rf_actual =  self.circuit_parameters[_DER_rating]['Rf_actual'] #Filter resistance
        self.Lf_actual = self.circuit_parameters[_DER_rating]['Lf_actual']  #Filter inductance
        self.C_actual =  self.circuit_parameters[_DER_rating]['C_actual']   #DC link capacitance
        self.Zf_actual =  self.Rf_actual + 1j*self.Lf_actual*Grid.wbase
        
        self.Rf = self.Rf_actual/Grid.Zbase        # VSC Filter resistance 
        self.Lf = self.Lf_actual/Grid.Lbase        # VSC Filter inductance
        self.Xf = (self.Lf_actual*Grid.wbase)/Grid.Zbase # VSC Filter reactance
        self.C = self.C_actual/Grid.Cbase          #DC link capacitor capacitance
        self.Zf = self.Zf_actual/Grid.Zbase
        
        #Interconnection to grid
        #Actual values
        self.Z1_actual = self.circuit_parameters[_DER_rating]['Z1_actual']
        self.R1_actual = self.Z1_actual.real
        self.L1_actual = self.Z1_actual.imag/(2*math.pi*60.0)
        
        #Per-unit values
        self.R1 = self.R1_actual/self.grid_model.Zbase  #Line/transformer resistance
        self.L1 = self.L1_actual/self.grid_model.Lbase  #Line/transformer inductance
        self.Z1 =self.Z1_actual/self.grid_model.Zbase    #Line/transformer impedance
        
        self.transformer_name = 'transformer_'+str(SolarPV_DER.DER_count)
        
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
    
    def show_PV_DER_states(self,quantity='voltage'):
        """Display values of voltage and current quantities."""
        if quantity not in {'voltage','current','power'}:
            raise ValueError('Unknown quantity: ' + str(plot_type))
        print('______{}_____\n'.format(self.name))
        if quantity ==  'voltage':
            print('Vdc:{:.2f}\nVta:{:.2f},Vtb:{:.2f},Vtb:{:.2f}\nVtrms:{:.2f}'.format(self.Vdc*self.Vbase,self.vta*self.Vbase,self.vtb*self.Vbase,self.vtc*self.Vbase,self.Vtrms*self.Vbase))
        elif quantity ==  'current':
            print('ia:{:.2f},ib:{:.2f},ic:{:.2f}\nIrms:{:.2f}'.format(self.ia*self.Ibase,self.ib*self.Ibase,self.ic*self.Ibase,self.Irms*self.Ibase))
        elif quantity ==  'power':
            print('Ppv:{:.2f}\nS:{:.2f}\nS_PCC:{:.2f}'.format(self.Ppv*self.Sbase,self.S*self.Sbase,self.S_PCC*self.Sbase))
    
    
    def show_PV_DER_parameters(self,quantity='voltages'):
        """Display rated values."""
        
        if quantity not in {'inverter_ratings','controller_gains','circuit_parameters'}:
            raise ValueError('Unknown quantity: ' + str(plot_type))
        
        if quantity ==  'inverter_ratings':
            print('Vdcrated:{:.3f}'.format(self.Vdcrated))
            print('Vtrated (L-G peak):{:.3f} V\nVrated (L-G peak):{:.3f} V'.format((self.Vdcrated/2)*self.m_steady_state,self.Varated))
        
        elif quantity == 'circuit_parameters':
            print('Cdc:{:.9f} F\nLf:{:.6f} H\nRf:{:.3f} Ohm'.format(self.C*self.Cbase,self.Lf*self.Lbase,self.Rf*self.Zbase))
        
        elif quantity == 'controller_parameters':
            print('Current controller:\nKp_GCC:{:.3f}, Ki_GCC:{:.3f}, wp:{:.3f}'.format(self.Kp_GCC,self.Ki_GCC,self.wp))
            print('DC link voltage controller:\nKp_DC:{:.3f}, Ki_DC:{:.3f}'.format(self.Kp_DC,self.Ki_DC))
            print('Reactive power controller:\nKp_Q:{:.3f}, Ki_Q:{:.3f}'.format(self.Kp_Q,self.Ki_Q))
            print('PLL controller:\nKp_PLL:{:.3f}, Ki_PLL:{:.3f}'.format(self.Kp_PLL,self.Ki_PLL))
    
    def check_PV_DER_parameters(self):
        """Method to check whether DER parameter's are feasible."""
        
        _del_I1max = 0.1*self.Iarated
        _Lf_min = self.Vdcrated/(16*self.fswitching*_del_I1max)
        _del_I1max_actual = self.Vdcrated/(16*self.fswitching*self.Lf_actual)
        if _del_I1max_actual > _del_I1max:   #Check if ripple current is less than 10 %
           print('Filter inductance {:.5} H is acceptable since AC side current ripple is {:.2}% (< 10%)'.format(_Lf_min,_del_I1max_actual/self.Iarated))
        else:
           print('Warning:Filter inductance {:.5} H results in AC side current ripple of {:.2}% (> 10%)'.format(_Lf_min,_del_I1max_actual/self.Iarated))
        
        _I_ripple = (0.25*self.Vdcrated)/(self.Lf_actual*self.fswitching)  #Maximum ripple voltage (p-p) at DC link
        _V_ripple = self.Vdcrated/(32*self.Lf_actual*self.C_actual*(self.fswitching**2))  #Maximum ripple voltage (p-p) at DC link
        _V_ripple_percentage = (_V_ripple/self.Vdcrated)*100
        if _V_ripple_percentage <= 1.0:   #Check if voltage ripple on DC link is less than 1%
           print('DC link capacitance of {} F is acceptable since voltage ripple is only {:.2}% (< 1%)'.format(self.C_actual,_V_ripple_percentage))
        
        else:
           _V_ripple_ideal = self.Vdcrated*0.01  #1% ripple is acceptable
           _C = self.Vdcrated/(32*self.Lf_actual*_V_ripple_ideal*(self.fswitching**2))
           #warnings.warn('Warning:DC link capacitance of {} F results in DC link voltage ripple of {:.3}% (> 1%)!Please use at least {} F.'.format(self.C_actual,_V_ripple_percentage,_C))
           print('Warning:DC link capacitance of {} F results in DC link voltage ripple of {:.3}% (> 1%)!Please use at least {} F.'.format(self.C_actual,_V_ripple_percentage,_C))
    
    def validate_model(self,PRINT_ERROR = True):
        """Compare error between RMS quantities and Phasor quantities."""
        
        #Calculation with phasor quantities
        self.Pf_phasor = self.S_calc().real-self.S_PCC_calc().real  #Active power consumed by filter resistor
        self.Qf_phasor = self.S_calc().imag-self.S_PCC_calc().imag  #Reactive power consumed by filter inductor 
        
        #Caculation with RMS quantities        
        self.Pf_RMS = 3*((self.Irms)**2)*self.Rf   #Active power consumed by filter resistor
        self.Qf_RMS = 3*((self.Irms)**2)*self.Xf   #Reactive power consumed by filter inductor 
        
        #Calculation with phasor quantities
        self.Pt_phasor = self.S_calc().real   #Active power output at inverter terminal
        self.Qt_phasor = self.S_calc().imag   #Reactive power output at inverter terminal
                
        #Caculation with RMS quantities 
        ra1,pha1 = cmath.polar(self.vta)
        ra2,pha2 = cmath.polar(self.ia)
        rb1,phb1 = cmath.polar(self.vtb)
        rb2,phb2 = cmath.polar(self.ib) 
        rc1,phc1 = cmath.polar(self.vtc)
        rc2,phc2 = cmath.polar(self.ic)
        
        self.Pt_RMS = (abs(self.vta)/math.sqrt(2))*(abs(self.ia)/math.sqrt(2))*math.cos(pha1-pha2) +\
                      (abs(self.vtb)/math.sqrt(2))*(abs(self.ib)/math.sqrt(2))*math.cos(phb1-phb2) +\
                      (abs(self.vtc)/math.sqrt(2))*(abs(self.ic)/math.sqrt(2))*math.cos(phc1-phc2)#Active power output at inverter terminal
        
        self.Qt_RMS = (abs(self.vta)/math.sqrt(2))*(abs(self.ia)/math.sqrt(2))*math.sin(pha1-pha2) +\
                      (abs(self.vtb)/math.sqrt(2))*(abs(self.ib)/math.sqrt(2))*math.sin(phb1-phb2) +\
                      (abs(self.vtc)/math.sqrt(2))*(abs(self.ic)/math.sqrt(2))*math.sin(phc1-phc2)#Reactive power output
        
        #self.Pt_RMS = 3*(self.Vtrms)*(self.Irms)*math.cos(ph1-ph2) #Active power output at inverter terminal
        #self.Qt_RMS = 3*(self.Vtrms)*(self.Irms)*math.sin(ph1-ph2) #Reactive power output at inverter terminal
        
        if PRINT_ERROR:
            print('Active power output error:{:.4f}\nReactive power output error:{:.4f}'.format(abs(self.Pt_phasor-self.Pt_RMS),abs(self.Qt_phasor-self.Qt_RMS)))    
            print('Inverter filter active power loss error:{:.4f}\nInverter filter reactive power loss error:{:.4f}'.format(abs(self.Pf_phasor-self.Pf_RMS),abs(self.Qf_phasor-self.Qf_RMS)))   
    
    #PLL equation (inverter frequency)
    def we_calc(self):
        """Calculate inverter frequency from PLL."""
        
        #return  self.Kp_PLL*(self.vdg) + self.xPLL + 2*math.pi*60.0
        return  self.Kp_PLL*(self.vd) + self.xPLL + 2*math.pi*60.0
    
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
    
    @property                         #Decorator used for auto updating
    def Vdc_actual(self):
        """Actual DC link voltage"""
        return min(self.Vdcrated,self.Vdc*self.Vdcbase)  #Calculate actual voltage
    
    #Apparent power output at inverter terminal
    def S_calc(self):
        """Inverter apparent power output"""
        
        return (1/2)*(self.vta*self.ia.conjugate() + self.vtb*self.ib.conjugate() + self.vtc*self.ic.conjugate())*1.0
        #return utility_functions.S_calc(self.vta,self.vtb,self.vtc,self.ia,self.ib,self.ic)
        
    #Apparent power output at PCC - LV side
    def S_PCC_calc(self):
        """Power output at PCC LV side"""
        return (1/2)*(self.va*self.ia.conjugate() + self.vb*self.ib.conjugate() + self.vc*self.ic.conjugate())*1.0
        #return utility_functions.S_calc(self.va,self.vb,self.vc,self.ia,self.ib,self.ic)
        
    def S_load1_calc(self):
        """Power absorbed by load at PCC LV side."""
        
        return (1/2)*(self.va*(-(self.va/self.Zload1)).conjugate() + self.vb*(-(self.vb/self.Zload1)).conjugate() + self.vc*(-(self.vc/self.Zload1)).conjugate())
    
    def S_G_calc(self):
        """Power absorbed/produced by grid voltage source."""
    
        return (1/2)*((-(self.ia-(self.va/self.Zload1))/self.a).conjugate()*self.grid_model.vag+(-(self.ib-(self.vb/self.Zload1))/self.a).conjugate()*self.grid_model.vbg+(-(self.ic-(self.vc/self.Zload1))/self.a).conjugate()*self.grid_model.vcg)
    
    def S_PCCa_calc(self):
        """Inverter apparent power output - phase a"""
        
        return (1/2)*(self.va*self.ia.conjugate())*1.0
    
    def S_PCCb_calc(self):
        """Inverter apparent power output - phase b"""
        
        return (1/2)*(self.vb*self.ib.conjugate())*1.0
    
    def S_PCCc_calc(self):
        """Inverter apparent power output - phase c"""
        
        return (1/2)*(self.vc*self.ic.conjugate())*1.0	
    
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
        
            val = ((self.grid_model.vbg+(self.ib/self.a)*self.grid_model.Z2)/(self.a) +self.ib*self.Z1)*((self.Zload1*self.a*self.a)/((self.a*self.a*(self.Z1+self.Zload1))+self.grid_model.Z2))
        else:
            val=self.gridVoltagePhaseB
        
        return val
    
    #@property
    def vc_calc(self):
        """PCC - LV side - Phase C"""
        
        if self.standAlone:
            
            val = ((self.grid_model.vcg+(self.ic/self.a)*self.grid_model.Z2)/(self.a) +self.ic*self.Z1)*((self.Zload1*self.a*self.a)/((self.a*self.a*(self.Z1+self.Zload1))+self.grid_model.Z2))
        else:
            val=self.gridVoltagePhaseC
        
        return val
    def iaload1_calc(self):
        """Current counsumed by load connected at PCC LV side -  Phase A."""
        
        return self.va/self.Zload1
    
    def ibload1_calc(self):
        """Current counsumed by load connected at PCC LV side -  Phase B."""
        
        return self.vb/self.Zload1
    
    def icload1_calc(self):
        """Current counsumed by load connected at PCC LV side -  Phase C."""
        
        return self.vc/self.Zload1
    #@property
    def Vtrms_calc(self):
        """Inverter terminal voltage -  RMS"""
        
        return utility_functions.Urms_calc(self.vta,self.vtb,self.vtc)
    
    def Vrms_calc(self):
        """PCC LV side voltage - RMS"""
        
        return utility_functions.Urms_calc(self.va,self.vb,self.vc)
    
    def Irms_calc(self):
        """Inverter current - RMS"""
        
        return utility_functions.Urms_calc(self.ia,self.ib,self.ic)
        
    def Vtabrms_calc(self):
        """Inverter terminal voltage - line to line  RMS"""
        
        return abs(self.vta-self.vtb)/math.sqrt(2)
    
    def Vabrms_calc(self):
        """PCC LV side voltage - line to line  RMS"""
        
        return abs(self.va-self.vb)/math.sqrt(2)
    
    def MPP_table(self):
        """Method to output Vdc reference corresponding to MPP at different insolation levels values."""
        
        if self.USE_POLYNOMIAL_MPP == True:
             _Vdcmpp = np.polyval(self.z , self.Sinsol)
        else:
            _Vdcmpp = self.Vdcrated

        _Vdcmpp=  max(min(_Vdcmpp,self.Vdcmpp_max),self.Vdcmpp_min)
        
        return _Vdcmpp/self.Vdcbase
    
    def update_inverter_states(self,ia,xa,ua,ib,xb,ub,ic,xc,uc,Vdc,xDC,xQ,xPLL,wte):
        """Update inverter states"""
        
        self.ia = ia
        self.xa = xa
        self.ua = ua
        
        self.ib = ib
        self.xb = xb
        self.ub = ub
        
        self.ic = ic
        self.xc = xc
        self.uc = uc
        
        self.Vdc = Vdc
        self.xDC = xDC
        self.xQ = xQ
        
        self.xPLL = xPLL
        self.wte = wte
        
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
    
    def ODE_model(self,y,t):
        """Derivatives for the equation."""
        
        iaR, iaI, xaR, xaI, uaR, uaI,\
        ibR, ibI, xbR, xbI, ubR, ubI,\
        icR, icI, xcR, xcI, ucR, ucI,\
        Vdc, xDC, xQ, xPLL, wte = y   # unpack current values of y
        
        self.update_inverter_states(iaR + 1j*iaI, xaR + 1j*xaI,uaR + 1j*uaI,\
                                    ibR + 1j*ibI, xbR + 1j*xbI,ubR + 1j*ubI,\
                                    icR + 1j*icI, xcR + 1j*xcI,ucR + 1j*ucI,\
                                    Vdc,xDC,xQ,\
                                    xPLL,wte)
                
        #self.xb = utility_functions.Ub_calc(self.xa)
        #self.xc = utility_functions.Uc_calc(self.xa)
        
        #PV conditions
        Sinsol_new,Tactual_new = self.events.solar_events(t)
        if abs(self.Sinsol- Sinsol_new) or abs(self.Tactual- Tactual_new) > 0.0:
             
            self.Sinsol = Sinsol_new
            self.Tactual = Tactual_new
            utility_functions.print_to_terminal("{}:PV module current output changed from {:.3f} A to {:.3f} A at {:.3f}".format(self.name,self.Iph,self.Iph_calc(),t))
            self.Iph = self.Iph_calc()
        
        self.Ppv = self.Ppv_calc(self.Vdc_actual)
        
        #Load at PCC LV side
        if self.standAlone:
            Zload1_actual_new =  self.events.load_events(t)
            Zload1_new = Zload1_actual_new/self.grid_model.Zbase
        
            if abs(self.Zload1- Zload1_new)> 0.0:
               self.Zload1 = Zload1_new
               utility_functions.print_to_terminal("Load at PCC LV side changed from {:.3f} VA to {:.3f} VA at {:.3f}".format(self.S_load1,self.S_load1_calc(),t))
        else:
            self.Zload1 =  self.Zload1    #Load at PCC-LV side is always constant
        
        #Update inverter terminal voltage
        self.vta = self.vta_calc()
        self.vtb = self.vtb_calc()
        self.vtc = self.vtc_calc()
        
        #Update PCC LV side voltage
        self.va = self.va_calc()
        self.vb = self.vb_calc()
        self.vc = self.vc_calc()
        self.Vrms = self.Vrms_calc()
        
        #Update power output
        self.S = self.S_calc()
        self.S_PCC = self.S_PCC_calc()
        
        #Update RMS values
        if self.DO_EXTRA_CALCULATIONS:
            self.Vtrms = self.Vtrms_calc()
            self.Vtabrms = self.Vtabrms_calc()        
            self.Vabrms = self.Vabrms_calc()        
            self.Irms = self.Irms_calc()
            self.S_PCCa = self.S_PCCa_calc()
            self.S_PCCb = self.S_PCCb_calc()
            self.S_PCCc = self.S_PCCc_calc()
        
        #Update load current in stand alone mode
        if self.standAlone:
            self.iaload1 = self.iaload1_calc()
            self.ibload1 = self.ibload1_calc()
            self.icload1 = self.icload1_calc()
            
            self.S_G = self.S_G_calc()
            self.S_load1 = self.S_load1_calc()
        
        #Get reactive power set-point
        if self.VOLT_VAR_ENABLE == True:
            self.Q_ref = self.Volt_VAR_logic(t)
        else:
            self.Q_ref = 0.0 
        
        #Get DC link voltage set point
        if self.MPPT_ENABLE == True:
            self.Vdc_ref_new = self.MPP_table()
        else:
            self.Vdc_ref_new = self.Vdcnominal
        
        if self.RAMP_ENABLE:
            self.Vdc_ramp(t)
        
        #Get current controller setpoint
        self.ia_ref = self.ia_ref_calc()
        self.ib_ref = self.ib_ref_calc()
        self.ic_ref = self.ic_ref_calc()
        
        #d-q transformation
        
        #Convert PCC LV side voltage from phasor to time domain
        self.vat,self.vbt,self.vct = utility_functions.phasor_to_time(upha = self.va,uphb = self.vb,uphc = self.vc,w=self.grid_model.wgrid,t=t)
        
        #Convert from 3ph time domain to d-q using Parks transformation
        self.vd,self.vq,self.v0 = utility_functions.abc_to_dq0(self.vat,self.vbt,self.vct,self.wte) #PCC LV side voltage
        
        #Convert grid voltage from phasor to time domain to d-q in stand alone mode
        if self.standAlone:
           self.vagt,self.vbgt,self.vcgt = utility_functions.phasor_to_time(upha = self.grid_model.vag,uphb = self.grid_model.vbg,uphc = self.grid_model.vcg,w=self.grid_model.wgrid,t=t)
           self.vdg,self.vqg,self.v0g = utility_functions.abc_to_dq0(self.vagt,self.vbgt,self.vcgt,self.wte) #Grid voltage
        
        #Calculate inverter frequency from PLL equation
        self.we = self.we_calc()
        self.winv = self.we

        if self.LVRT_ENABLE == True:
            self.LVRT(t)
            if self.LVRT_TRIP == True and self.LVRT_RECONNECT == False:
               self.PV_DER_disconnect()
        #LFRT trip logic
        if self.LFRT_ENABLE == True:
            self.FRT(t)
            if self.LFRT_TRIP == True and self.LFRT_RECONNECT == False:
               self.PV_DER_disconnect() 
        
        #Phase a inverter output current
        diaR = (1/self.Lf)*(-self.Rf*self.ia.real - self.va.real + self.vta.real) + (self.winv/self.wbase)*self.ia.imag 
        diaI = (1/self.Lf)*(-self.Rf*self.ia.imag - self.va.imag + self.vta.imag) - (self.winv/self.wbase)*self.ia.real  
       
        #Current controller dynamics
        if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit*1e1:
            if np.sign(self.Ki_GCC*self.ua.real) == np.sign(self.xa.real):
                dxaR = 0.0
            else:
                dxaR = self.Ki_GCC*self.ua.real
            if np.sign(self.Ki_GCC*self.ua.imag) == np.sign(self.xa.imag):
                dxaI = 0.0
            else:
                dxaI = self.Ki_GCC*self.ua.imag
                #six.print_(dxaR+1j*dxaI,np.sign(self.Ki_GCC*self.ua))
        else:
            dxaR = self.Ki_GCC*self.ua.real
            dxaI = self.Ki_GCC*self.ua.imag
            
        if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit*1e1:
            if np.sign( (self.wp)*(-self.ua.real +  self.ia_ref.real - self.ia.real)) == np.sign(self.ua.real):
                duaR = 0.0
            else:
                duaR = (self.wp)*(-self.ua.real +  self.ia_ref.real - self.ia.real)
                
            if np.sign((self.wp)*(-self.ua.imag +  self.ia_ref.imag - self.ia.imag)) == np.sign(self.ua.imag):
                duaI = 0.0
            else:
                duaI = (self.wp)*(-self.ua.imag +  self.ia_ref.imag - self.ia.imag)
                
        else:
            duaR = (self.wp)*(-self.ua.real +  self.ia_ref.real - self.ia.real)
            duaI = (self.wp)*(-self.ua.imag +  self.ia_ref.imag - self.ia.imag)
        
        #Phase b inverter output current
        dibR = (1/self.Lf)*(-self.Rf*self.ib.real - self.vb.real + self.vtb.real) + (self.winv/self.wbase)*self.ib.imag 
        dibI = (1/self.Lf)*(-self.Rf*self.ib.imag - self.vb.imag + self.vtb.imag) - (self.winv/self.wbase)*self.ib.real  
        
        #Current controller dynamics - Phase b
        #dxbR = self.Ki_GCC*self.ub.real
        #dxbI = self.Ki_GCC*self.ub.imag
        #dubR = (self.wp)*(-self.ub.real +  self.ib_ref.real - self.ib.real)
        #dubI = (self.wp)*(-self.ub.imag +  self.ib_ref.imag - self.ib.imag)
        
        if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
            if np.sign(self.Ki_GCC*self.ub.real) == np.sign(self.xb.real):
                dxbR = 0.0
            else:
                dxbR = self.Ki_GCC*self.ub.real
            if np.sign(self.Ki_GCC*self.ub.imag) == np.sign(self.xb.imag):
                dxbI = 0.0
            else:
                dxbI = self.Ki_GCC*self.ub.imag
                #print(dxbR+1j*dxbI,np.sign(self.Ki_GCC*self.ub))
        else: 
            dxbR = self.Ki_GCC*self.ub.real
            dxbI = self.Ki_GCC*self.ub.imag
        
        if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
            if np.sign( (self.wp)*(-self.ub.real +  self.ib_ref.real - self.ib.real)) == np.sign(self.ub.real):
                dubR = 0.0
            else:
                dubR = (self.wp)*(-self.ub.real +  self.ib_ref.real - self.ib.real)
                
            if np.sign((self.wp)*(-self.ub.imag +  self.ib_ref.imag - self.ib.imag)) == np.sign(self.ub.imag):
                dubI = 0.0
            else:
                dubI = (self.wp)*(-self.ub.imag +  self.ib_ref.imag - self.ib.imag)
        else:
            dubR = (self.wp)*(-self.ub.real +  self.ib_ref.real - self.ib.real)
            dubI = (self.wp)*(-self.ub.imag +  self.ib_ref.imag - self.ib.imag)
        
        #Phase c inverter output current
        dicR = (1/self.Lf)*(-self.Rf*self.ic.real - self.vc.real + self.vtc.real) + (self.winv/self.wbase)*self.ic.imag 
        dicI = (1/self.Lf)*(-self.Rf*self.ic.imag - self.vc.imag + self.vtc.imag) - (self.winv/self.wbase)*self.ic.real 
        
        #Current controller dynamics - Phase c
        #dxcR = self.Ki_GCC*self.uc.real
        #dxcI = self.Ki_GCC*self.uc.imag
        #ducR = (self.wp)*(-self.uc.real +  self.ic_ref.real - self.ic.real)
        #ducI = (self.wp)*(-self.uc.imag +  self.ic_ref.imag - self.ic.imag)
        
        if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
            if np.sign(self.Ki_GCC*self.uc.real) == np.sign(self.xc.real):
                dxcR = 0.0
            else:
                dxcR = self.Ki_GCC*self.uc.real
                
            if np.sign(self.Ki_GCC*self.uc.imag) == np.sign(self.xc.imag):
                dxcI = 0.0
            else:
                dxcI = self.Ki_GCC*self.uc.imag
                #print(dxaR+1j*dxaI,np.sign(self.Ki_GCC*self.ua))
        else: 
            dxcR = self.Ki_GCC*self.uc.real
            dxcI = self.Ki_GCC*self.uc.imag
        
        if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
            if np.sign( (self.wp)*(-self.uc.real +  self.ic_ref.real - self.ic.real)) == np.sign(self.uc.real):
                ducR = 0.0
            else:
                ducR = (self.wp)*(-self.uc.real +  self.ic_ref.real - self.ic.real)
            
            if np.sign((self.wp)*(-self.uc.imag +  self.ic_ref.imag - self.ic.imag)) == np.sign(self.uc.imag):
                ducI = 0.0
            else:
                ducI = (self.wp)*(-self.uc.imag +  self.ic_ref.imag - self.ic.imag)
        else:
            ducR = (self.wp)*(-self.uc.real +  self.ic_ref.real - self.ic.real)
            ducI = (self.wp)*(-self.uc.imag +  self.ic_ref.imag - self.ic.imag)
        
        #DC link voltage dynamics
        dVdc = (self.Ppv - self.S.real)/(self.Vdc*self.C)
        
        if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ  - self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
            if np.sign(self.Ki_DC*(self.Vdc_ref - self.Vdc)) == np.sign(self.xDC):
                dxDC = 0.0
            else:
                dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
        else:
            dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
            
        # Reactive power controller dynamics
        if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ  - self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
            
            if np.sign(-self.Ki_Q*(self.Q_ref - self.S_PCC.imag)) == np.sign(self.xQ):
                dxQ = 0.0
            else:
                dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
        else:
            dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
        
        #SRF-PLL dynamics
        #dxPLL = self.Ki_PLL*(self.vdg)
        dxPLL = self.Ki_PLL*(self.vd)
        
        #Frequency integration to get angle
        dwte = self.we
        
        result =     [ diaR,# list of dy/dt=f functions
                       diaI,
                       dxaR,
                       dxaI,
                       duaR,
                       duaI,
                       dibR,
                       dibI,
                       dxbR,
                       dxbI,
                       dubR,
                       dubI,
                       dicR,
                       dicI,
                       dxcR,
                       dxcI,
                       ducR,
                       ducI,
                       dVdc,
                       dxDC,
                       dxQ,
                       dxPLL,
                       dwte]
        
        return result

    def jac_ODE_model(self,y,t):
        """Jacobian for the system of ODE's."""
        
        iaR, iaI, xaR, xaI, uaR, uaI,\
        ibR, ibI, xbR, xbI, ubR, ubI,\
        icR, icI, xcR, xcI, ucR, ucI,\
        Vdc, xDC, xQ, xPLL, wte = y   # unpack current values of y
        
        self.update_inverter_states(iaR + 1j*iaI, xaR + 1j*xaI,uaR + 1j*uaI,\
                                    ibR + 1j*ibI, xbR + 1j*xbI,ubR + 1j*ubI,\
                                    icR + 1j*icI, xcR + 1j*xcI,ucR + 1j*ucI,\
                                    Vdc,xDC,xQ,\
                                    xPLL,wte)
        
        #self.xb = utility_functions.Ub_calc(self.xa)
        #self.xc = utility_functions.Uc_calc(self.xa)
        
        J=np.zeros((self.n_total_ODE,self.n_total_ODE))
        varInd={}; n=0
        for entry in ['iaR','iaI','xaR','xaI','uaR','uaI','ibR','ibI','xbR','xbI','ubR','ubI','icR','icI',\
            'xcR','xcI','ucR','ucI','Vdc','xDC','xQ','xPLL','wte']:
            varInd[entry]=n
            n+=1
        
        #PV conditions
        Sinsol_new,Tactual_new = self.events.solar_events(t)
        if abs(self.Sinsol- Sinsol_new) or abs(self.Tactual- Tactual_new) > 0.0:
             
            self.Sinsol = Sinsol_new
            self.Tactual = Tactual_new
            utility_functions.print_to_terminal("PV panel current output changed from {:.3f} A to {:.3f} A at {:.3f}".format(self.Iph,self.Iph_calc(),t))
            self.Iph = self.Iph_calc()
        
        self.Ppv = self.Ppv_calc(self.Vdc_actual)
        
        """
        #Load at PCC LV side
        if self.standAlone:
            Zload1_actual_new =  self.events.load_events(t)
            Zload1_new = Zload1_actual_new/self.grid_model.Zbase
        
            if abs(self.Zload1- Zload1_new)> 0.0:
               self.Zload1 = Zload1_new
               utility_functions.print_to_terminal("Load at PCC LV side changed from {:.3f} VA to {:.3f} VA at {:.3f}".format(self.S_load1,self.S_load1_calc(),t))
        else:
            self.Zload1 =  self.Zload1    #Load at PCC-LV side is always constant
        """
        
        #Update inverter terminal voltage
        self.vta = self.vta_calc()
        self.vtb = self.vtb_calc()
        self.vtc = self.vtc_calc()
        
        #Update PCC LV side voltage
        self.va = self.va_calc()
        self.vb = self.vb_calc()
        self.vc = self.vc_calc()
        self.Vrms = self.Vrms_calc()
        
        #Update RMS values
        if self.DO_EXTRA_CALCULATIONS:
            self.Vtrms = self.Vtrms_calc()
            self.Vtabrms = self.Vtabrms_calc()        
            self.Vabrms = self.Vabrms_calc()        
            self.Irms = self.Irms_calc()        
        
        #Update power output
        self.S = self.S_calc()
        self.S_PCC = self.S_PCC_calc()
        
        if self.standAlone:
            #Update load current in stand alone mode
            self.iaload1 = self.iaload1_calc()
            self.ibload1 = self.ibload1_calc()
            self.icload1 = self.icload1_calc()
            
            self.S_G = self.S_G_calc()
            self.S_load1 = self.S_load1_calc()
        
        #Get reactive power set-point
        if self.VOLT_VAR_ENABLE == True:
            self.Q_ref = self.Volt_VAR_logic(t)
        else:
            self.Q_ref = 0.0 
        
        #Get DC link voltage set point
        if self.MPPT_ENABLE == True:
            self.Vdc_ref_new = self.MPP_table()
        else:
            self.Vdc_ref_new = self.Vdcnominal
        
        if self.RAMP_ENABLE:
            self.Vdc_ramp(t)
        
        #Get current controller setpoint
        self.ia_ref = self.ia_ref_calc()
        self.ib_ref = self.ib_ref_calc()
        self.ic_ref = self.ic_ref_calc()
        
        #d-q transformation
        
        #Convert PCC LV side voltage from phasor to time domain
        self.vat,self.vbt,self.vct = utility_functions.phasor_to_time(upha = self.va,uphb = self.vb,uphc = self.vc,w=self.grid_model.wgrid,t=t)
        
        #Convert from 3ph time domain to d-q using Parks transformation
        self.vd,self.vq,self.v0 = utility_functions.abc_to_dq0(self.vat,self.vbt,self.vct,self.wte) #PCC LV side voltage
        
        #Calculate inverter frequency from PLL equation
        self.we = self.we_calc()
        self.winv = self.we
       
        #LVRT trip logic
        if self.LVRT_TRIP == True and self.LVRT_RECONNECT == False:
               self.PV_DER_disconnect()
        #LFRT trip logic
        if self.LFRT_TRIP == True and self.LFRT_RECONNECT == False:
               self.PV_DER_disconnect() 
        
        #Phase a inverter output current
        
        ra,theta_a = cmath.polar(self.va)
        rb,theta_b = cmath.polar(self.vb)
        rc,theta_c = cmath.polar(self.vc)
        theta_a = self.grid_model.wgrid*t + theta_a - math.pi/2
        theta_b = self.grid_model.wgrid*t + theta_b - math.pi/2
        theta_c = self.grid_model.wgrid*t + theta_c - math.pi/2
            
        J[varInd['iaR'],varInd['iaR']] = -self.Rf/self.Lf            
        J[varInd['iaR'],varInd['iaI']] = (self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
        J[varInd['iaR'],varInd['xaR']] = self.Vdc/(2*self.Lf)
        J[varInd['iaR'],varInd['uaR']] = (self.Vdc*self.Kp_GCC)/(2*self.Lf)
        J[varInd['iaR'],varInd['Vdc']] = (self.xa.real+self.ua.real*self.Kp_GCC)/(2*self.Lf)
        J[varInd['iaR'],varInd['xPLL']] = self.ia.imag/self.wbase
        J[varInd['iaR'],varInd['wte']] = (2/3)*((self.Kp_PLL*self.ia.imag)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte)
                                         + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
                                         + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
        
        J[varInd['iaI'],varInd['iaR']]= -(self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
        J[varInd['iaI'],varInd['iaI']]= -self.Rf/self.Lf
        J[varInd['iaI'],varInd['xaI']]= self.Vdc/(2*self.Lf) 
        J[varInd['iaI'],varInd['uaI']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
        J[varInd['iaI'],varInd['Vdc']]= (self.xa.imag+self.ua.imag*self.Kp_GCC)/(2*self.Lf)
        J[varInd['iaI'],varInd['xPLL']]= -self.ia.real/self.wbase
        J[varInd['iaI'],varInd['wte']] = -(2/3)*((self.Kp_PLL*self.ia.real)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte) 
                                         + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
                                         + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
            
        #Current controller dynamics
        if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit*1e1:
            if np.sign(self.Ki_GCC*self.ua.real) == np.sign(self.xa.real):
                J[varInd['xaR'],varInd['uaR']]=0.0
            else:
                J[varInd['xaR'],varInd['uaR']]=self.Ki_GCC
            if np.sign(self.Ki_GCC*self.ua.imag) == np.sign(self.xa.imag):
                J[varInd['xaI'],varInd['uaI']]=0.0
            else:
                J[varInd['xaI'],varInd['uaI']]=self.Ki_GCC
                
        else:
                J[varInd['xaR'],varInd['uaR']]=self.Ki_GCC
                J[varInd['xaI'],varInd['uaI']]=self.Ki_GCC
        
        if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit*1e1:
            if np.sign( (self.wp)*(-self.ua.real +  self.ia_ref.real - self.ia.real)) == np.sign(self.ua.real):
                J[varInd['uaR'],varInd['iaR']]= 0.0
                J[varInd['uaR'],varInd['uaR']]= 0.0
                J[varInd['uaR'],varInd['Vdc']]= 0.0
                J[varInd['uaR'],varInd['xDC']]= 0.0   
            else:
                #duaR = (self.wp)*(-self.ua.real +  self.ia_ref.real - self.ia.real)
                J[varInd['uaR'],varInd['iaR']]= -self.wp
                J[varInd['uaR'],varInd['uaR']]= -self.wp
                J[varInd['uaR'],varInd['Vdc']]= -self.wp*self.Kp_DC
                J[varInd['uaR'],varInd['xDC']]= self.wp                    
                    
            if np.sign((self.wp)*(-self.ua.imag +  self.ia_ref.imag - self.ia.imag)) == np.sign(self.ua.imag):
                #duaI = 0.0
                J[varInd['uaI'],varInd['iaR']]= 0.0
                J[varInd['uaI'],varInd['iaI']]= 0.0
                J[varInd['uaI'],varInd['uaI']]= 0.0
                J[varInd['uaI'],varInd['xQ']]= 0.0
                    
                J[varInd['uaI'],varInd['ibR']]= 0.0
                J[varInd['uaI'],varInd['ibI']]= 0.0
                J[varInd['uaI'],varInd['icR']]= 0.0
                J[varInd['uaI'],varInd['icI']]= 0.0
            else:
                #duaI = (self.wp)*(-self.ua.imag +  self.ia_ref.imag - self.ia.imag)
                J[varInd['uaI'],varInd['iaR']]= (self.Kp_Q*self.wp*self.va.imag/2)
                J[varInd['uaI'],varInd['iaI']]= -self.wp - (self.Kp_Q*self.wp*self.va.real/2)
                J[varInd['uaI'],varInd['uaI']]= -self.wp
                J[varInd['uaI'],varInd['xQ']]= self.wp
                    
                J[varInd['uaI'],varInd['ibR']]= (self.Kp_Q*self.wp*self.vb.imag/2)
                J[varInd['uaI'],varInd['ibI']]= - (self.Kp_Q*self.wp*self.vb.real/2)
                J[varInd['uaI'],varInd['icR']]= (self.Kp_Q*self.wp*self.vc.imag/2)
                J[varInd['uaI'],varInd['icI']]= - (self.Kp_Q*self.wp*self.vc.real/2)
                
        else:
            #duaR = (self.wp)*(-self.ua.real +  self.ia_ref.real - self.ia.real)
            #duaI = (self.wp)*(-self.ua.imag +  self.ia_ref.imag - self.ia.imag)
            J[varInd['uaR'],varInd['iaR']]= -self.wp
            J[varInd['uaR'],varInd['uaR']]= -self.wp
            J[varInd['uaR'],varInd['Vdc']]= -self.wp*self.Kp_DC
            J[varInd['uaR'],varInd['xDC']]= self.wp   
                
            J[varInd['uaI'],varInd['iaR']]= (self.Kp_Q*self.wp*self.va.imag/2)
            J[varInd['uaI'],varInd['iaI']]= -self.wp - (self.Kp_Q*self.wp*self.va.real/2)
            J[varInd['uaI'],varInd['uaI']]= -self.wp
            J[varInd['uaI'],varInd['xQ']]= self.wp
                
            J[varInd['uaI'],varInd['ibR']]= (self.Kp_Q*self.wp*self.vb.imag/2)
            J[varInd['uaI'],varInd['ibI']]= - (self.Kp_Q*self.wp*self.vb.real/2)
            J[varInd['uaI'],varInd['icR']]= (self.Kp_Q*self.wp*self.vc.imag/2)
            J[varInd['uaI'],varInd['icI']]= - (self.Kp_Q*self.wp*self.vc.real/2)
        
        #Phase b inverter output current
        #dibR = (1/self.Lf)*(-self.Rf*self.ib.real - self.vb.real + self.vtb.real) + (self.winv/self.wbase)*self.ib.imag 
        #dibI = (1/self.Lf)*(-self.Rf*self.ib.imag - self.vb.imag + self.vtb.imag) - (self.winv/self.wbase)*self.ib.real  
            
        J[varInd['ibR'],varInd['ibR']]= -self.Rf/self.Lf            
        J[varInd['ibR'],varInd['ibI']]= (self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
        J[varInd['ibR'],varInd['xbR']]= self.Vdc/(2*self.Lf)
        J[varInd['ibR'],varInd['ubR']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
        J[varInd['ibR'],varInd['Vdc']]= (self.xb.real+self.ub.real*self.Kp_GCC)/(2*self.Lf)
        J[varInd['ibR'],varInd['xPLL']]= self.ib.imag/self.wbase     
        J[varInd['ibR'],varInd['wte']] = (2/3)*((self.Kp_PLL*self.ib.imag)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte)
                                         + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
                                         + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
            
        J[varInd['ibI'],varInd['ibR']]= -(self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
        J[varInd['ibI'],varInd['ibI']]= -self.Rf/self.Lf
        J[varInd['ibI'],varInd['xbI']]= self.Vdc/(2*self.Lf)
        J[varInd['ibI'],varInd['ubI']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
        J[varInd['ibI'],varInd['Vdc']]= (self.xb.imag+self.ub.imag*self.Kp_GCC)/(2*self.Lf)
        J[varInd['ibI'],varInd['xPLL']]= -self.ib.real/self.wbase
        J[varInd['ibI'],varInd['wte']] = -(2/3)*((self.Kp_PLL*self.ib.real)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte) 
                                         + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
                                         + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
        
        #Current controller dynamics - Phase b
        
        if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
            if np.sign(self.Ki_GCC*self.ub.real) == np.sign(self.xb.real):
                #dxbR = 0.0
                J[varInd['xbR'],varInd['ubR']]=0.0
            else:
                #dxbR = self.Ki_GCC*self.ub.real
                J[varInd['xbR'],varInd['ubR']]=self.Ki_GCC
            if np.sign(self.Ki_GCC*self.ub.imag) == np.sign(self.xb.imag):
                #dxbI = 0.0
                J[varInd['xbI'],varInd['ubI']]=0.0
            else:
                #dxbI = self.Ki_GCC*self.ub.imag
                J[varInd['xbI'],varInd['ubI']]=self.Ki_GCC

        else: 
            #dxbR = self.Ki_GCC*self.ub.real
            #dxbI = self.Ki_GCC*self.ub.imag
            J[varInd['xbR'],varInd['ubR']]=self.Ki_GCC
            J[varInd['xbI'],varInd['ubI']]=self.Ki_GCC
        
        if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
            if np.sign( (self.wp)*(-self.ub.real +  self.ib_ref.real - self.ib.real)) == np.sign(self.ub.real):
                #dubR = 0.0
                J[varInd['ubR'],varInd['ibR']]= 0.0
                J[varInd['ubR'],varInd['ubR']]= 0.0
                J[varInd['ubR'],varInd['Vdc']]= 0.0
                J[varInd['ubR'],varInd['xDC']]= 0.0
            else:
                #dubR = (self.wp)*(-self.ub.real +  self.ib_ref.real - self.ib.real)
                J[varInd['ubR'],varInd['ibR']]= -self.wp
                J[varInd['ubR'],varInd['ubR']]= -self.wp
                J[varInd['ubR'],varInd['Vdc']]= -self.wp*self.Kp_DC
                J[varInd['ubR'],varInd['xDC']]= self.wp  
                    
            if np.sign((self.wp)*(-self.ub.imag +  self.ib_ref.imag - self.ib.imag)) == np.sign(self.ub.imag):
                #dubI = 0.0
                J[varInd['ubI'],varInd['ibR']]= 0.0
                J[varInd['ubI'],varInd['ibI']]= 0.0
                J[varInd['ubI'],varInd['ubI']]= 0.0
                J[varInd['ubI'],varInd['xQ']]= 0.0
                    
                J[varInd['ubI'],varInd['iaR']]= 0.0
                J[varInd['ubI'],varInd['iaI']]= 0.0
                J[varInd['ubI'],varInd['icR']]= 0.0
                J[varInd['ubI'],varInd['icI']]= 0.0
                
            else:
                #dubI = (self.wp)*(-self.ub.imag +  self.ib_ref.imag - self.ib.imag)
                J[varInd['ubI'],varInd['ibR']]= (self.Kp_Q*self.wp*self.vb.imag/2)
                J[varInd['ubI'],varInd['ibI']]= -self.wp - (self.Kp_Q*self.wp*self.vb.real/2)
                J[varInd['ubI'],varInd['ubI']]= -self.wp
                J[varInd['ubI'],varInd['xQ']]= self.wp
                    
                J[varInd['ubI'],varInd['iaR']]= (self.Kp_Q*self.wp*self.va.imag/2)
                J[varInd['ubI'],varInd['iaI']]= - (self.Kp_Q*self.wp*self.va.real/2)
                J[varInd['ubI'],varInd['icR']]= (self.Kp_Q*self.wp*self.vc.imag/2)
                J[varInd['ubI'],varInd['icI']]= - (self.Kp_Q*self.wp*self.vc.real/2)
                    
        else:
            #dubR = (self.wp)*(-self.ub.real +  self.ib_ref.real - self.ib.real)
            #dubI = (self.wp)*(-self.ub.imag +  self.ib_ref.imag - self.ib.imag)
            J[varInd['ubR'],varInd['ibR']]= -self.wp
            J[varInd['ubR'],varInd['ubR']]= -self.wp
            J[varInd['ubR'],varInd['Vdc']]= -self.wp*self.Kp_DC
            J[varInd['ubR'],varInd['xDC']]= self.wp
                    
            J[varInd['ubI'],varInd['ibR']]= (self.Kp_Q*self.wp*self.vb.imag/2)
            J[varInd['ubI'],varInd['ibI']]= -self.wp - (self.Kp_Q*self.wp*self.vb.real/2)
            J[varInd['ubI'],varInd['ubI']]= -self.wp
            J[varInd['ubI'],varInd['xQ']]= self.wp
                    
            J[varInd['ubI'],varInd['iaR']]= (self.Kp_Q*self.wp*self.va.imag/2)
            J[varInd['ubI'],varInd['iaI']]= - (self.Kp_Q*self.wp*self.va.real/2)
            J[varInd['ubI'],varInd['icR']]= (self.Kp_Q*self.wp*self.vc.imag/2)
            J[varInd['ubI'],varInd['icI']]= - (self.Kp_Q*self.wp*self.vc.real/2)
        
        #Phase c inverter output current
        #dicR = (1/self.Lf)*(-self.Rf*self.ic.real - self.vc.real + self.vtc.real) + (self.winv/self.wbase)*self.ic.imag 
        #dicI = (1/self.Lf)*(-self.Rf*self.ic.imag - self.vc.imag + self.vtc.imag) - (self.winv/self.wbase)*self.ic.real 
        
        J[varInd['icR'],varInd['icR']]= -self.Rf/self.Lf            
        J[varInd['icR'],varInd['icI']]= (self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
        J[varInd['icR'],varInd['xcR']]= self.Vdc/(2*self.Lf)
        J[varInd['icR'],varInd['ucR']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
        J[varInd['icR'],varInd['Vdc']]= (self.xc.real+self.uc.real*self.Kp_GCC)/(2*self.Lf)
        J[varInd['icR'],varInd['xPLL']]= self.ic.imag/self.wbase 
        J[varInd['icR'],varInd['wte']] = (2/3)*((self.Kp_PLL*self.ic.imag)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte)
                                          + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
                                          + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
            
        J[varInd['icI'],varInd['icR']]= -(self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
        J[varInd['icI'],varInd['icI']]= -self.Rf/self.Lf
        J[varInd['icI'],varInd['xcI']]= self.Vdc/(2*self.Lf)
        J[varInd['icI'],varInd['ucI']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
        J[varInd['icI'],varInd['Vdc']]= (self.xc.imag+self.uc.imag*self.Kp_GCC)/(2*self.Lf)
        J[varInd['icI'],varInd['xPLL']]= -self.ic.real/self.wbase
        J[varInd['icI'],varInd['wte']] = -(2/3)*((self.Kp_PLL*self.ic.real)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte) 
                                         + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
                                         + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
            
        #Current controller dynamics - Phase c
        
        if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
            if np.sign(self.Ki_GCC*self.uc.real) == np.sign(self.xc.real):
                #dxcR = 0.0
                J[varInd['xcR'],varInd['ucR']]=0.0
            else:
                #dxcR = self.Ki_GCC*self.uc.real
                J[varInd['xcR'],varInd['ucR']]=self.Ki_GCC
            if np.sign(self.Ki_GCC*self.uc.imag) == np.sign(self.xc.imag):
                #dxcI = 0.0
                J[varInd['xcI'],varInd['ucI']]=0.0
            else:
                #dxcI = self.Ki_GCC*self.uc.imag
                J[varInd['xcI'],varInd['ucI']]=self.Ki_GCC
        else: 
            #dxcR = self.Ki_GCC*self.uc.real
            #dxcI = self.Ki_GCC*self.uc.imag
            J[varInd['xcR'],varInd['ucR']]=self.Ki_GCC
            J[varInd['xcI'],varInd['ucI']]=self.Ki_GCC
        
        if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
            if np.sign( (self.wp)*(-self.uc.real +  self.ic_ref.real - self.ic.real)) == np.sign(self.uc.real):
                #ducR = 0.0
                J[varInd['ucR'],varInd['icR']]= 0.0
                J[varInd['ucR'],varInd['ucR']]= 0.0
                J[varInd['ucR'],varInd['Vdc']]= 0.0
                J[varInd['ucR'],varInd['xDC']]= 0.0
            else:
                ducR = (self.wp)*(-self.uc.real +  self.ic_ref.real - self.ic.real)
                J[varInd['ucR'],varInd['icR']]= -self.wp
                J[varInd['ucR'],varInd['ucR']]= -self.wp
                J[varInd['ucR'],varInd['Vdc']]= -self.wp*self.Kp_DC
                J[varInd['ucR'],varInd['xDC']]= self.wp
                    
            if np.sign((self.wp)*(-self.uc.imag +  self.ic_ref.imag - self.ic.imag)) == np.sign(self.uc.imag):
                #ducI = 0.0
                J[varInd['ucI'],varInd['icR']]= 0.0
                J[varInd['ucI'],varInd['icI']]= 0.0
                J[varInd['ucI'],varInd['ucI']]= 0.0
                    
                J[varInd['ucI'],varInd['iaR']]= 0.0
                J[varInd['ucI'],varInd['iaI']]= 0.0
                J[varInd['ucI'],varInd['ibR']]= 0.0
                J[varInd['ucI'],varInd['ibI']]= 0.0
                    
                J[varInd['ucI'],varInd['Vdc']]= 0.0
                J[varInd['ucI'],varInd['xDC']]= 0.0
                J[varInd['ucI'],varInd['xQ']]= 0.0
                    
            else:
                #ducI = (self.wp)*(-self.uc.imag +  self.ic_ref.imag - self.ic.imag)
                J[varInd['ucI'],varInd['icR']]= -(self.Kp_Q*self.wp*self.vc.imag/4)
                J[varInd['ucI'],varInd['icI']]= -self.wp + (self.Kp_Q*self.wp*self.vc.real/4)
                J[varInd['ucI'],varInd['ucI']]= -self.wp
                    
                J[varInd['ucI'],varInd['iaR']]= -self.Kp_Q*self.wp*self.va.imag/4
                J[varInd['ucI'],varInd['iaI']]= self.Kp_Q*self.wp*self.va.real/4
                J[varInd['ucI'],varInd['ibR']]= -self.Kp_Q*self.wp*self.vb.imag/4
                J[varInd['ucI'],varInd['ibI']]= self.Kp_Q*self.wp*self.vb.real/4
                    
                J[varInd['ucI'],varInd['Vdc']]= -0.8660254037*self.Kp_DC*self.wp
                J[varInd['ucI'],varInd['xDC']]= 0.8660254037*self.wp
                J[varInd['ucI'],varInd['xQ']]= -0.5*self.wp                    
        else:
            #ducR = (self.wp)*(-self.uc.real +  self.ic_ref.real - self.ic.real)
            #ducI = (self.wp)*(-self.uc.imag +  self.ic_ref.imag - self.ic.imag)
            J[varInd['ucR'],varInd['icR']]= -self.wp
            J[varInd['ucR'],varInd['ucR']]= -self.wp
            J[varInd['ucR'],varInd['Vdc']]= -self.wp*self.Kp_DC
            J[varInd['ucR'],varInd['xDC']]= self.wp
                
            J[varInd['ucI'],varInd['icR']]= -(self.Kp_Q*self.wp*self.vc.imag/4)
            J[varInd['ucI'],varInd['icI']]= -self.wp + (self.Kp_Q*self.wp*self.vc.real/4)
            J[varInd['ucI'],varInd['ucI']]= -self.wp
                
            J[varInd['ucI'],varInd['iaR']]= -self.Kp_Q*self.wp*self.va.imag/4
            J[varInd['ucI'],varInd['iaI']]= self.Kp_Q*self.wp*self.va.real/4
            J[varInd['ucI'],varInd['ibR']]= -self.Kp_Q*self.wp*self.vb.imag/4
            J[varInd['ucI'],varInd['ibI']]= self.Kp_Q*self.wp*self.vb.real/4
                
            J[varInd['ucI'],varInd['Vdc']]= -0.8660254037*self.Kp_DC*self.wp
            J[varInd['ucI'],varInd['xDC']]= 0.8660254037*self.wp
            J[varInd['ucI'],varInd['xQ']]= -0.5*self.wp
        
        #DC link voltage dynamics
        dVdc = (self.Ppv - self.S.real)/(self.Vdc*self.C)
        J[varInd['Vdc'],varInd['iaR']]= -(self.xa.real+self.Kp_GCC*self.ua.real)/(4*self.C)
        J[varInd['Vdc'],varInd['iaI']]= -(self.xa.imag+self.Kp_GCC*self.ua.imag)/(4*self.C)
        J[varInd['Vdc'],varInd['xaR']]= -self.ia.real/(4*self.C)
        J[varInd['Vdc'],varInd['xaI']]= -self.ia.imag/(4*self.C)
        J[varInd['Vdc'],varInd['uaR']]= -(self.Kp_GCC*self.ia.real)/(4*self.C)
        J[varInd['Vdc'],varInd['uaI']]= -(self.Kp_GCC*self.ia.imag)/(4*self.C)

        J[varInd['Vdc'],varInd['ibR']]= -(self.xb.real+self.Kp_GCC*self.ub.real)/(4*self.C)
        J[varInd['Vdc'],varInd['ibI']]= -(self.xb.imag+self.Kp_GCC*self.ub.imag)/(4*self.C)
        J[varInd['Vdc'],varInd['xbR']]= -self.ib.real/(4*self.C)
        J[varInd['Vdc'],varInd['xbI']]= -self.ib.imag/(4*self.C)
        J[varInd['Vdc'],varInd['ubR']]= -(self.Kp_GCC*self.ib.real)/(4*self.C)
        J[varInd['Vdc'],varInd['ubI']]= -(self.Kp_GCC*self.ib.imag)/(4*self.C)

        J[varInd['Vdc'],varInd['icR']]= -(self.xc.real+self.Kp_GCC*self.uc.real)/(4*self.C)
        J[varInd['Vdc'],varInd['icI']]= -(self.xc.imag+self.Kp_GCC*self.uc.imag)/(4*self.C)
        J[varInd['Vdc'],varInd['xcR']]= -self.ic.real/(4*self.C)
        J[varInd['Vdc'],varInd['xcI']]= -self.ic.imag/(4*self.C)
        J[varInd['Vdc'],varInd['ucR']]= -(self.Kp_GCC*self.ic.real)/(4*self.C)
        J[varInd['Vdc'],varInd['ucI']]= -(self.Kp_GCC*self.ic.imag)/(4*self.C)                
            
        J[varInd['Vdc'],varInd['Vdc']]= (-(self.q*self.Np*self.Irs*(self.Vdcbase**2))/(self.C*self.k*self.A*self.Ns*self.Tactual*self.Sbase))*math.exp((self.q*self.Vdc*self.Vdcbase)/(self.k*self.A*self.Ns*self.Tactual))
        #DC link voltage controller dynamics
        if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ  - self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
            if np.sign(self.Ki_DC*(self.Vdc_ref - self.Vdc)) == np.sign(self.xDC):
                #dxDC = 0.0
                J[varInd['xDC'],varInd['Vdc']]= 0.0
            else:
                #dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
                J[varInd['xDC'],varInd['Vdc']]=-self.Ki_DC
        else:
            #dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
            J[varInd['xDC'],varInd['Vdc']]=-self.Ki_DC
            
        # Reactive power controller dynamics
        if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ  - self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
            
            if np.sign(-self.Ki_Q*(self.Q_ref - self.S_PCC.imag)) == np.sign(self.xQ):
                #dxQ = 0.0
                J[varInd['xQ'],varInd['iaR']]= 0.0
                J[varInd['xQ'],varInd['iaI']]= 0.0
                J[varInd['xQ'],varInd['ibR']]= 0.0
                J[varInd['xQ'],varInd['ibI']]= 0.0
                J[varInd['xQ'],varInd['icR']]= 0.0
                J[varInd['xQ'],varInd['icI']]= 0.0
            else:
                #dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
                J[varInd['xQ'],varInd['iaR']]= (self.Ki_Q*self.va.imag/2)
                J[varInd['xQ'],varInd['iaI']]= -(self.Ki_Q*self.va.real/2)
                J[varInd['xQ'],varInd['ibR']]= (self.Ki_Q*self.vb.imag/2)
                J[varInd['xQ'],varInd['ibI']]= -(self.Ki_Q*self.vb.real/2)
                J[varInd['xQ'],varInd['icR']]= (self.Ki_Q*self.vc.imag/2)
                J[varInd['xQ'],varInd['icI']]= -(self.Ki_Q*self.vc.real/2)
        
        else:
            #dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
            J[varInd['xQ'],varInd['iaR']]= (self.Ki_Q*self.va.imag/2)
            J[varInd['xQ'],varInd['iaI']]= -(self.Ki_Q*self.va.real/2)
            J[varInd['xQ'],varInd['ibR']]= (self.Ki_Q*self.vb.imag/2)
            J[varInd['xQ'],varInd['ibI']]= -(self.Ki_Q*self.vb.real/2)
            J[varInd['xQ'],varInd['icR']]= (self.Ki_Q*self.vc.imag/2)
            J[varInd['xQ'],varInd['icI']]= -(self.Ki_Q*self.vc.real/2)
        
        #SRF-PLL dynamics
        #dxPLL = self.Ki_PLL*(self.vd)
        J[varInd['xPLL'],varInd['wte']] = (2/3)*self.Ki_PLL*(-ra*math.cos(theta_a)*math.sin(self.wte) +
                                             0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)+
                                             0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
        
        #Frequency integration to get angle
        #dwte = self.we
        J[varInd['wte'],varInd['xPLL']]= 1
        J[varInd['wte'],varInd['wte']] = (2/3)*self.Kp_PLL*(-ra*math.cos(theta_a)*math.sin(self.wte) +
                                          0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)+
                                          0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
        
        return J
    
    def power_error_calc(self,x):
        """Function for power."""
        maR = x[0]
        maI = x[1]
    
        ma = maR + 1j*maI
        mb = utility_functions.Ub_calc(ma)
        mc = utility_functions.Uc_calc(ma)
    
        vta = ma*(self.Vdc/2)
        vtb = mb*(self.Vdc/2)
        vtc = mc*(self.Vdc/2)
    
        ia = (vta - self.gridVoltagePhaseA)/(self.Rf + 1j*(self.winv/self.wbase)*self.Lf)
        ib = (vtb - self.gridVoltagePhaseB)/(self.Rf + 1j*(self.winv/self.wbase)*self.Lf)
        ic = (vtc - self.gridVoltagePhaseC)/(self.Rf + 1j*(self.winv/self.wbase)*self.Lf)
    
        St = (vta*ia.conjugate() + vtb*ib.conjugate() + vtc*ic.conjugate())/2
        S_PCC = (self.gridVoltagePhaseA*ia.conjugate() + self.gridVoltagePhaseB*ib.conjugate() + self.gridVoltagePhaseC*ic.conjugate())/2

        Ploss_filter = ((abs(ia)/math.sqrt(2))**2)*self.Rf + ((abs(ib)/math.sqrt(2))**2)*self.Rf + ((abs(ic)/math.sqrt(2))**2)*self.Rf
    
        P_PCC_error = ((S_PCC.real + Ploss_filter)   - self.Ppv)**2
        Q_PCC_error = (S_PCC.imag - self.Q_ref)**2   
    
        return P_PCC_error  + Q_PCC_error
    
    def steady_state_calc(self):
        """Return steady state values."""
        
        #Find duty cycle that minimize steady state error
        print('Solving for steady state at current operating point.')
        x0 = np.array([0.89,0.0])
        result = minimize(self.power_error_calc, x0, method='nelder-mead',options={'xtol': 1e-8, 'disp': True})
        
        if result.success == False:
           raise ValueError('Steady state solution did not converge! Change operating point or disable steady state flag and try again.')
        
        ma0 = result.x[0] + 1j*result.x[1]
        mb0 = utility_functions.Ub_calc(ma0)
        mc0 = utility_functions.Uc_calc(ma0)
        
        ua0 = 0.0+0.0j
        ub0 = 0.0+0.0j
        uc0 = 0.0+0.0j
        
        xa0 = ma0
        xb0 = mb0
        xc0 = mc0
        
        vta0 = ma0*(self.Vdc/2)
        vtb0 = mb0*(self.Vdc/2)
        vtc0 = mc0*(self.Vdc/2)
        
        ia0 = (vta0 - self.gridVoltagePhaseA)/(self.Rf + 1j*(self.winv/self.wbase)*self.Lf)
        ib0 = (vtb0 - self.gridVoltagePhaseB)/(self.Rf + 1j*(self.winv/self.wbase)*self.Lf)
        ic0 = (vtc0 - self.gridVoltagePhaseC)/(self.Rf + 1j*(self.winv/self.wbase)*self.Lf)
        
        xDC0 = ia0.real
        xQ0 = ia0.imag
        
        xPLL0 = 0.0
        wte0 = 2*math.pi
        
        St0 = (vta0*ia0.conjugate() + vtb0*ib0.conjugate() + vtc0*ic0.conjugate())/2
        S_PCC0 = (self.gridVoltagePhaseA*ia0.conjugate() + self.gridVoltagePhaseB*ib0.conjugate() + self.gridVoltagePhaseC*ic0.conjugate())/2
        
        print('Steady state values for operating point defined by Ppv:{:.2f} W, Vdc:{:.2f} V, va:{:.2f} V found at:'.format(self.Ppv*self.Sbase,self.Vdc*self.Vdcbase,self.gridVoltagePhaseA*self.Vbase))
        print('Pt:{:.2f} W,Qt:{:.2f} VAR'.format(St0.real*self.Sbase,St0.imag*self.Sbase))
        print('P_PCC:{:.2f} W,Q_PCC:{:.2f} VAR'.format(S_PCC0.real*self.Sbase,S_PCC0.imag*self.Sbase))
        print('ma:{:.2f}'.format(ma0))
        print('vta:{:.2f}'.format(vta0*self.Vbase))
        print('ia:{:.2f}'.format(ia0*self.Ibase))
        print('Voltage sum:{:.4f}'.format(vta0+vtb0+vtc0))
        
        return [ia0,xa0,ua0,\
                ib0,xb0,ub0,\
                ic0,xc0,uc0,\
                self.Vdc,xDC0,xQ0,xPLL0,wte0]
    
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
        
        self.t_reconnect_limit = self.t_LV2_limit  #Time lag before reconnecting
        
        if pvderConfig['LVRT_INSTANTANEOUS_TRIP']:
            self.t_LV0_limit = self.t_LV1_limit = self.t_LV1_limit = 1/60 #Disconnect within one cycle

        assert (self.V_LV1 < self.V_LV2 and self.t_LV1_limit < self.t_LV2_limit) == True, "Voltage level 2 should be greater than Voltage level 1"

        #V1 to V2 - zone 2,V1 < - zone 1
    
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
