"""Code for initializing and validating PV-DER model instances."""

from __future__ import division
import logging

import math
import cmath
import numpy as np
from scipy.optimize import fsolve, minimize

from pvder.utility_classes import Logging
from pvder.grid_components import BaseValues
from pvder import utility_functions
from pvder import defaults,templates,specifications

PHASE_DIFFERENCE_120 = 120.0*(math.pi/180.0)

class PVDER_SetupUtilities(BaseValues,Logging):
    """
       Utility class for error checking during model initialization.
    """
    module_parameters = {}
    inverter_ratings = {}
    circuit_parameters = {}
    controller_gains ={}
    steadystate_values= {}    
    
    solver_spec = specifications.steadystate_solver_spec
    steadystate_solver = defaults.STEADYSTATE_SOLVER
    
    def creation_message(self):
        """Message after PV-DER instance was created."""        
       
        self.logger.info('{}:Instance created with DER parameter ID: {}; Specifications - Srated:{:.1f} kVA, Ppv:{:.1f} kW, Vrms:{:.1f} V, Steady state:{},LVRT Enable:{},HVRT Enable:{}'.format(self.name,self.parameter_ID,self.Sinverter_rated/1e3,(self.Ppv*BaseValues.Sbase)/1e3,self.Vrms_rated,self.steady_state_initialization,self.LVRT_ENABLE,self.HVRT_ENABLE))
        
    def update_DER_config(self,DER_config,DER_arguments,DER_id):
        """Update PV-DER design."""  
        
        self.check_DER_config(DER_config,DER_id)
        
        #parameters_from_arguments = {}
        #parameters_from_config = {}
        #parameters_from_default = {}
        
        for DER_component in self.DER_design_template:
            if DER_component in DER_arguments['derConfig']: 
                self.DER_config[DER_component] = DER_arguments['derConfig'][DER_component]
            else:
                for DER_parameter in self.DER_design_template[DER_component]:
                    _ = self.update_DER_parameter(DER_config,DER_arguments,DER_id,DER_component,DER_parameter)

        self.basic_specs = {DER_id:self.DER_config['basic_specs']}
        self.module_parameters = {DER_id:self.DER_config['module_parameters']}
        self.inverter_ratings = {DER_id:self.DER_config['inverter_ratings']}
        self.circuit_parameters = {DER_id:self.DER_config['circuit_parameters']}
        self.circuit_parameters[DER_id].update({'Z1_actual':self.DER_config['circuit_parameters']['R1_actual'] + 1j*self.DER_config['circuit_parameters']['X1_actual']})
        
        self.controller_gains = {DER_id:self.DER_config['controller_gains']}
        self.steadystate_values = {DER_id:self.DER_config['steadystate_values']}   
        
    def update_DER_parameter(self,DER_config,DER_arguments,DER_id,DER_component,DER_parameter):
        """Update DER config using both config file and arguments."""
        
        parameters ={}
        
        if DER_parameter in DER_arguments: #Check if parameter exists in DER key word arguments
            self.logger.debug('{}:Parameter {} in {} found in arguments with value {}.'.format(self.name,DER_parameter,DER_component,DER_arguments[DER_parameter])) 
            self.DER_config[DER_component].update({DER_parameter:DER_arguments[DER_parameter]})
            parameters.update({DER_parameter:DER_arguments[DER_parameter]})
            source = 'DER arguments'
            
        elif DER_parameter in DER_config[DER_component]: #Check if parameter exists in config file
            if isinstance(DER_config[DER_component][DER_parameter],(int,float)):
                self.logger.debug('{}:Parameter {} in {} found for ID:{} - with value {}.'.format(self.name,DER_parameter,DER_component,DER_id,DER_config[DER_component][DER_parameter])) 
                self.DER_config[DER_component].update({DER_parameter:DER_config[DER_component][DER_parameter]})
                parameters.update({DER_parameter:DER_config[DER_component][DER_parameter]})
                source = 'DER config'
            else:
                raise ValueError('Found {} to have type {} - expected type:(int,float)!'.format(DER_parameter,type(DER_config[DER_component][DER_parameter])))
        
        elif DER_parameter in self.DER_design_template[DER_component]:  #Check if parameter exists in default DER config
            self.logger.debug('{}:Parameter {} in {} not found for ID:{} - updating with default value {}.'.format(self.name,DER_parameter,DER_component,DER_id,self.DER_design_template[DER_component][DER_parameter])) 
            self.DER_config[DER_component].update({DER_parameter:self.DER_design_template[DER_component][DER_parameter]})
            parameters.update({DER_parameter:self.DER_design_template[DER_component][DER_parameter]})
            source = 'DER default values'
        else:
            raise ValueError('{}: Parameter {} in {} not found for ID:{} - update config file and try again.'.format(self.name,DER_parameter,DER_component,DER_id))
        
        #self.logger.debug('{}:Following parameters were obtained from {}:{}'.format(self.name,source,parameters))
        
        return parameters
             
    def check_DER_config(self,DER_config,DER_id):
        """Check DER config."""
        
        available_keys = sorted(list(DER_config.keys()))
        required_keys = sorted(self.DER_design_template.keys())
        missing_keys = set(required_keys)-set(available_keys)
        if not bool(missing_keys):
            self.logger.debug('{}:DER configuration with ID {} contained all required keys.'.format(self.name,DER_id))        
        else:
            raise KeyError('{}:DER configuration with ID:{} did not contain the required config keys:{}!'.format(self.name, DER_id,missing_keys))              
    
    def check_config_file(self,config_file,config_dict):
        """Check whether DER config file contains necessary fields."""
        
        self.logger.debug('Checking config file {}.'.format(config_file))                 
        for DER_id,DER_config in config_dict.items():
            self.check_DER_config(DER_config,DER_id)           
                               
    def modify_DER_parameters(self,parameter_ID):
        """Modify the DER parameters to parameters corresponding to the given parameter ID.
        
        Args:
             parameter_ID (str): User specified parameter ID (can be None). 
        """
        
        self.parameter_ID = self.get_DER_id(derId=parameter_ID)
        self.initialize_module_parameters()
        self.initialize_DER_model()
        
        self.initialize_states(ia0 = 0+0j, xa0 = 0+0j, ua0 = 0+0j,\
                               xDC0 = 0, xQ0 = 0, xPLL0 = 0.0,wte0 = 2*math.pi)
        
        self.initialize_derived_quantities()
        
        self.logger.info('{}:PV-DER parameters updated with parameters from  parameter dictionary {}!'.format(self.name,self.parameter_ID))
    
    def initialize_basic_specs(self):
        """Initialize number of ODEs and phases"""
        
        self.n_ODE = templates.DER_design_template[self.DER_model_type]['basic_specs']['n_ODE'] #23  #Number of ODE's
        self.n_phases = templates.DER_design_template[self.DER_model_type]['basic_specs']['n_phases'] #3 #Number of phases  
        
    def initialize_basic_options(self):
        """Initialize basic options"""  

        self.t_stable = self.DER_config['basic_options']['t_stable'] #0.5  #Time delay before activating logic for MPP, Volt-VAR control,  LVRT/LFRT 
        self.m_steady_state = self.DER_config['basic_options']['m_steady_state'] #0.96 #Expected duty cycle at steady state    
    
    def initialize_grid_measurements(self,DER_arguments):
        """Initialize inverter states.
        Args:
             gridVoltagePhaseA (complex): Value of gridVoltagePhaseA
             gridVoltagePhaseB (complex): Value of gridVoltagePhaseB
             gridVoltagePhaseC (complex): Value of gridVoltagePhaseC
             gridFrequency (float): Value of gridFrequency
        
        """        

        if not self.standAlone:
            if 'gridFrequency' in DER_arguments:
                self.gridFrequency = DER_arguments['gridFrequency']
            else:
                raise ValueError('Grid voltage source Frequency need to be supplied if model is not stand alone!')
            
            model_name = type(self).__name__            
            
            if templates.DER_design_template[model_name]['basic_specs']['n_phases'] >=1: #Check if model has one phase
                if 'gridVoltagePhaseA' in DER_arguments:
                    self.gridVoltagePhaseA = DER_arguments['gridVoltagePhaseA']/self.Vbase
                else:
                    raise ValueError('Grid voltage source phase A need to be supplied if model is not stand alone!')
               
            if templates.DER_design_template[model_name]['basic_specs']['n_phases'] >=2: #Check if model has 2 phases
                if templates.DER_design_template[model_name]['basic_specs']['unbalanced']: #Check if model is unbalanced
                    if 'gridVoltagePhaseB' in DER_arguments:
                        self.gridVoltagePhaseB = DER_arguments['gridVoltagePhaseB']/self.Vbase
                    else:
                        raise ValueError('Grid voltage source phase B need to be supplied if model is not stand alone!')
                else:
                    self.gridVoltagePhaseB = utility_functions.Ub_calc(self.gridVoltagePhaseA)
                    
            if templates.DER_design_template[model_name]['basic_specs']['n_phases'] >=3: #Check if model has 3 phases
                if templates.DER_design_template[model_name]['basic_specs']['unbalanced']: #Check if model is unbalanced
                    if 'gridVoltagePhaseC' in DER_arguments:
                        self.gridVoltagePhaseC = DER_arguments['gridVoltagePhaseC']/self.Vbase
                    else:
                        raise ValueError('Grid voltage source phase C need to be supplied if model is not stand alone!')
                else:
                    self.gridVoltagePhaseC = utility_functions.Uc_calc(self.gridVoltagePhaseA)
            
            if templates.DER_design_template[model_name]['basic_specs']['n_phases'] >=4: #Check if model has 3 phases
                raise ValueError('Model has more than 3 phases!')
                
    def initialize_states(self,DER_arguments):
        """Initialize inverter states.

        Args:
            ia0 (float): Initial current
            xa0 (float): Initial controller state
            ua0 (float): Initial controller state

        """        
        
        if 'ia' in DER_arguments:
            ia0 = DER_arguments['ia']
        else:
            ia0 = self.DER_config['initial_states']['iaR'] + 1j*self.DER_config['initial_states']['iaI']
        if 'xa' in DER_arguments:
            xa0 = DER_arguments['xa']
        else:
            xa0 = self.DER_config['initial_states']['xaR'] + 1j*self.DER_config['initial_states']['xaI']
        if 'ua' in DER_arguments:
            ua0 = DER_arguments['ua']
        else:
            ua0 = self.DER_config['initial_states']['uaR'] + 1j*self.DER_config['initial_states']['uaI']
        if 'xDC' in DER_arguments:
            xDC0 = DER_arguments['xDC']
        else:
            xDC0 = self.DER_config['initial_states']['xDC']
        if 'xQ' in DER_arguments:
            xQ0 = DER_arguments['xQ']
        else:
            xQ0 = self.DER_config['initial_states']['xQ']
        
        xPLL0 = self.DER_config['initial_states']['xPLL']
        wte0 = self.DER_config['initial_states']['wte']        
        
        self.Vdc = self.Vdc_ref  #DC link voltage        
        self.Ppv = self.Ppv_calc(self.Vdc_actual) #PV module power output    
        
         #Initialize all states with steady state values at current operating point
        if self.steady_state_initialization:
            #ia0,xa0,ua0,ib0,xb0,ub0,ic0,xc0,uc0,Vdc0,xDC0,xQ0,xPLL0,wte0
            self.steady_state_calc()
        else:
            #Phase a
            self.ia = ia0
            self.xa = xa0
            self.ua = ua0
            
            #DC link voltage and reactive power controller
            self.xDC = xDC0
            self.xQ = xQ0

            #PLL
            self.xPLL = xPLL0
            self.wte = wte0

            if type(self).__name__ == 'SolarPV_DER_ThreePhase':
                ib0 = utility_functions.Ub_calc(ia0)
                xb0 = utility_functions.Ub_calc(xa0)
                ub0 = utility_functions.Ub_calc(ua0)
                
                ic0 = utility_functions.Uc_calc(ia0)
                xc0 = utility_functions.Uc_calc(xa0)
                uc0 = utility_functions.Uc_calc(ua0)
        
                #Phase b
                self.ib = ib0
                self.xb = xb0  #Shift by -120 degrees
                self.ub = ub0

                #Phase c
                self.ic = ic0
                self.xc = xc0   #Shift by +120 degrees
                self.uc = uc0

    
    def initialize_derived_quantities(self):
        """Initialize quantities other than states."""
        
        self.update_voltages()
        self.update_power()        
        self.update_RMS()        
       
        self.update_iref() #Reference currents
        self.update_inverter_frequency(t=0.0)
    
    def initialize_DER_model(self,DER_arguments):
        """Initialize DER ratings.

        Args:
             DER_arguments (dict): Key word arguments.

        Raises:
             ValueError: If specified parameters correponding to `parameter_ID` is not available.
        """
        
        self.logger.debug('Creating inverter instance for DER with parameter ID:{}!'.format(self.parameter_ID))
        
        self.initialize_module_parameters()
        self.initialize_inverter_ratings(DER_arguments) #Initialize inverter ratings according to DER rating
        self.check_voltage() #Check if voltage is feasible
               
        
        self.initialize_circuit_parameters()
        self.check_circuit_parameters()  #Check PV-DER circuit parameters
            
        self.initialize_controller_gains()  #Initialize control loop gains according to DER rating           
    
    def initialize_inverter_ratings(self,DER_arguments):
        """Initialize inverter voltage and power ratings."""
        
        if self.parameter_ID in self.inverter_ratings:
            
            self.initialize_Sinverter(DER_arguments)
            self.initialize_Vdc(DER_arguments)
            self.initialize_Vac(DER_arguments)    
            self.initialize_Iac()
            
        else:
            raise ValueError('Inverter voltage, current, power ratings not available for parameter ID {}!'.format(self.parameter_ID))       
    
    def initialize_Sinverter(self,DER_arguments):
        """Initialize inverter power rating."""
        
        if 'powerRating' in DER_arguments:
            self.Sinverter_rated = DER_arguments['powerRating'] #Inverter rating in kVA
        elif 'Srated' in self.inverter_ratings[self.parameter_ID]:
            self.Sinverter_rated = self.inverter_ratings[self.parameter_ID]['Srated']
        else:
            raise ValueError('Sinverter rated not found for {}'.format(self.parameter_ID))
        
        self.Sinverter_nominal = (self.Sinverter_rated/BaseValues.Sbase) #Converting to p.u. value           
        
    def initialize_Vdc(self,DER_arguments):
        """Initialize DC side voltages."""
        
        if 'Vdcrated' in DER_arguments:
            self.Vdcrated = DER_arguments['Vdcrated'] #Rated DC link voltage in V
        else:
            self.Vdcrated = self.inverter_ratings[self.parameter_ID]['Vdcrated'] #Rated DC link voltage in V
        
        self.Vdcnominal = self.Vdcrated/self.Vdcbase   #Converting to p.u. value 
        self.Vdc_ref = self.set_Vdc_ref()
        self.Vdc_ref_new = self.set_Vdc_ref()        
    
    def initialize_Vac(self,DER_arguments):
        """Initialize AC side voltages."""
        
        if 'VrmsRating' in DER_arguments: #First check if rated voltage is given in arguments
            self.Vrms_rated = DER_arguments['VrmsRating'] #Inverter rating in kVA
        elif 'Vrmsrated' in self.inverter_ratings[self.parameter_ID]:
            self.Vrms_rated = self.inverter_ratings[self.parameter_ID]['Vrmsrated']
        else:
            raise ValueError('Vrms rated not found for {}'.format(self.parameter_ID))
        self.Varated = self.Vrms_rated*math.sqrt(2) #L-G peak to peak 
        
        self.Vanominal = self.Varated/BaseValues.Vbase #Converting to p.u. value
        self.Vrms_ref =  self.Vanominal/math.sqrt(2)   #Reference voltage          
    
    def initialize_Iac(self):
        """Initialize AC side currents."""
        
        self.Iarated = (self.Sinverter_rated/(self.n_phases*(self.Varated/math.sqrt(2))))*math.sqrt(2)
        #if type(self).__name__ == 'SolarPV_DER_SinglePhase':
        #    self.Iarated = (self.Sinverter_rated/(self.Varated/math.sqrt(2)))*math.sqrt(2)
        #elif type(self).__name__ == 'SolarPV_DER_ThreePhase' or type(self).__name__ == 'SolarPV_DER_ThreePhaseBalanced':
        #    self.Iarated = (self.Sinverter_rated/(3*(self.Varated/math.sqrt(2))))*math.sqrt(2)            
        
        self.Ioverload = self.inverter_ratings[self.parameter_ID]['Ioverload']  #Inverter current overload rating (Max 10s)            
        self.iref_limit = (self.Iarated/self.Ibase)*self.Ioverload #Maximum current reference
        
    def initialize_circuit_parameters(self):
        """Initialize C, Lf, and Rf parameters."""
        
        if self.parameter_ID in self.inverter_ratings:
        
            if self.standAlone:
                self.a = self.grid_model.Vgridrated/self.Varated  #Transformer turns ratio

            self.Rf_actual =  self.circuit_parameters[self.parameter_ID]['Rf_actual'] #Filter resistance
            self.Lf_actual = self.circuit_parameters[self.parameter_ID]['Lf_actual']  #Filter inductance
            self.C_actual =  self.circuit_parameters[self.parameter_ID]['C_actual']   #DC link capacitance
            self.Zf_actual =  self.Rf_actual + 1j*self.Lf_actual*BaseValues.wbase

            self.Rf = self.Rf_actual/BaseValues.Zbase        # VSC Filter resistance 
            self.Lf = self.Lf_actual/BaseValues.Lbase        # VSC Filter inductance
            self.Xf = (self.Lf_actual*BaseValues.wbase)/BaseValues.Zbase # VSC Filter reactance
            self.C = self.C_actual/BaseValues.Cbase          #DC link capacitor capacitance
            self.Zf = self.Zf_actual/BaseValues.Zbase

            #Interconnection to PCC - HV side
            #Actual values
            self.Z1_actual = self.circuit_parameters[self.parameter_ID]['Z1_actual']
            self.R1_actual = self.Z1_actual.real
            self.L1_actual = self.Z1_actual.imag/(2*math.pi*60.0)

            #Per-unit values
            self.R1 = self.R1_actual/BaseValues.Zbase  #Line/transformer resistance
            self.L1 = self.L1_actual/BaseValues.Lbase  #Line/transformer inductance
            self.Z1 =self.Z1_actual/BaseValues.Zbase    #Line/transformer impedance

            #Exhibit different behavior when used in standalone mode
            if self.standAlone:
                self.n_total_ODE = self.n_ODE + self.grid_model.n_ODE
                self.Zload1_actual =  self.events.load_events(t=0.0)          #Load at PCC   

            else:
                self.n_total_ODE = self.n_ODE
                self.Zload1_actual =  10e6+1j*0.0     #Using a large value of impedance to represent a no load condition

            self.Rload1_actual = self.Zload1_actual.real
            self.Lload1_actual = self.Z1_actual.imag/(2*math.pi*60.0)

            #Per-unit values
            self.Rload1 = self.Rload1_actual/BaseValues.Zbase  #resistance in per unit value
            self.Lload1 = self.Lload1_actual/BaseValues.Lbase  #inductance in per unit value
            self.Zload1 = self.Zload1_actual/BaseValues.Zbase

            self.transformer_name = 'transformer_'+str(self.ID)            
            
        else:
            raise ValueError('Inverter circuit_parameters not available for parameter ID: {}!'.format(self.parameter_ID))        
        
    def initialize_controller_gains(self):
        """Initialize controller settings."""
        
        if self.parameter_ID in self.controller_gains:
            
            #Current controller gains
            self.Kp_GCC = self.controller_gains[self.parameter_ID]['Kp_GCC'] #Current controller Proportional constant
            self.Ki_GCC = self.controller_gains[self.parameter_ID]['Ki_GCC']  #Current controller Integral constant
            self.wp =  self.controller_gains[self.parameter_ID]['wp']      #First order filter gain for GCC

            #Power (active and reactive) controller gains
            self.Kp_DC = self.controller_gains[self.parameter_ID]['Kp_DC']   #Active power controller Proportional constant
            self.Ki_DC = self.controller_gains[self.parameter_ID]['Ki_DC'] #Active power controller Integral constant
            self.Kp_Q = self.controller_gains[self.parameter_ID]['Kp_Q']  #Reactive power controller Proportional constant
            self.Ki_Q = self.controller_gains[self.parameter_ID]['Ki_Q']   #Reactive power controller Integral constant   
            
        else:
            raise ValueError('Controller gains not available for parameter ID: {}!'.format(self.parameter_ID))    
            
    def initialize_jacobian(self):
        """Create a Jacobian matrix with zero values."""
        
        self.J = np.zeros((self.n_total_ODE,self.n_total_ODE))
        self.varInd={}
        n=0
        
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
    
    def attach_grid_model(self,DER_arguments):
        """Attach a grid model to the PV-DER instance.

        Args:
             grid_model: An instance of `GridModel`.

        Returns:
             bool: Description of return value
        """
        
        #Connect grid instance only if working in stand alone mode
        if self.standAlone:
            if 'gridModel' in DER_arguments:
                self.grid_model = DER_arguments['gridModel']
                self.logger.debug('{}:Grid model {} attached since PV-DER instance is stand alone!'.format(self.name,DER_arguments['gridModel'].name))
            else:
                raise ValueError('`Grid` instance need to be provided in stand alone mode for creating `SolarPV_DER` instance`!')
        else: #Grid model is not connected
            self.logger.debug('{}:No grid model attached since PV-DER instance is not stand alone!'.format(self.name))
                             
    def check_voltage(self):
        """Method to check whether inverter voltage ratings are feasible.
        
        Raises:
             ValueError: If any of the specificied voltage ratings is infeasible.
        """
                             
        if not self.standAlone and abs(abs(self.gridVoltagePhaseA) - self.Vanominal)/self.Vanominal > 0.1:
            raise ValueError('The rated PCC-LV voltage {:.2f} V has more than 10% deviation from the voltage input from external program {:.2f} V!'.format(self.Vanominal, abs(self.gridVoltagePhaseA)))
                             
        if self.m_steady_state*(self.Vdcrated/2) < self.Varated:
            raise ValueError('The nominal DC link voltage {:.1f} V is not sufficient for the inverter to generate the nominal voltage at PCC - LV side {:.1f} V (L-G peak). Increase nominal DC link voltage to {:.1f} V.'.format(self.Vdcrated,self.Varated,math.ceil((self.Varated/self.m_steady_state)*2)))
        
        if self.m_steady_state*(self.Vdcmpp_min/2) < self.Varated:
            raise ValueError('The minimum DC link voltage {:.1f} V is not sufficient for the inverter to generate the nominal voltage at PCC - LV side {:.1f} V (L-G peak). Increase minimum DC link voltage to {:.1f} V.'.format(self.Vdcmpp_min,self.Varated,math.ceil((self.Varated/self.m_steady_state)*2)))
            
    def check_circuit_parameters(self):
        """Method to check whether inverter circuit parameter's are feasible."""
        
        _del_I1max = 0.1*self.Iarated
        _Lf_min = self.Vdcrated/(16*self.fswitching*_del_I1max)
        _del_I1max_actual = self.Vdcrated/(16*self.fswitching*self.Lf_actual)
        if _del_I1max_actual > _del_I1max:   #Check if ripple current is less than 10 %
            self.logger.debug('{}:Filter inductance {:.4f} H is acceptable since AC side current ripple is {:.2f}% (< 10%)'.format(self.name,_Lf_min,_del_I1max_actual/self.Iarated))
        else:
            self.logger.debug('{}:Warning:Filter inductance {:.4f} H results in AC side current ripple of {:.2f}% (> 10%)'.format(self.name,_Lf_min,_del_I1max_actual/self.Iarated))
        
        _I_ripple = (0.25*self.Vdcrated)/(self.Lf_actual*self.fswitching)  #Maximum ripple voltage (p-p) at DC link
        _V_ripple = self.Vdcrated/(32*self.Lf_actual*self.C_actual*(self.fswitching**2))  #Maximum ripple voltage (p-p) at DC link
        _V_ripple_percentage = (_V_ripple/self.Vdcrated)*100
        if _V_ripple_percentage <= 1.0:   #Check if voltage ripple on DC link is less than 1%
            self.logger.debug('{}:DC link capacitance of {:.4f} F is acceptable since voltage ripple is only {:.2f}% (< 1%)'.format(self.name,self.C_actual,_V_ripple_percentage))
        
        else:
            _V_ripple_ideal = self.Vdcrated*0.01  #1% ripple is acceptable
            _C = self.Vdcrated/(32*self.Lf_actual*_V_ripple_ideal*(self.fswitching**2))
            self.logger.debug('{}:Warning:DC link capacitance of {:.4f} F results in DC link voltage ripple of {:.2}% (> 1%)!Please use at least {} F.'.format(self.name,self.C_actual,_V_ripple_percentage,_C))
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

        #Sloss_filter = ((abs(ia)/math.sqrt(2))**2)*self.Zf 
    
        if type(self).__name__ == 'SolarPV_DER_ThreePhase':
            
            if not self.allow_unbalanced_m:
                mb = utility_functions.Ub_calc(ma)
                mc = utility_functions.Uc_calc(ma)
                
                ib = self.ib = utility_functions.Ub_calc(ia)
                ic = self.ic = utility_functions.Uc_calc(ia)
            
            elif self.allow_unbalanced_m:
                mbR = x[4]
                mbI = x[5]
                ibR = x[6]
                ibI = x[7]

                mcR = x[8]
                mcI = x[9]
                icR = x[10]
                icI = x[11]

                mb = mbR + 1j*mbI
                mc = mcR + 1j*mcI                

                ib = self.ib = ibR + 1j*ibI
                ic = self.ic = icR + 1j*icI
            
            vtb = mb*(self.Vdc/2)
            vtc = mc*(self.Vdc/2)
            
            vb = self.vb_calc()
            vc = self.vc_calc()
            
            St = St + (vtb*ib.conjugate() + vtc*ic.conjugate())/2
            S_PCC = S_PCC + (vb*ib.conjugate() + vc*ic.conjugate())/2
            
            #Sloss_filter = Sloss_filter + ((abs(ib)/math.sqrt(2))**2 + (abs(ic)/math.sqrt(2))**2)*self.Zf         
        
        P_error = (St.real - self.Ppv)**2
        Q_error = (S_PCC.imag - self.Q_ref)**2
        
        #Qloss_filter_expected = self.n_phases*((self.Ppv/(self.n_phases*self.Vrms_ref))**2)*self.Xf
        #Ploss_filter_expected = self.n_phases*((self.Ppv/(self.n_phases*self.Vrms_ref))**2)*self.Rf
                
        #P_PCC_error = ((S_PCC.real + Sloss_filter.real)   - self.Ppv)**2         
        #S_error = (abs(St -(S_PCC+Sloss_filter)))**2
        #Q_error_filter_expected =  (St.imag - Qloss_filter_expected)**2 
        
        
        if type(self).__name__ == 'SolarPV_DER_ThreePhase':
            del_1 = utility_functions.relative_phase_calc(ma,mb)
            del_2 = utility_functions.relative_phase_calc(ma,mc)
            del_3 = utility_functions.relative_phase_calc(mb,mc)
            
            if del_1 > math.pi:
                del_1 = abs(del_1 - 2*math.pi)

            if del_2 > math.pi:
                del_2 = abs(del_2 - 2*math.pi)

            if del_3> math.pi:
                del_3 = abs(del_3 - 2*math.pi)
            
            if self.allow_unbalanced_m:
                m_error = (abs((va+vb+vc)- (vta+vtb+vtc)))**2 
            else:
                m_error = (del_1 - PHASE_DIFFERENCE_120)**2 + (del_2 - PHASE_DIFFERENCE_120)**2 + (del_3 - PHASE_DIFFERENCE_120)**2 +\
                          (abs(ma) - abs(mb))**2 + (abs(ma) - abs(mc))**2 + (abs(mb) - abs(mc))**2 
            
            m_error = m_error
                              
            i_error = (abs(ia + ib+ic) +  abs(vta-va - ia*self.Zf)+  abs(vtb-vb - ib*self.Zf)+abs(vtc-vc - ic*self.Zf))**2
        else:
            m_error = 0.0
            i_error = (abs(vta-va - ia*self.Zf))**2
        
        return  P_error + Q_error + m_error + i_error #+S_error + P_PCC_error

    def steady_state_calc(self):
        """Find duty cycle and inverter current that minimize steady state error and return steady state values."""        
        
        self.logger.debug('Solving for steady state at current operating point.') 
        
        if self.standAlone:
            x0 = [self.steadystate_values[self.parameter_ID]['maR0'],self.steadystate_values[self.parameter_ID]['maI0'],
                  self.steadystate_values[self.parameter_ID]['iaR0'],self.steadystate_values[self.parameter_ID]['iaI0']]        
        
        else:
            va = self.va_calc()
            x0 = [va.real,va.imag,
                  self.steadystate_values[self.parameter_ID]['iaR0'],self.steadystate_values[self.parameter_ID]['iaI0']]    
        
        if self.allow_unbalanced_m:
            self.logger.info('Using unbalanced option - steady state duty cycles may be unbalanced.') 
            if self.standAlone:
                x0.extend([self.steadystate_values[self.parameter_ID]['mbR0'],self.steadystate_values[self.parameter_ID]['mbI0'],
                           self.steadystate_values[self.parameter_ID]['ibR0'],self.steadystate_values[self.parameter_ID]['ibI0'],
                           self.steadystate_values[self.parameter_ID]['mcR0'],self.steadystate_values[self.parameter_ID]['mcI0'],
                           self.steadystate_values[self.parameter_ID]['icR0'],self.steadystate_values[self.parameter_ID]['icI0']])       
            else:
                vb = self.vb_calc()
                vc = self.vc_calc()
            
                x0.extend([vb.real,vb.imag,
                           self.steadystate_values[self.parameter_ID]['ibR0'],self.steadystate_values[self.parameter_ID]['ibI0'],
                           vc.real,vc.imag,
                           self.steadystate_values[self.parameter_ID]['icR0'],self.steadystate_values[self.parameter_ID]['icI0']])      
            
        x0 = np.array(x0)
        
        disp = bool(self.verbosity == 'DEBUG')
        self.solver_spec[self.steadystate_solver].update({'disp':bool(self.verbosity == 'DEBUG')})                   
        result = minimize(self.power_error_calc, x0,
                          method=self.steadystate_solver,options=self.solver_spec[self.steadystate_solver])
        
        if not result.success:
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
            
            if not self.allow_unbalanced_m:
                mb0 = utility_functions.Ub_calc(ma0)
                mc0 = utility_functions.Uc_calc(ma0)
                
                self.ib = utility_functions.Ub_calc(self.ia)
                self.ic = utility_functions.Uc_calc(self.ia)            
            
            else:
                mb0 = result.x[4] + 1j*result.x[5]
                mc0 = result.x[8] + 1j*result.x[9]
                self.ib = result.x[6] + 1j*result.x[7]
                self.ic = result.x[10] + 1j*result.x[11]                
            
            self.xb = mb0
            self.xc = mc0
            
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
        
        self.logger.debug('{}:Steady state values for operating point defined by Ppv:{:.2f} W, Vdc:{:.2f} V, va:{:.2f} V found at:'.format(self.name,self.Ppv*self.Sbase,self.Vdc*self.Vdcbase,self.va*self.Vbase))
            
        if self.verbosity == 'DEBUG':
            self.show_PV_DER_states(quantity='power')
            self.show_PV_DER_states(quantity='duty cycle')
            self.show_PV_DER_states(quantity='voltage')
            self.show_PV_DER_states(quantity='current')
            
    def check_jacobian(self,t=0.0):
        """Compare analytical and numerical Jacobian of the ODE model."""
        
        #Calculate numerical Jacobian using finite differences
        x = self.y0
        x0 = x.copy()

        eps = 1e-6
        Jn = np.zeros([len(x), len(x)], dtype = np.float)

        for i in range(len(x)):
            x1 = x.copy()
            x2 = x.copy()

            x1[i] += eps
            x2[i] -= eps

            f1 = self.ODE_model(x1,t)
            f2 = self.ODE_model(x2,t)

            Jn[ : , i] = (f1 - f2) / (2 * eps)
        
        Ja = self.jac_ODE_model(x0,t) #Calculate analytical Jacobian
        
        #Compare numerical Jacobian with analytical Jacobian
        err=np.abs(Ja-Jn)
        err[err<1e-2]=0
        ind=np.where(err>0)

        nErr=len(ind[0])
        err=[]
        
        if nErr > 0:
            print('Differences in analytical and numerical Jacobian found in {} entries!'.format(nErr))
            jac_CHECK = False
            for r,c in zip(ind[0],ind[1]):
                err.append([r,c,Jn[r,c],Ja[r,c]])
                #print(r,c,Jn[r,c],Ja[r,c])
                r_char = list(self.varInd.keys())[list(self.varInd.values()).index(r)]
                c_char = list(self.varInd.keys())[list(self.varInd.values()).index(c)]
                print('J[{}][{}]--Jn:{},Ja:{}'.format(r_char,c_char,Jn[r,c],Ja[r,c]))
        else:
            jac_CHECK = True
            print('No differences in analytical and numerical Jacobian!')
        
        return jac_CHECK,Jn,Ja
