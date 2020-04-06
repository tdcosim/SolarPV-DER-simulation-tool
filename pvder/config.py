"""Store configuration options."""

#Results options

# Changing these adjusts the size and layout of the plot
FIGURE_WIDTH = 8
FIGURE_HEIGHT = 8
FIGURE_DPI = 1200

#Logging
DEFAULT_LOGGING_LEVEL = 'INFO'

#Model options

#Changing these impacts the how the parameter dictionaries for DER models are initialized
DEFAULT_PARAMETER_ID_1PH = '10'
DEFAULT_PARAMETER_ID_3PH = '50'

DEFAULT_Ioverload = 1.3

#Smart inverter feature options
DEFAULT_Vthreshold_high_limit = 1.5 #percentage
DEFAULT_tthreshold_high_limit = 100.0 #seconds

DEFAULT_tdisconnect_low_limit = 1/120.0 #seconds #Minimum time required to initiate DER disconnection (output cessation)
DEFAULT_treconnect_low_limit = 0.4 #seconds #Minimum time required to initiate DER reconnection (output restoration)

DEFAULT_freconnect_LF = 58.8
DEFAULT_freconnect_HF = 61.2

#DC link voltage control
DEFAULT_del_Vdc_ref = 2.0
DEFAULT_del_t_Vdc_ref = 0.5


#Frequency estimation
DEFAULT_USE_FREQUENCY_ESTIMATE = False

#Solver options
DEFAULT_DELTA_T = 0.001 #Simulation time step
DEFAULT_max_steps = 1000 #Max steps to be used by solver before producing error

#Steady state solver options
DEFAULT_STEADYSTATE_SOLVER = 'SLSQP'


#Default DER parameters
DEFAULT_Z1_actual= 0.0019 + 1j*0.0561
DEFAULT_basic_specs = {"Sinsol0":1.0}
DEFAULT_module_parameters = {'Np':2,'Ns':1000,'Vdcmpp0':750.0,'Vdcmpp_min': 650.0,'Vdcmpp_max': 800.0}
DEFAULT_inverter_ratings = {"Vdcrated":550.0,"Ioverload":DEFAULT_Ioverload,"Vrmsrated":177.0}
DEFAULT_circuit_parameters = {"Rf_actual":0.002,"Lf_actual":25.0e-6,
                              "C_actual":300.0e-6,"R1_actual":DEFAULT_Z1_actual.real,"X1_actual":DEFAULT_Z1_actual.imag}
DEFAULT_controller_gains ={"Kp_GCC":12000.0,"Ki_GCC":4000.0,
                               "Kp_DC":-4.0,"Ki_DC":-20.0,
                               "Kp_Q":0.4,"Ki_Q":20.0,"wp": 20e4}
DEFAULT_steadystate_values = {"maR0":0.7,"maI0":0.0,"iaR0":0.5,"iaI0":0.01}
DEFAULT_initial_values = {"iaR0":0,"iaI0":0.0,"xaR0":0.0,"xaI0":0.0,"uaR0":0.0,"uaI0":0.0,
                          "xDC0":0.0,"xQ0":0.0,"xPLL0":0.0,"wte0":6.28}

#Default VRT config
default_RT_config =  {'LVRT':{'0':{'V_threshold':0.5,
                                    't_threshold':1.0,
                                    'mode':'mandatory_operation',
                                    't_start':0.0,
                                    'threshold_breach':False},
                                '1':{'V_threshold':0.7,
                                     't_threshold':3.5,
                                     'mode':'momentary_cessation', #'momentary_cessation'
                                     't_start':0.0,
                                     'threshold_breach':False},
                                '2':{'V_threshold':0.88,
                                     't_threshold':5.0,
                                     'mode':'mandatory_operation',
                                     't_start':0.0,
                                     'threshold_breach':False},
                              },                       
                       'HVRT':{'0':{'V_threshold':1.12,
                                    't_threshold':0.5,
                                    'mode':'mandatory_operation',
                                    't_start':0.0,
                                    'threshold_breach':False},
                               '1':{'V_threshold':1.06,
                                     't_threshold':1.0,
                                     'mode':'mandatory_operation',
                                     't_start':0.0,
                                     'threshold_breach':False},
                              },
                       'OUTPUT_CESSATION_DELAY':0.01,
                       'OUTPUT_RESTORE_DELAY':1.75,
                       'RESTORE_Vdc':False}