"""Store configuration options."""

#Results options

# Changing these adjusts the size and layout of the plot
FIGURE_WIDTH = 8
FIGURE_HEIGHT = 8
FIGURE_DPI = 1200


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

#Steady state solver options
DEFAULT_STEADYSTATE_SOLVER = 'SLSQP'