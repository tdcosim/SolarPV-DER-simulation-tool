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
DEFAULT_Vreconnect_LV =  0.95  #Lower limit of voltage required for reconnection after momentary cessation
DEFAULT_Vreconnect_HV =  1.02  #Uppder limit of voltage required for reconnection after momentary cessation

#Frequency estimation
DEFAULT_USE_FREQUENCY_ESTIMATE = False

#Solver options
DEFAULT_DELTA_T = 0.001 #Simulation time step