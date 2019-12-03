# PV-DER Simulation Utility initialization, variables, flags, and methods

This document is a supplement to the [Jupyter notebook usage example](../examples/PV-DER_usage_example.ipynb). It provides information on important arguments during object creation, functionalities of various binary flags/variables (object attributes), and attached methods.

### Common initialization arguments

1. **verbosity (string):** Logging level for objects (available: 'DEBUG','INFO,'WARNING').
2. **identifier (string):** User defined identification string that will be added to the object name. 
### DER model objects
Object type name: *SolarPV_DER_ThreePhase*,*SolarPV_DER_ThreePhaseBalanced*,*SolarPV_DER_SinglePhase*
#### Essential initialization arguments
1. **Sinverter_rated (float):** Inverter power rating in W (Available:  10.0e3, 50.0e3, 250.0e3, default: 50.0e3).
1. **Vrms_rated (float):** Inverter voltage rating (L-G RMS) in V (default: None).
1. **standAlone (Boolean):** If this flag is **True**, the DER model is expected to work as a stand alone model connected to a stif voltage source. If **False**, the DER model expects to recieve voltage  values from an external program (default: **True**).
1. **STEADY_STATE_INITIALIZATION (Boolean):** If this flag is **True**, the states in the  dynamic DER model will be initialized with values corresponding to its steady state operating point. If **False** the states will be initialized with zeroes (default: **False**).

#### Essential variables and flags
1. **LVRT_ENABLE (Boolean):** If this flag is **True**, the low voltage ride through and protection logic will be enabled. If **False** the DER instance will neither trip nor enter momentary cessation when abnormal low voltage conditions are encountered  (default: **True**).
2. **HVRT_ENABLE (Boolean):** If this flag is **True**, the high voltage ride through and protection logic will be enabled. If **False** the DER instance will neither trip nor enter momentary cessation when abnormal high voltage conditions are encountered (default: **True**).
2. **LFRT_ENABLE (Boolean):** If this flag is **True**, the low frequency ride through and protection logic will be enabled. If **False** the DER instance will never trip when abnormal low frequency conditions are encountered (default: **False**).
3. **VOLT_VAR_ENABLE (Boolean):** If this flag is **True**, Volt-VAR control is enabled during voltage sags within a specified voltage range. If **False** the DER will neither supply or absorb reactive power when voltage sags are encountered (default: **False**).
4. **use_frequency_estimate (Boolean):** If this flag is **True**, grid frequency is estimated using the difference between phase angles at two consecutive time steps (default: **False**).

#### Essential methods
1. **show_PV_DER_states(quantity):** Show the values for the specified DER variable (default: 'voltage'). 
1. **show_PV_DER_parameters(parameter_type):** Show the values for the specified DER parameter (default: 'inverter_ratings'). 
2. **initialize_parameter_dict(parameter_ID,source_parameter_ID):** Initialize a new parameter dictionary. 
2. **update_parameter_dict(parameter_ID,parameter_type,parameter_dict):** Update an existing parameter dictionary with new values.

### Dynamic simulation model objects
Object types: *DynamicSimulation*

#### Essential initialization arguments
1. **LOOP_MODE (Boolean):** If this flag is **True**,  the simulation is expected to run in a loop and get updated grid voltages at every iteration in the loop (default: **False**).
2. **COLLECT_SOLUTION (Boolean):** If this flag is **True**,  the time series data for states and other variables are collected at the end of every call to the ODE solver (default: **True**).

#### Essential variables and flags
1. **tStart (float):** Start time for simulation in seconds (default: 0.0 s).
2. **tEnd (float):** End time for simulation in seconds (default: 0.5 s).
2. **tInc (float):** Time step for simulation in seconds (default: 0.001 s).
1. **jacFlag (Boolean):**  If this flag is **True**,  the analytical Jacobian will be passed to the SciPy ODE solver which may improve solution time. If this flag is **False** the solver will have to numerically calculate the Jacobian (default: False).
2. **DEBUG_SIMULATION (Boolean):** If this flag is **True**, the value of the model variables at each time step will be printed to the terminal at each time step. If this flag is **False** only information from ride through logic will be printed (default: False).
2. **DEBUG_SOLVER (Boolean):** If this flag is **True**, solution status from ODE solver is printed during each call to solver. If  **False**, solution status will only be printed if there is an exception (default: False).
2. **PER_UNIT (Boolean):** If this flag is **True**, all the displayed electrical quantities will be in per unit values. If **False**, all the displayed quantities will be in actual values (default: True).
#### Essential methods
1a. **run_simulation():** If LOOP_MODE is True the simulation is run from **tStart** to **tEnd** with time step of **tInc**. 
1b. **run_simulation(gridVoltagePhaseA, gridVoltagePhaseB, gridVoltagePhaseC, y0, t):** If LOOP_MODE is True the voltages, states, and time steps need to to be provided at every iteration.


### Simulation events object
Object type name: *SimulationEvents*
**Note:** Grid events are introduced during the simulation only if **standAlone** flag is True in *DER model* and **LOOP_MODE** is False in *DynamicSimulation*.

#### Essential initialization arguments
1. **SOLAR_EVENT_ENABLE (Boolean):** If this flag is **True**, any solar events that were added will be introduced during the simulation (default: **True**).
1. **GRID_EVENT_ENABLE  (Boolean):** If this flag is **True**, any grid events that were added will be introduced during the simulation (default: **True**).

#### Essential methods

1. **add_grid_event(T, Vgrid, Vgrid_angle, fgrid):** Add grid event at T (default: Vgrid=1.0, Vgrid_angle = 0.0, fgrid = 60). 
2. **add_solar_event(T, Sinsol):** Add solar event at T (default: Sinsol=100). 
3. **remove_grid_event(T):** Remove all grid events (voltage, phase angle, or frequency) at T.
4. **remove_solar_event(T):** Remove solar event at T.
5. **show_events():** Show list of all simulation events in chronological order.

### Results object
Object type name: *SimulationResults*

#### Essential initialization arguments
1. **PER_UNIT (Boolean):** If this flag is **True**, all the plots for electrical quantities will be in per unit values. If **False**, all the displayed quantities will be in actual values (default: **True**).

#### Essential variables
1. **font_size (int):** Font size of text within the plots (default: 18).

#### Essential methods
1. **plot_DER_simulation(plot_type):** Time series plot for specified electrical quantity (Available: 'power','active_power','active_power_Ppv_Pac_PCC','reactive_power','reactive_power_Q_PCC','voltage','voltage_Vdc','voltage_HV','voltage_HV_imbalance','voltage_LV','voltage_Vpcclv','current','duty_cycle','voltage_Vpcclv_all_phases'; default:'power', )
