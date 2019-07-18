**Status:** Expect regular updates and bug fixes.
# Utility for simulating dynamics of PV-DER

[![Build Status](https://travis-ci.org/sibyjackgrove/SolarPV-DER-simulation-utility.svg?branch=master)](https://travis-ci.org/sibyjackgrove/SolarPV-DER-simulation-utility)
[![CodeFactor](https://www.codefactor.io/repository/github/sibyjackgrove/solarpv-der-simulation-utility/badge)](https://www.codefactor.io/repository/github/sibyjackgrove/solarpv-der-simulation-utility)

Solar photovoltaic distributed energy resources (PV-DER) are power electronic inverter based generation (IBG) connected to the electric power distribution system (eg. roof top solar PV systems). This utility can be used to simulate the dynamics of a single DER connected to a stiff voltage source as shown in the following schematic:

![schematic of PV-DER](PVDER_schematic.png)

## Basics
The dynamics of the DER are modelled using dynamic phasors. Detailed description of the concepts behind this utility can be found in the IEEE publication [Dynamic Modeling of Solar PV Systems for Distribution System Stability Analysis](https://www.researchgate.net/publication/333985171_Dynamic_Modeling_of_Solar_PV_Systems_for_Distribution_System_Stability_Analysis) and detailed list of equations can be found in the [Model specification document.](docs/PV_DER_model_specification_rev3.docx)

## Links
* Source code repository: https://github.com/sibyjackgrove/SolarPV-DER-simulation-utility
* API Documentation: https://solarpv-der-simulation-utility.readthedocs.io/en/latest/

## Installation
You can install the module directly from GitHub with following commands:
```
git clone https://github.com/sibyjackgrove/SolarPV-DER-simulation-utility.git
cd SolarPV-DER-simulation-utility
pip install -e .
```

## Use cases
Following projects are using Solar PV-DER simulation utility:
1. [Argonne Transmission and Distribution systems Co-Simulation tool (TDcoSim)](https://github.com/tdcosim/TDcoSim)
2. [OpenAI Gym Distributed Energy Resource Environment  (Gym-DER)](https://github.com/sibyjackgrove/gym-SolarPVDER-environment)

## Using the module
The module can be imported as a normal python module:

```python
import pvder
```
The following features are available currently:
1. Single phase or three phase DER models (phase voltages may be unbalanced).
2. Run simulation in stand alone mode with internal grid voltage source model.
3. Run simulation in loop mode where grid voltage and frequency is supplied every time step by outside program.
4. Visualize simulation results for voltages, current, active, and reactive power.
5. Introduce solar insolation events (in all modes), grid voltage, and frequency change events (in stand alone mode).
6. Enable Low voltage ride through (LVRT) and Volt-VAR control logic.

### Using the stand alone single phase DER model with 10 kW power rating
The following steps are required:
1. First import the following classes:
```
from pvder.DER_components_single_phase import SolarPV_DER_SinglePhase
from pvder.grid_components import Grid
from pvder.dynamic_simulation import DynamicSimulation
from pvder.simulation_events import SimulationEvents
from pvder.simulation_utilities import SimulationResults
```
1. Create a **_SimulationEvents_** object: This object is used to add or remove disturbance events occurs during the simulation.
```
events = SimulationEvents()
```
2. Create a **Grid** object: This object describes the steady state model for the grid voltage source. It needs to be supplied with an **_SimulationEvents_** object.
```
grid = Grid(events=events)
```
3. Create a **SolarPV_DER_SinglePhase** or **SolarPV_DER_ThreePhase** object: This object describes the dynamic DER model. It needs both an **_SimulationEvents_** object and a **Grid** object. The power rating of the DER are should also be provided.
```
PV_DER = SolarPV_DER_SinglePhase(grid_model=grid,events=events,Sinverter_rated = 10.0e3,standAlone = True)
```
4. Create a **DynamicSimulation** object: This object runs the simulation and stores the solution. It takes **_SimulationEvents_**, **Grid** and, **SolarPV_DER_SinglePhase** objects as arguments.
```
sim = DynamicSimulation(grid_model=grid,PV_model=PV_DER,events = events)
```
5. Create a **SimulationResults** object: This object is used to visualize the simulation results.
```
results = SimulationResults(simulation = sim)
```
6. Add an event (for e.g. solar insolation change at 10.0 s):
```
events.add_solar_event(10,90)
```
7. Specify simulation flags (for e.g. set the DEBUG_SIMULATION and DEBUG_POWER flag to true to observe the power at each time step.):
```
sim.DEBUG_SIMULATION = False
sim.DEBUG_POWER = False
```
8. Specify simulation stop time (for e.g. 20.0 s):
```
sim.tStop = 20.0
```
9. Run the simulation:
```
sim.run_simulation()
```
10. Visualize the results (for e.g. the power output at PCC-LV side):
```
results.PER_UNIT = False
results.plot_DER_simulation(plot_type='active_power_Ppv_Pac_PCC')
```
Try it out in Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/sibyjackgrove/SolarPV-DER-simulation-utility/blob/master/examples/PV-DER_usage_example.ipynb)

## Module details
A schematic of the relationship between differen classes in the module is shown in the figure below:
![schematic of software architecture](docs/software_architecture.png)

Dependencies: SciPy, Numpy, Matlplotlib

## Issues
Please feel free to raise an issue when bugs are encountered or if you are need further documentation.

## Who is responsible?

**Core developer:**
- Siby Jose Plathottam splathottam@anl.gov

**Support:**

- Karthikeyan Balasubramaniam kbalasubramaniam@anl.gov

## Acknowledgement
The authors would like to acknowledge [Shrirang Abhyankar](https://github.com/abhyshr) and Puspal Hazra for their contribution.

## Citation
If you use this code please cite it as:
```
@misc{pvder,
  title = {{SolarPV-DER-simulation-utility}: A simulation utility for or solar photovoltaic distributed energy resources},
  author = "{Siby Jose Plathottam,Karthikeyan Balasubramaniam}",
  howpublished = {\url{https://github.com/sibyjackgrove/SolarPV-DER-simulation-utility}},
  url = "https://github.com/sibyjackgrove/SolarPV-DER-simulation-utility",
  year = 2019,
  note = "[Online; accessed 19-March-2019]"
}
```
### Copyright and License
Copyright © 2019, UChicago Argonne, LLC

Photovoltaic Distributed Energy Resource (PV-DER) Simulation Utility is distributed under the terms of [BSD-3 OSS License.](LICENSE)
