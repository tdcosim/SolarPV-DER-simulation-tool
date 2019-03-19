**Status:** Maintenance (expect bug fixes and minor updates)
# Gym environment for PV-DER

Solar photovoltaic distributed energy resources (PV-DER) are power electronic inverter based generation (IBG) connected to the electric power distribution system (eg. roof top solar PV systems). This environment consists of a single DER connected to a stiff voltage source as shown in the following schematic:

![schematic of PV-DER](PVDER_schematic.png)

## Basics
The dynamics of the DER are modelled using dynamic phasors. One step is equivalent to one cycle (0.016 s).

## Installation
You can install the environment using following commands:
```
git clone https://github.com/sibyjackgrove/gym-SolarPVDER-environment.git
cd gym-SolarPVDER-environment
pip install -e .
```
## Using the environment
The environment can be instantiated just like any other OpenAI Gym environment as show below:
```
import gym
import gym_PVDER
env = gym.make('PVDER-v0')
```
## Citation
If you use this code please cite it as:
```
@misc{gym-PVDER,
  title = {{gym-SolarPVDER-environment}: A environment for solar photovoltaic distributed energy resources},
  author = "{Siby Jose Plathottam}",
  howpublished = {\url{https://github.com/sibyjackgrove/gym-SolarPVDER-environment}},
  url = "https://github.com/tensorflow/agents",
  year = 2019,
  note = "[Online; accessed 18-March-2019]"
}
```
