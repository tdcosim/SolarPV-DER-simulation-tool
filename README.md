**Status:** Maintenance (expect bug fixes and minor updates)
# Utilitiy for simulating dynamics of PV-DER

Solar photovoltaic distributed energy resources (PV-DER) are power electronic inverter based generation (IBG) connected to the electric power distribution system (eg. roof top solar PV systems). This utility can be used to simulate the behaviour a single DER connected to a stiff voltage source as shown in the following schematic:

![schematic of PV-DER](PVDER_schematic.png)

## Basics
The dynamics of the DER are modelled using dynamic phasors.

## Installation
You can install the module using following commands:
```
git clone https://github.com/sibyjackgrove/SolarPV-DER-simulation-utility.git
cd SolarPV-DER-simulation-utility
pip install -e .
```
## Using the module
The module can be imported as shown below:
```
import pvder
```
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
