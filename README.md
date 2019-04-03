
**Status:** Maintenance (expect bug fixes and minor updates)
# Utilitiy for simulating dynamics of PV-DER

Solar photovoltaic distributed energy resources (PV-DER) are power electronic inverter based generation (IBG) connected to the electric power distribution system (eg. roof top solar PV systems). This utility can be used to simulate the behaviour a single DER connected to a stiff voltage source as shown in the following schematic:

![schematic of PV-DER](PVDER_schematic.png)

## Basics
The dynamics of the DER are modelled using dynamic phasors. Detailed description of the concepts behind this utility can be found in the IEEE publication **Dynamic Modeling of Solar PV Systems for Distribution System Stability Analysis** and detailed list of equations can be found in the [Model specification document.](docs/PV_DER_model_specification_rev3.docx)

## Installation
You can install the module directly from github with following commands:
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
Try it out in Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/sibyjackgrove/SolarPV-DER-simulation-utility/blob/master/examples/PV-DER_usage_example.ipynb)

## Module details
A schematic of the relationship between differen classes in the module is shown in the figure below:
![schematic of software architecture](examples/software_architecture.png)

Dependencies: SciPy, Numpy, Matlplotlib

## Who is responsible?
- Siby Jose Plathottam splathottam@anl.gov
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
