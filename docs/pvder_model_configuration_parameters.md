# DER model parameters

The parameters used in the DER model and available for customization through the config file are are described here. These parameters can impact both the dynamic and steady state behaviour of a DER model. If the user doesn't specify any parameter for a given model type, they are automatically populated from the default values.

* **Module parameters**
  * *Np/Ns (int):* Parallel/series connected solar cells. Determines rated power output from the model.
  * *Vdcmpp0 (float):* Default maximum power point DC link voltage.

* * **Inverter ratings**
  * *Srated (float):*  Rated apparent power rating of the inverter (unit: kVA).
  * *Vrmsrated (float):* Rated voltage (RMS) of the inverter (unit: Volts).
  * *Vdcrated (float):* Rated DC link voltage of the inverter (unit: Volts).
  * *Ioverload (float):* Overload capacity of the inverter.

* * **Circuit parameters**
  * *Rf_actual (float):*  Inverter output filter resistance (unit: Ohm).
  * *Lf_actual (float):*  Inverter output filter inductance (unit: Henry).
  * *C_actual (float):*  Inverter DC link capacitor capacitance (unit: Farad).
  
* * **Controller gains**
  * *Kp_GCC (float):* Current controller propotional gain.
  * *Ki_GCC (float):* Current controller integral gain.
  * *wp (float):* Current controller input filter time constant.
  * *Kp_P (float):* Active power controller propotional gain (only available in constant Vdc models).
  * *Ki_P (float):* Active power controller integral gain (only available in constant Vdc models).
  * *Kp_DC (float):* DC link voltage controller propotional gain .
  * *Ki_DC (float):* DC link voltage controller integral gain .
  * *Kp_Q (float):* Reactive power controller propotional gain. 
  * *Ki_Q (float):* Reactive power controller integral gain.
