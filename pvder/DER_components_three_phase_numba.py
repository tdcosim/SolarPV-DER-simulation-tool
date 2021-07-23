"""Three phase PV-DER components."""

from __future__ import division
import numpy as np
import math
import cmath
import scipy
import six
import pdb
import warnings
import numba

from pvder.DER_components import SolarPVDER,PVModule
from pvder.DER_check_and_initialize import PVDER_SetupUtilities
from pvder.DER_features import PVDER_SmartFeatures
from pvder.DER_utilities import PVDER_ModelUtilities
from pvder.grid_components import BaseValues
from pvder import utility_functions_numba,utility_functions
from pvder import defaults,templates,constants
from pvder.logutil import LogUtil

debug_flag =  False
cache_flag =  True

iaR_id = 0
iaI_id = 1
xaR_id = 2
xaI_id = 3
uaR_id = 4
uaI_id = 5
ibR_id = 6
ibI_id = 7
xbR_id = 8
xbI_id = 9
ubR_id = 10
ubI_id = 11
icR_id = 12
icI_id = 13
xcR_id = 14
xcI_id = 15
ucR_id = 16
ucI_id = 17
Vdc_id = 18
xDC_id = 19
xQ_id = 20
xPLL_id = 21
wte_id  = 22

@numba.njit(cache=cache_flag)
def update_inverter_states(iaR,iaI,xaR,xaI,uaR,uaI,
						   ibR,ibI,xbR,xbI,ubR,ubI,
						   icR,icI,xcR,xcI,ucR,ucI,
						   Vdc,xDC,xQ,xPLL,wte):
	"""Update inverter states"""

	ia = iaR + 1j*iaI
	xa = xaR + 1j*xaI
	ua = uaR + 1j*uaI
	
	ib = ibR + 1j*ibI
	xb = xbR + 1j*xbI
	ub = ubR + 1j*ubI
	
	ic = icR + 1j*icI
	xc = xcR + 1j*xcI
	uc = ucR + 1j*ucI
	
	Vdc = Vdc
	xDC = xDC
	xQ = xQ
	
	xPLL = xPLL
	wte = wte
	
	return ia,xa,ua,ib,xb,ub,ic,xc,uc,Vdc,xDC,xQ,xPLL,wte

@numba.njit(cache=cache_flag)
def Iph_calc(Tactual,Sinsol):
	"""Photocurrent from a single cell for given insolation and temperature."""
	return (constants.Iscr+(constants.Kv*(Tactual-constants.T0)))*(Sinsol/100.0)

@numba.njit(cache=cache_flag)
def Ppv_calc(Tactual,Iph,Vdc_actual,Np,Ns,Sbase):
	"""PV panel power output from  solar insolation.
		
		Args:
		  Vdc_actual (float): DC link voltage in volts
		  
		Returns:
			 float: Power output from PV module in p.u.
	"""
	Ipv = (Np*Iph)-(Np*constants.Irs*(math.exp((constants.q*Vdc_actual)/(constants.k*Tactual*constants.A*Ns))-1)) #Faster  with Pure Python functions
	
	return max(0,(Ipv*Vdc_actual))/Sbase
	#return utility_functions.Ppv_calc(self.Iph,self.Np,self.Ns,Vdc_actual,self.Tactual,Grid.Sbase)

@numba.njit(cache=cache_flag)
def diR_calc(i,vt,v,winv,Rf,Lf,wbase):
	"""Calcluate derivative of inverter branch current - real component - iR"""
	diR = (1/Lf)*(-Rf*i.real - v.real + vt.real) + (winv/wbase)*i.imag 
	
	return diR

@numba.njit(cache=cache_flag)
def diI_calc(i,vt,v,winv,Rf,Lf,wbase):
	"""Calcluate derivative of iaR"""
	diI = (1/Lf)*(-Rf*i.imag - v.imag + vt.imag) - (winv/wbase)*i.real
	
	return diI

@numba.njit(cache=cache_flag)
def duR_calc(i,u,i_ref,wp):
	"""Calculate derivate of 1st order filter- real component"""	
	duR = (wp)*(-u.real +	i_ref.real - i.real)
	
	return duR

@numba.njit(cache=cache_flag)
def duI_calc(i,u,i_ref,wp):
	"""Calculate derivate of 1st order filter - imag component"""
	duI = (wp)*(-u.imag +	i_ref.imag - i.imag)
	
	return duI

@numba.njit(cache=cache_flag)
def dxR_calc(u,Ki_GCC):
	"""Calculate derivative of current controller integral state - real component"""
	dxR = Ki_GCC*u.real
	
	return dxR

@numba.njit(cache=cache_flag)
def dxI_calc(u,Ki_GCC):
	"""Calculate derivative of current controller integral state - imag component"""
	dxI = Ki_GCC*u.imag
	
	return dxI

@numba.njit(cache=cache_flag)
def dVdc_calc(Vdc,Ppv,S,C):
	"""Calculate derivative of Vdc"""
	dVdc = (Ppv - S.real)/(Vdc*C)
	
	return dVdc

@numba.njit(cache=cache_flag)
def dxDC_calc(Vdc,Vdc_ref,Ki_DC):
	"""Calculate derivative of xDC"""
	dxDC = Ki_DC*(Vdc_ref - Vdc)
	
	return dxDC

@numba.njit(cache=cache_flag)
def dxQ_calc(S,Q_ref,Ki_Q):
	"""Calculate derivative of xQ"""
	dxQ = -Ki_Q*(Q_ref - S.imag)
	
	return dxQ

@numba.njit(cache=cache_flag)
def m_abs_calc(u,x,Kp_GCC):
	"""Calculate absolute value of duty cycle"""
	m_abs = abs(Kp_GCC*u + x)
	
	return m_abs


@numba.njit(cache=cache_flag)
def du_calc(u,x,i,i_ref,wp,Kp_GCC,m_limit):
	"""Calculate absolute value of duty cycle"""
	if m_abs_calc(u,x,Kp_GCC) > m_limit:
		if np.sign(duR_calc(i,u,i_ref,wp)) == np.sign(u.real):
			duR = 0.0
		else:
			duR = duR_calc(i,u,i_ref,wp)
				
		if np.sign(duI_calc(i,u,i_ref,wp)) == np.sign(u.imag):
			duI = 0.0
		else:
			duI = duI_calc(i,u,i_ref,wp)
				
	else:
		duR = duR_calc(i,u,i_ref,wp)
		duI = duI_calc(i,u,i_ref,wp)
		
	return duR,duI

@numba.njit(cache=cache_flag)
def dx_logic(u,x,Ki_GCC,Kp_GCC,m_limit):
	"""Calculate absolute value of duty cycle"""
	
	if m_abs_calc(u,x,Kp_GCC)>m_limit:
		if np.sign(dxR_calc(u,Ki_GCC)) == np.sign(x.real):
			dxR = 0.0
		else:
			#dxR = Ki_GCC*u.real
			dxR = dxR_calc(u,Ki_GCC)
		if np.sign(dxI_calc(u,Ki_GCC)) == np.sign(x.imag):
			dxI = 0.0
		else:
			#dxbI = Ki_GCC*ub.imag
			dxI = dxI_calc(u,Ki_GCC)
				
	else:
		#dxbR = Ki_GCC*ub.real
		#dxbI = Ki_GCC*ub.imag
		dxR = dxR_calc(u,Ki_GCC)
		dxI = dxI_calc(u,Ki_GCC)
	
	return dxR,dxI

@numba.njit(cache=cache_flag)
def dxDC_logic(xDC,xQ,Vdc,S,Vdc_ref,Q_ref,Ki_DC,Kp_DC,Kp_Q,iref_limit):
	"""Calculate absolute value of duty cycle"""
	
	if abs(utility_functions_numba.ia_ref_Vdc_Q_control_calc(xDC,xQ,Vdc,S,Vdc_ref,Q_ref,Kp_DC,Kp_Q))>iref_limit:
		if np.sign(dxDC_calc(Vdc,Vdc_ref,Ki_DC)) == np.sign(xDC):
			dxDC = 0.0
		else:
			dxDC = dxDC_calc(Vdc,Vdc_ref,Ki_DC)
	else:
		#dxDC = Ki_DC*(Vdc_ref - Vdc)
		dxDC = dxDC_calc(Vdc,Vdc_ref,Ki_DC)
	
	return dxDC

@numba.njit(cache=cache_flag)
def dxQ_logic(xDC,xQ,Vdc,S,Vdc_ref,Q_ref,Ki_Q,Kp_DC,Kp_Q,iref_limit):
	"""Reactive power controller dynamics"""
	
	if abs(utility_functions_numba.ia_ref_Vdc_Q_control_calc(xDC,xQ,Vdc,S,Vdc_ref,Q_ref,Kp_DC,Kp_Q))>iref_limit:
		if np.sign(dxQ_calc(S,Q_ref,Ki_Q)) == np.sign(xQ):
			dxQ = 0.0
		else:
			dxQ = dxQ_calc(S,Q_ref,Ki_Q)
	else:
		#dxQ = -Ki_Q*(Q_ref - S_PCC.imag)
		dxQ = dxQ_calc(S,Q_ref,Ki_Q)
	
	return dxQ

@numba.njit(cache=cache_flag)
def dxPLL_calc(vd,Ki_PLL):
	#Calculate absolute value of duty cycle

	dxPLL = Ki_PLL*vd
	
	return dxPLL

class SolarPVDERThreePhaseNumba(PVModule,SolarPVDER):
	"""
		Class for describing a Solar Photo-voltaic Distributed Energy Resource consisting of panel, converters, and
		control systems.
		Attributes:
		count (int): Number of instances of `SolarPVDERThreePhase`.
		n_ODE (int): Number of ODE's.
	"""
	count = 0 #Object count

	def __init__(self,events,configFile=None,**kwargs):		
		"""Creates an instance of `SolarPV_DER_ThreePhase`.
		Args:
			events (SimulationEvents): An instance of `SimulationEvents`.
			gridModel (Grid): An instance of `Gridl`(only need to be suppled for stand alone simulation).
			powerRating (float): A scalar specifying the rated power (VA) of the DER.
			VrmsRating (float): A scalar specifying the rated RMS L-G voltage (V) of the DER.
			ia0,xa0,ua0 (complex): Initial value of inverter states in p.u. of the DER instance.
			xDC0,xQ0,xPLL0,wte0 (float): Initial value of inverter states in the DER instance.
			gridVoltatePhaseA,gridVoltatePhaseA,gridVoltatePhaseA (float): Initial voltage phasor (V) at PCC - LV side from external program (only need to be suppled if model is not stand alone).
			standAlone (bool): Specify if the DER instance is a stand alone simulation or part of a larger simulation.
			STEADY_STATE_INITIALIZATION (bool): Specify whether states in the DER instance will be initialized to steady state values.
			allow_unbalanced_m (bool): Allow duty cycles to take on unbalanced values during initialization (default: False).
			derConfig (dict): Configuration parameters that may be supplied from an external program.
			identifier (str): An identifier that can be used to name the instance (default: None).
		Raises:
			ValueError: If parameters corresponding to `Sinverter_rated` are not available.
			ValueError: If rated DC link voltage is not sufficient.
		"""
		try:
			SolarPVDERThreePhaseNumba.count = SolarPVDERThreePhaseNumba.count+1 #Increment count to keep track of number of PV-DER model instances
			DER_arguments = self.setup_DER(events,configFile,**kwargs)		 
		
			if six.PY3:
				super().__init__(self.DER_config['basic_options']['Sinsol'])	#Initialize PV module class (base class)
			elif six.PY2:
				super(SolarPVDERThreePhaseNumba,self).__init__(self.DER_config['basic_options']['Sinsol'])
			self.initialize_DER(DER_arguments)
			self.creation_message()
		except:
			LogUtil.exception_handler()

	@property#Decorator used for auto updating
	def y0(self):
		"""List of initial states"""
		try:
			return	[self.ia.real, self.ia.imag, self.xa.real, self.xa.imag, self.ua.real,self.ua.imag,\
					 self.ib.real, self.ib.imag, self.xb.real, self.xb.imag, self.ub.real,self.ub.imag,\
					 self.ic.real, self.ic.imag, self.xc.real, self.xc.imag, self.uc.real,self.uc.imag,\
					 self.Vdc,self.xDC,self.xQ,self.xPLL,self.wte]
		except:
			LogUtil.exception_handler()

	def update_Ppv(self,t):
		"""Update PV module power output based on solar events and DC link voltage."""
		try:
			Sinsol_new,Tactual_new = self.events.solar_events(t)
			if abs(self.Sinsol- Sinsol_new) or abs(self.Tactual- Tactual_new) > 0.0:	#Update Iph only if solar insolation changes
				self.Sinsol = Sinsol_new
				self.Tactual = Tactual_new
				utility_functions.print_to_terminal("{}:PV module current output changed from {:.3f} A to {:.3f} A at {:.3f} s".format(self.name,self.Iph,self.Iph_calc(),t))
				self.Iph = Iph_calc(self.Tactual,self.Sinsol)
			
			self.Ppv = Ppv_calc(self.Tactual,self.Iph,self.Vdc_actual,self.Np,self.Ns,self.Sbase)
		except:
			LogUtil.exception_handler()

	def update_voltages(self):
		"""Update voltages."""
		try:
			#Update inverter terminal voltage
			ma = utility_functions_numba.m_calc(self.Kp_GCC,self.ua,self.xa)
			mb = utility_functions_numba.m_calc(self.Kp_GCC,self.ub,self.xb)
			mc = utility_functions_numba.m_calc(self.Kp_GCC,self.uc,self.xc)
			
			self.vta = utility_functions_numba.vt_calc(ma,self.Vdc)
			self.vtb = utility_functions_numba.vt_calc(mb,self.Vdc)
			self.vtc = utility_functions_numba.vt_calc(mc,self.Vdc)			
			
			#Update PCC LV side voltage
			self.va = self.va_calc()
			self.vb = self.vb_calc()
			self.vc = self.vc_calc()
		
			#self.vtunbalance = utility_functions_numba.Uunbalance_calc(self.vta,self.vtb,self.vtc)
			#self.vunbalance = utility_functions_numba.Uunbalance_calc(self.va,self.vb,self.vc)
		except:
			LogUtil.exception_handler()


	def update_RMS(self):
		"""Update RMS voltages."""
		try:
			
			self.Vtrms = utility_functions_numba.Urms_calc(self.vta,self.vtb,self.vtc)
			self.Vrms = utility_functions_numba.Urms_calc(self.va,self.vb,self.vc)
			self.Vrms_min = utility_functions_numba.Urms_min_calc(self.va,self.vb,self.vc)
			self.Irms = utility_functions_numba.Urms_calc(self.ia,self.ib,self.ic)

			#Update RMS values
			if self.DO_EXTRA_CALCULATIONS:
				self.Vtabrms = utility_functions_numba.Vabrms_calc(self.vta,self.vtb)
				self.Vabrms = utility_functions_numba.Vabrms_calc(self.va,self.vb)
		except:
			LogUtil.exception_handler()

	def update_power(self):
		"""Update RMS voltages."""
		try:
			#Update power output
			self.S = utility_functions_numba.S_calc(self.vta,self.vtb,self.vtc,self.ia,self.ib,self.ic)
			self.S_PCC = utility_functions_numba.S_calc(self.va,self.vb,self.vc,self.ia,self.ib,self.ic)
			
			if self.DO_EXTRA_CALCULATIONS:
				self.S_PCCa = self.S_PCCph_calc(self.va,self.ia)
				self.S_PCCb = self.S_PCCph_calc(self.vb,self.ib)
				self.S_PCCc = self.S_PCCph_calc(self.vc,self.ic)
		
			if self.standAlone:	#Update load current and grid voltage source power only in stand alone mode
				
				self.iaload1 = utility_functions_numba.i_load_calc(self.va,self.Zload1)
				self.ibload1 = utility_functions_numba.i_load_calc(self.vb,self.Zload1)
				self.icload1 = utility_functions_numba.i_load_calc(self.vc,self.Zload1)
				
				self.S_G = utility_functions_numba.S_G_calc(self.va,self.vb,self.vc,
															self.grid_model.vag,self.grid_model.vbg,self.grid_model.vcg,
															self.ia,self.ib,self.ic,self.Zload1,self.a)
				self.S_load1 = utility_functions_numba.S_load_calc(self.va,self.vb,self.vc,self.Zload1)
		except:
			LogUtil.exception_handler()

	def update_iref(self,t):
		"""Update reference reference current."""
		try:
			#Get current controller setpoint
			if not self.current_gradient_limiter:
				self.ia_ref = utility_functions_numba.ia_ref_Vdc_Q_control_calc(self.xDC,self.xQ,self.Vdc,self.S_PCC,self.Vdc_ref,self.Q_ref,self.Kp_DC,self.Kp_Q)				  
			else:
				self.ia_ref = self.get_ramp_limited_iref(t,utility_functions_numba.ia_ref_Vdc_Q_control_calc(self.xDC,self.xQ,self.Vdc,self.S_PCC,self.Vdc_ref,self.Q_ref,self.Kp_DC,self.Kp_Q))
			self.t_iref = t	
			self.ib_ref = utility_functions_numba.Ub_calc(self.ia_ref)
			self.ic_ref = utility_functions_numba.Uc_calc(self.ia_ref)
			
		except:
			LogUtil.exception_handler()
	def update_inverter_frequency(self,t):
		"""Update d-q quantities."""
		try:
			self.wgrid_measured = self.wgrid_calc(t) #Update grid frequency
			#d-q transformation
			#Convert PCC LV side voltage from phasor to time domain
			self.vat,self.vbt,self.vct = utility_functions_numba.phasor_to_time(upha = self.va,uphb = self.vb,uphc = self.vc,w=self.wgrid_measured,t=t)
		
			#Convert from 3ph time domain to d-q using Parks transformation
			self.vd,self.vq,self.v0 = utility_functions_numba.abc_to_dq0(self.vat,self.vbt,self.vct,self.wte) #PCC LV side voltage
			
			self.we = utility_functions_numba.we_calc(self.xPLL,self.vd,self.Kp_PLL)
			self.winv = self.we
		except:
			LogUtil.exception_handler()


	def ODE_model(self,y,t):
		"""Derivatives for the equation."""
		try:
			"""
			iaR, iaI, xaR, xaI, uaR, uaI,\
			ibR, ibI, xbR, xbI, ubR, ubI,\
			icR, icI, xcR, xcI, ucR, ucI,\
			Vdc, xDC, xQ, xPLL, wte = y	 # unpack current values of y
			"""
			
			self.ia,self.xa,self.ua,\
			self.ib,self.xb,self.ub,\
			self.ic,self.xc,self.uc,\
			self.Vdc,self.xDC,self.xQ,\
			self.xPLL,self.wte = update_inverter_states( y[iaR_id],y[iaI_id],y[xaR_id],y[xaI_id],y[uaR_id],y[uaI_id],
														y[ibR_id],y[ibI_id],y[xbR_id],y[xbI_id],y[ubR_id],y[ubI_id],
														y[icR_id],y[icI_id],y[xcR_id],y[xcI_id],y[ucR_id],y[ucI_id],
														y[Vdc_id],y[xDC_id],y[xQ_id],y[xPLL_id],y[wte_id])
			self.update_Ppv(t)
			self.update_Zload1(t)
		
			self.update_voltages()
			self.update_power()
			self.update_RMS()
		
			self.update_Qref(t)
			self.update_Vdc_ref(t)	
			self.update_iref(t)
		
			self.update_inverter_frequency(t)
		
			self.update_ridethrough_flags(t)
			self.disconnect_or_reconnect(t)
		
			#Phase a inverter output current
			diaR =  diR_calc(self.ia,self.vta,self.va,self.winv,self.Rf,self.Lf,self.wbase)
			diaI =  diI_calc(self.ia,self.vta,self.va,self.winv,self.Rf,self.Lf,self.wbase)
			#Current controller dynamics
			dxaR,dxaI = dx_logic(self.ua,self.xa,self.Ki_GCC,self.Kp_GCC,self.m_limit)
			duaR,duaI= du_calc(self.ua,self.xa,self.ia,self.ia_ref,self.wp,self.Kp_GCC,self.m_limit)
			
			#Phase b inverter output current
			dibR =  diR_calc(self.ib,self.vtb,self.vb,self.winv,self.Rf,self.Lf,self.wbase)
			dibI =  diI_calc(self.ib,self.vtb,self.vb,self.winv,self.Rf,self.Lf,self.wbase)
			#Current controller dynamics - Phase b
			#dxbR = self.Ki_GCC*self.ub.real
			#dxbI = self.Ki_GCC*self.ub.imag
			#dubR = (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)
			#dubI = (self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)
			dxbR,dxbI = dx_logic(self.ub,self.xb,self.Ki_GCC,self.Kp_GCC,self.m_limit)	
			dubR,dubI= du_calc(self.ub,self.xb,self.ib,self.ib_ref,self.wp,self.Kp_GCC,self.m_limit)
			
			#Phase c inverter output current
			dicR =  diR_calc(self.ic,self.vtc,self.vc,self.winv,self.Rf,self.Lf,self.wbase)
			dicI =  diI_calc(self.ic,self.vtc,self.vc,self.winv,self.Rf,self.Lf,self.wbase)
			#Current controller dynamics - Phase c
			#dxcR = self.Ki_GCC*self.uc.real
			#dxcI = self.Ki_GCC*self.uc.imag
			#ducR = (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)
			#ducI = (self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)
			dxcR,dxcI = dx_logic(self.uc,self.xc,self.Ki_GCC,self.Kp_GCC,self.m_limit)	
			ducR,ducI= du_calc(self.uc,self.xc,self.ic,self.ic_ref,self.wp,self.Kp_GCC,self.m_limit)
			
			#DC link voltage dynamics
			#dVdc = (self.Ppv - self.S.real)/(self.Vdc*self.C)
			dVdc =  dVdc_calc(self.Vdc,self.Ppv,self.S,self.C)
			
			dxDC = dxDC_logic(self.xDC,self.xQ,self.Vdc,self.S_PCC,self.Vdc_ref,self.Q_ref,self.Ki_DC,self.Kp_DC,self.Kp_Q,self.iref_limit)
			dxQ = dxQ_logic(self.xDC,self.xQ,self.Vdc,self.S_PCC,self.Vdc_ref,self.Q_ref,self.Ki_Q,self.Kp_DC,self.Kp_Q,self.iref_limit)
			
			#SRF-PLL dynamics
			#dxPLL = self.Ki_PLL*(self.vd)
			dxPLL = dxPLL_calc(self.vd,self.Ki_PLL)
		
			#Frequency integration to get angle
			dwte = self.we
		
			result =	 [ diaR,# list of dy/dt=f functions
							 diaI,
							 dxaR,
							 dxaI,
							 duaR,
							 duaI,
							 dibR,
							 dibI,
							 dxbR,
							 dxbI,
							 dubR,
							 dubI,
							 dicR,
							 dicI,
							 dxcR,
							 dxcI,
							 ducR,
							 ducI,
							 dVdc,
							 dxDC,
							 dxQ,
							 dxPLL,
							 dwte]

			return np.array(result)
		except:
			LogUtil.exception_handler()


	def jac_ODE_model(self,y,t):
		"""Jacobian for the system of ODE's."""
		try:
			"""
			iaR, iaI, xaR, xaI, uaR, uaI,\
			ibR, ibI, xbR, xbI, ubR, ubI,\
			icR, icI, xcR, xcI, ucR, ucI,\
			Vdc, xDC, xQ, xPLL, wte = y	 # unpack current values of y
			
			self.update_inverter_states(iaR + 1j*iaI, xaR + 1j*xaI,uaR + 1j*uaI,\
										ibR + 1j*ibI, xbR + 1j*xbI,ubR + 1j*ubI,\
										icR + 1j*icI, xcR + 1j*xcI,ucR + 1j*ucI,\
										Vdc,xDC,xQ,\
										xPLL,wte)
			"""
			self.ia,self.xa,self.ua,\
			self.ib,self.xb,self.ub,\
			self.ic,self.xc,self.uc,\
			self.Vdc,self.xDC,self.xQ,\
			self.xPLL,self.wte = update_inverter_states( y[iaR_id],y[iaI_id],y[xaR_id],y[xaI_id],y[uaR_id],y[uaI_id],
														y[ibR_id],y[ibI_id],y[xbR_id],y[xbI_id],y[ubR_id],y[ubI_id],
														y[icR_id],y[icI_id],y[xcR_id],y[xcI_id],y[ucR_id],y[ucI_id],
														y[Vdc_id],y[xDC_id],y[xQ_id],y[xPLL_id],y[wte_id])

			J = self.J
			varInd = self.varInd 
			self.update_Ppv(t)
		
			self.update_voltages()
			self.update_power()
			self.update_RMS()
		
			self.update_Qref(t)
			self.update_Vdc_ref(t)	
			self.update_iref(t)
		
			#d-q transformation
			self.update_inverter_frequency(t)
		
			self.update_ridethrough_flags(t)
			self.disconnect_or_reconnect(t)
		
			#Phase a inverter output current
		
			ra,theta_a = cmath.polar(self.va)
			rb,theta_b = cmath.polar(self.vb)
			rc,theta_c = cmath.polar(self.vc)
			theta_a = self.wgrid_measured*t + theta_a - math.pi/2
			theta_b = self.wgrid_measured*t + theta_b - math.pi/2
			theta_c = self.wgrid_measured*t + theta_c - math.pi/2
			
			J[varInd['iaR'],varInd['iaR']] = -self.Rf/self.Lf			
			J[varInd['iaR'],varInd['iaI']] = (self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
			J[varInd['iaR'],varInd['xaR']] = self.Vdc/(2*self.Lf)
			J[varInd['iaR'],varInd['uaR']] = (self.Vdc*self.Kp_GCC)/(2*self.Lf)
			J[varInd['iaR'],varInd['Vdc']] = (self.xa.real+self.ua.real*self.Kp_GCC)/(2*self.Lf)
			J[varInd['iaR'],varInd['xPLL']] = self.ia.imag/self.wbase
			J[varInd['iaR'],varInd['wte']] = (2/3)*((self.Kp_PLL*self.ia.imag)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte)
											 + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
											 + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
		
			J[varInd['iaI'],varInd['iaR']]= -(self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
			J[varInd['iaI'],varInd['iaI']]= -self.Rf/self.Lf
			J[varInd['iaI'],varInd['xaI']]= self.Vdc/(2*self.Lf) 
			J[varInd['iaI'],varInd['uaI']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
			J[varInd['iaI'],varInd['Vdc']]= (self.xa.imag+self.ua.imag*self.Kp_GCC)/(2*self.Lf)
			J[varInd['iaI'],varInd['xPLL']]= -self.ia.real/self.wbase
			J[varInd['iaI'],varInd['wte']] = -(2/3)*((self.Kp_PLL*self.ia.real)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte) 
											 + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
											 + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
			
			#Current controller dynamics
			if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit:
				if np.sign(self.Ki_GCC*self.ua.real) == np.sign(self.xa.real):
					J[varInd['xaR'],varInd['uaR']]=0.0
				else:
					J[varInd['xaR'],varInd['uaR']]=self.Ki_GCC
				if np.sign(self.Ki_GCC*self.ua.imag) == np.sign(self.xa.imag):
					J[varInd['xaI'],varInd['uaI']]=0.0
				else:
					J[varInd['xaI'],varInd['uaI']]=self.Ki_GCC
				
			else:
					J[varInd['xaR'],varInd['uaR']]=self.Ki_GCC
					J[varInd['xaI'],varInd['uaI']]=self.Ki_GCC
		
			if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit:
				if np.sign( (self.wp)*(-self.ua.real +	self.ia_ref.real - self.ia.real)) == np.sign(self.ua.real):
					J[varInd['uaR'],varInd['iaR']]= 0.0
					J[varInd['uaR'],varInd['uaR']]= 0.0
					J[varInd['uaR'],varInd['Vdc']]= 0.0
					J[varInd['uaR'],varInd['xDC']]= 0.0	 
				else:
					#duaR = (self.wp)*(-self.ua.real +	self.ia_ref.real - self.ia.real)
					J[varInd['uaR'],varInd['iaR']]= -self.wp
					J[varInd['uaR'],varInd['uaR']]= -self.wp
					J[varInd['uaR'],varInd['Vdc']]= -self.wp*self.Kp_DC
					J[varInd['uaR'],varInd['xDC']]= self.wp					
					
				if np.sign((self.wp)*(-self.ua.imag +	self.ia_ref.imag - self.ia.imag)) == np.sign(self.ua.imag):
					#duaI = 0.0
					J[varInd['uaI'],varInd['iaR']]= 0.0
					J[varInd['uaI'],varInd['iaI']]= 0.0
					J[varInd['uaI'],varInd['uaI']]= 0.0
					J[varInd['uaI'],varInd['xQ']]= 0.0
					
					J[varInd['uaI'],varInd['ibR']]= 0.0
					J[varInd['uaI'],varInd['ibI']]= 0.0
					J[varInd['uaI'],varInd['icR']]= 0.0
					J[varInd['uaI'],varInd['icI']]= 0.0
				else:
					#duaI = (self.wp)*(-self.ua.imag +	self.ia_ref.imag - self.ia.imag)
					J[varInd['uaI'],varInd['iaR']]= (self.Kp_Q*self.wp*self.va.imag/2)
					J[varInd['uaI'],varInd['iaI']]= -self.wp - (self.Kp_Q*self.wp*self.va.real/2)
					J[varInd['uaI'],varInd['uaI']]= -self.wp
					J[varInd['uaI'],varInd['xQ']]= self.wp
					
					J[varInd['uaI'],varInd['ibR']]= (self.Kp_Q*self.wp*self.vb.imag/2)
					J[varInd['uaI'],varInd['ibI']]= - (self.Kp_Q*self.wp*self.vb.real/2)
					J[varInd['uaI'],varInd['icR']]= (self.Kp_Q*self.wp*self.vc.imag/2)
					J[varInd['uaI'],varInd['icI']]= - (self.Kp_Q*self.wp*self.vc.real/2)
				
			else:
				#duaR = (self.wp)*(-self.ua.real +	self.ia_ref.real - self.ia.real)
				#duaI = (self.wp)*(-self.ua.imag +	self.ia_ref.imag - self.ia.imag)
				J[varInd['uaR'],varInd['iaR']]= -self.wp
				J[varInd['uaR'],varInd['uaR']]= -self.wp
				J[varInd['uaR'],varInd['Vdc']]= -self.wp*self.Kp_DC
				J[varInd['uaR'],varInd['xDC']]= self.wp	 
				
				J[varInd['uaI'],varInd['iaR']]= (self.Kp_Q*self.wp*self.va.imag/2)
				J[varInd['uaI'],varInd['iaI']]= -self.wp - (self.Kp_Q*self.wp*self.va.real/2)
				J[varInd['uaI'],varInd['uaI']]= -self.wp
				J[varInd['uaI'],varInd['xQ']]= self.wp
				
				J[varInd['uaI'],varInd['ibR']]= (self.Kp_Q*self.wp*self.vb.imag/2)
				J[varInd['uaI'],varInd['ibI']]= - (self.Kp_Q*self.wp*self.vb.real/2)
				J[varInd['uaI'],varInd['icR']]= (self.Kp_Q*self.wp*self.vc.imag/2)
				J[varInd['uaI'],varInd['icI']]= - (self.Kp_Q*self.wp*self.vc.real/2)
		
			#Phase b inverter output current
			#dibR = (1/self.Lf)*(-self.Rf*self.ib.real - self.vb.real + self.vtb.real) + (self.winv/self.wbase)*self.ib.imag 
			#dibI = (1/self.Lf)*(-self.Rf*self.ib.imag - self.vb.imag + self.vtb.imag) - (self.winv/self.wbase)*self.ib.real	
			
			J[varInd['ibR'],varInd['ibR']]= -self.Rf/self.Lf			
			J[varInd['ibR'],varInd['ibI']]= (self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
			J[varInd['ibR'],varInd['xbR']]= self.Vdc/(2*self.Lf)
			J[varInd['ibR'],varInd['ubR']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
			J[varInd['ibR'],varInd['Vdc']]= (self.xb.real+self.ub.real*self.Kp_GCC)/(2*self.Lf)
			J[varInd['ibR'],varInd['xPLL']]= self.ib.imag/self.wbase	 
			J[varInd['ibR'],varInd['wte']] = (2/3)*((self.Kp_PLL*self.ib.imag)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte)
											 + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
											 + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
			
			J[varInd['ibI'],varInd['ibR']]= -(self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
			J[varInd['ibI'],varInd['ibI']]= -self.Rf/self.Lf
			J[varInd['ibI'],varInd['xbI']]= self.Vdc/(2*self.Lf)
			J[varInd['ibI'],varInd['ubI']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
			J[varInd['ibI'],varInd['Vdc']]= (self.xb.imag+self.ub.imag*self.Kp_GCC)/(2*self.Lf)
			J[varInd['ibI'],varInd['xPLL']]= -self.ib.real/self.wbase
			J[varInd['ibI'],varInd['wte']] = -(2/3)*((self.Kp_PLL*self.ib.real)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte) 
											 + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
											 + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
		
			#Current controller dynamics - Phase b
		
			if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
				if np.sign(self.Ki_GCC*self.ub.real) == np.sign(self.xb.real):
					#dxbR = 0.0
					J[varInd['xbR'],varInd['ubR']]=0.0
				else:
					#dxbR = self.Ki_GCC*self.ub.real
					J[varInd['xbR'],varInd['ubR']]=self.Ki_GCC
				if np.sign(self.Ki_GCC*self.ub.imag) == np.sign(self.xb.imag):
					#dxbI = 0.0
					J[varInd['xbI'],varInd['ubI']]=0.0
				else:
					#dxbI = self.Ki_GCC*self.ub.imag
					J[varInd['xbI'],varInd['ubI']]=self.Ki_GCC

			else: 
				#dxbR = self.Ki_GCC*self.ub.real
				#dxbI = self.Ki_GCC*self.ub.imag
				J[varInd['xbR'],varInd['ubR']]=self.Ki_GCC
				J[varInd['xbI'],varInd['ubI']]=self.Ki_GCC
		
			if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
				if np.sign( (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)) == np.sign(self.ub.real):
					#dubR = 0.0
					J[varInd['ubR'],varInd['iaR']]= 0.0
					J[varInd['ubR'],varInd['iaI']]= 0.0 
					J[varInd['ubR'],varInd['ibR']]= 0.0
					J[varInd['ubR'],varInd['ibI']]= 0.0
					J[varInd['ubR'],varInd['ubR']]= 0.0
					J[varInd['ubR'],varInd['icR']]= 0.0
					J[varInd['ubR'],varInd['icI']]= 0.0
				
					J[varInd['ubR'],varInd['Vdc']]= 0.0
					J[varInd['ubR'],varInd['xDC']]= 0.0
					J[varInd['ubR'],varInd['xQ']]= 0.0
				else:
					#dubR = (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)
					J[varInd['ubR'],varInd['iaR']]= 0.866025403*(self.Kp_Q*self.wp*self.va.imag/2)
					J[varInd['ubR'],varInd['iaI']]= -0.866025403*(self.Kp_Q*self.wp*self.va.real/2)
					J[varInd['ubR'],varInd['ibR']]= -self.wp + 0.866025403*(self.Kp_Q*self.wp*self.vb.imag/2)
					J[varInd['ubR'],varInd['ibI']]= -0.866025403*(self.Kp_Q*self.wp*self.vb.real/2)
					J[varInd['ubR'],varInd['ubR']]= -self.wp
					J[varInd['ubR'],varInd['icR']]= 0.866025403*(self.Kp_Q*self.wp*self.vc.imag/2)
					J[varInd['ubR'],varInd['icI']]= -0.866025403*(self.Kp_Q*self.wp*self.vc.real/2)
				
					J[varInd['ubR'],varInd['Vdc']]= 0.5*self.wp*self.Kp_DC
					J[varInd['ubR'],varInd['xDC']]= -0.5*self.wp	
					J[varInd['ubR'],varInd['xQ']]= 0.866025403*self.wp
					
				if np.sign((self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)) == np.sign(self.ub.imag):
					#dubI = 0.0
					J[varInd['ubI'],varInd['iaR']]= 0.0
					J[varInd['ubI'],varInd['iaI']]= 0.0
					J[varInd['ubI'],varInd['ibR']]= 0.0
					J[varInd['ubI'],varInd['ibI']]= 0.0
					J[varInd['ubI'],varInd['ubI']]= 0.0
					J[varInd['ubI'],varInd['icR']]= 0.0
					J[varInd['ubI'],varInd['icI']]= 0.0
				
					J[varInd['ubI'],varInd['Vdc']]= 0.0
					J[varInd['ubI'],varInd['xDC']]= 0.0
					J[varInd['ubI'],varInd['xQ']]= 0.0
				
				else:
					#dubI = (self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)
					J[varInd['ubI'],varInd['iaR']]= -(self.Kp_Q*self.wp*self.va.imag/4)
					J[varInd['ubI'],varInd['iaI']]= (self.Kp_Q*self.wp*self.va.real/4)
					J[varInd['ubI'],varInd['ibR']]= -(self.Kp_Q*self.wp*self.vb.imag/4)
					J[varInd['ubI'],varInd['ibI']]= -self.wp + (self.Kp_Q*self.wp*self.vb.real/4)
					J[varInd['ubI'],varInd['ubI']]= -self.wp				
					J[varInd['ubI'],varInd['icR']]= -(self.Kp_Q*self.wp*self.vc.imag/4)
					J[varInd['ubI'],varInd['icI']]= (self.Kp_Q*self.wp*self.vc.real/4)
				
					J[varInd['ubI'],varInd['Vdc']]= 0.866025403*self.Kp_DC*self.wp
					J[varInd['ubI'],varInd['xDC']]= -0.866025403*self.wp
					J[varInd['ubI'],varInd['xQ']]= -0.5*self.wp
					
			else:
				#dubR = (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)
				#dubI = (self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)
			
				J[varInd['ubR'],varInd['iaR']]= 0.866025403*(self.Kp_Q*self.wp*self.va.imag/2)
				J[varInd['ubR'],varInd['iaI']]= -0.866025403*(self.Kp_Q*self.wp*self.va.real/2)
				J[varInd['ubR'],varInd['ibR']]= -self.wp + 0.866025403*(self.Kp_Q*self.wp*self.vb.imag/2)
				J[varInd['ubR'],varInd['ibI']]= -0.866025403*(self.Kp_Q*self.wp*self.vb.real/2)
				J[varInd['ubR'],varInd['ubR']]= -self.wp
				J[varInd['ubR'],varInd['icR']]= 0.866025403*(self.Kp_Q*self.wp*self.vc.imag/2)
				J[varInd['ubR'],varInd['icI']]= -0.866025403*(self.Kp_Q*self.wp*self.vc.real/2)
			
				J[varInd['ubR'],varInd['Vdc']]= 0.5*self.wp*self.Kp_DC
				J[varInd['ubR'],varInd['xDC']]= -0.5*self.wp	
				J[varInd['ubR'],varInd['xQ']]= 0.866025403*self.wp
					
				J[varInd['ubI'],varInd['iaR']]= -(self.Kp_Q*self.wp*self.va.imag/4)
				J[varInd['ubI'],varInd['iaI']]= (self.Kp_Q*self.wp*self.va.real/4)
				J[varInd['ubI'],varInd['ibR']]= -(self.Kp_Q*self.wp*self.vb.imag/4)
				J[varInd['ubI'],varInd['ibI']]= -self.wp + (self.Kp_Q*self.wp*self.vb.real/4)
				J[varInd['ubI'],varInd['ubI']]= -self.wp				
				J[varInd['ubI'],varInd['icR']]= -(self.Kp_Q*self.wp*self.vc.imag/4)
				J[varInd['ubI'],varInd['icI']]= (self.Kp_Q*self.wp*self.vc.real/4)
				
				J[varInd['ubI'],varInd['Vdc']]= 0.866025403*self.Kp_DC*self.wp
				J[varInd['ubI'],varInd['xDC']]= -0.866025403*self.wp
				J[varInd['ubI'],varInd['xQ']]= -0.5*self.wp					
			
			#Phase c inverter output current
			#dicR = (1/self.Lf)*(-self.Rf*self.ic.real - self.vc.real + self.vtc.real) + (self.winv/self.wbase)*self.ic.imag 
			#dicI = (1/self.Lf)*(-self.Rf*self.ic.imag - self.vc.imag + self.vtc.imag) - (self.winv/self.wbase)*self.ic.real 
		
			J[varInd['icR'],varInd['icR']]= -self.Rf/self.Lf			
			J[varInd['icR'],varInd['icI']]= (self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
			J[varInd['icR'],varInd['xcR']]= self.Vdc/(2*self.Lf)
			J[varInd['icR'],varInd['ucR']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
			J[varInd['icR'],varInd['Vdc']]= (self.xc.real+self.uc.real*self.Kp_GCC)/(2*self.Lf)
			J[varInd['icR'],varInd['xPLL']]= self.ic.imag/self.wbase 
			J[varInd['icR'],varInd['wte']] = (2/3)*((self.Kp_PLL*self.ic.imag)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte)
												+ 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
												+ 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
			
			J[varInd['icI'],varInd['icR']]= -(self.xPLL+self.Kp_PLL*self.vd+2*math.pi*60)/self.wbase
			J[varInd['icI'],varInd['icI']]= -self.Rf/self.Lf
			J[varInd['icI'],varInd['xcI']]= self.Vdc/(2*self.Lf)
			J[varInd['icI'],varInd['ucI']]= (self.Vdc*self.Kp_GCC)/(2*self.Lf)
			J[varInd['icI'],varInd['Vdc']]= (self.xc.imag+self.uc.imag*self.Kp_GCC)/(2*self.Lf)
			J[varInd['icI'],varInd['xPLL']]= -self.ic.real/self.wbase
			J[varInd['icI'],varInd['wte']] = -(2/3)*((self.Kp_PLL*self.ic.real)/self.wbase)*(-ra*math.cos(theta_a)*math.sin(self.wte) 
											 + 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)
											 + 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
			
			#Current controller dynamics - Phase c
		
			if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
				if np.sign(self.Ki_GCC*self.uc.real) == np.sign(self.xc.real):
					#dxcR = 0.0
					J[varInd['xcR'],varInd['ucR']]=0.0
				else:
					#dxcR = self.Ki_GCC*self.uc.real
					J[varInd['xcR'],varInd['ucR']]=self.Ki_GCC
				if np.sign(self.Ki_GCC*self.uc.imag) == np.sign(self.xc.imag):
					#dxcI = 0.0
					J[varInd['xcI'],varInd['ucI']]=0.0
				else:
					#dxcI = self.Ki_GCC*self.uc.imag
					J[varInd['xcI'],varInd['ucI']]=self.Ki_GCC
			else: 
				#dxcR = self.Ki_GCC*self.uc.real
				#dxcI = self.Ki_GCC*self.uc.imag
				J[varInd['xcR'],varInd['ucR']]=self.Ki_GCC
				J[varInd['xcI'],varInd['ucI']]=self.Ki_GCC
		
			if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
				if np.sign( (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)) == np.sign(self.uc.real):
					#ducR = 0.0
					J[varInd['ucR'],varInd['iaR']]= 0.0
					J[varInd['ucR'],varInd['iaI']]= 0.0
					J[varInd['ucR'],varInd['ibR']]= 0.0
					J[varInd['ucR'],varInd['ibI']]= 0.0
					J[varInd['ucR'],varInd['icR']]= 0.0
					J[varInd['ucR'],varInd['icI']]= 0.0
					J[varInd['ucR'],varInd['ucR']]= 0.0
					J[varInd['ucR'],varInd['Vdc']]= 0.0
					J[varInd['ucR'],varInd['xDC']]= 0.0
					J[varInd['ucR'],varInd['xQ']]= 0.0
				else:
					#ducR = (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)
					J[varInd['ucR'],varInd['iaR']]= -0.866025403*(self.Kp_Q*self.wp*self.va.imag/2)
					J[varInd['ucR'],varInd['iaI']]= 0.866025403*(self.Kp_Q*self.wp*self.va.real/2)
					J[varInd['ucR'],varInd['ibR']]= -0.866025403*(self.Kp_Q*self.wp*self.vb.imag/2)
					J[varInd['ucR'],varInd['ibI']]= 0.866025403*(self.Kp_Q*self.wp*self.vb.real/2)
					J[varInd['ucR'],varInd['icR']]= -self.wp	-0.866025403*(self.Kp_Q*self.wp*self.vc.imag/2)
					J[varInd['ucR'],varInd['icI']]= 0.866025403*(self.Kp_Q*self.wp*self.vc.real/2)
					J[varInd['ucR'],varInd['ucR']]= -self.wp
				
					J[varInd['ucR'],varInd['Vdc']]=	0.5*self.wp*self.Kp_DC
					J[varInd['ucR'],varInd['xDC']]=	-0.5*self.wp
					J[varInd['ucR'],varInd['xQ']]=	-0.866025403*self.wp
					
				if np.sign((self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)) == np.sign(self.uc.imag):
					#ducI = 0.0
					J[varInd['ucI'],varInd['icR']]= 0.0
					J[varInd['ucI'],varInd['icI']]= 0.0
					J[varInd['ucI'],varInd['ucI']]= 0.0
					
					J[varInd['ucI'],varInd['iaR']]= 0.0
					J[varInd['ucI'],varInd['iaI']]= 0.0
					J[varInd['ucI'],varInd['ibR']]= 0.0
					J[varInd['ucI'],varInd['ibI']]= 0.0
					
					J[varInd['ucI'],varInd['Vdc']]= 0.0
					J[varInd['ucI'],varInd['xDC']]= 0.0
					J[varInd['ucI'],varInd['xQ']]= 0.0
					
				else:
					#ducI = (self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)
					J[varInd['ucI'],varInd['iaR']]= -self.Kp_Q*self.wp*self.va.imag/4
					J[varInd['ucI'],varInd['iaI']]= self.Kp_Q*self.wp*self.va.real/4
					J[varInd['ucI'],varInd['ibR']]= -self.Kp_Q*self.wp*self.vb.imag/4
					J[varInd['ucI'],varInd['ibI']]= self.Kp_Q*self.wp*self.vb.real/4
					J[varInd['ucI'],varInd['icR']]= -(self.Kp_Q*self.wp*self.vc.imag/4)
					J[varInd['ucI'],varInd['icI']]= -self.wp + (self.Kp_Q*self.wp*self.vc.real/4)
					J[varInd['ucI'],varInd['ucI']]= -self.wp
					
					J[varInd['ucI'],varInd['Vdc']]= -0.8660254037*self.Kp_DC*self.wp
					J[varInd['ucI'],varInd['xDC']]= 0.8660254037*self.wp
					J[varInd['ucI'],varInd['xQ']]= -0.5*self.wp					
			else:
				#ducR = (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)
				#ducI = (self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)
				J[varInd['ucR'],varInd['iaR']]= -0.866025403*(self.Kp_Q*self.wp*self.va.imag/2)
				J[varInd['ucR'],varInd['iaI']]= 0.866025403*(self.Kp_Q*self.wp*self.va.real/2)
				J[varInd['ucR'],varInd['ibR']]= -0.866025403*(self.Kp_Q*self.wp*self.vb.imag/2)
				J[varInd['ucR'],varInd['ibI']]= 0.866025403*(self.Kp_Q*self.wp*self.vb.real/2)
				J[varInd['ucR'],varInd['icR']]= -self.wp	-0.866025403*(self.Kp_Q*self.wp*self.vc.imag/2)
				J[varInd['ucR'],varInd['icI']]= 0.866025403*(self.Kp_Q*self.wp*self.vc.real/2)
				J[varInd['ucR'],varInd['ucR']]= -self.wp
				
				J[varInd['ucR'],varInd['Vdc']]=	0.5*self.wp*self.Kp_DC
				J[varInd['ucR'],varInd['xDC']]=	-0.5*self.wp
				J[varInd['ucR'],varInd['xQ']]=	-0.866025403*self.wp
				
				J[varInd['ucI'],varInd['iaR']]= -self.Kp_Q*self.wp*self.va.imag/4
				J[varInd['ucI'],varInd['iaI']]= self.Kp_Q*self.wp*self.va.real/4
				J[varInd['ucI'],varInd['ibR']]= -self.Kp_Q*self.wp*self.vb.imag/4
				J[varInd['ucI'],varInd['ibI']]= self.Kp_Q*self.wp*self.vb.real/4
				J[varInd['ucI'],varInd['icR']]= -(self.Kp_Q*self.wp*self.vc.imag/4)
				J[varInd['ucI'],varInd['icI']]= -self.wp + (self.Kp_Q*self.wp*self.vc.real/4)
				J[varInd['ucI'],varInd['ucI']]= -self.wp
				
				J[varInd['ucI'],varInd['Vdc']]= -0.8660254037*self.Kp_DC*self.wp
				J[varInd['ucI'],varInd['xDC']]= 0.8660254037*self.wp
				J[varInd['ucI'],varInd['xQ']]= -0.5*self.wp
		
			#DC link voltage dynamics
			dVdc = (self.Ppv - self.S.real)/(self.Vdc*self.C)
			J[varInd['Vdc'],varInd['iaR']]= -(self.xa.real+self.Kp_GCC*self.ua.real)/(4*self.C)
			J[varInd['Vdc'],varInd['iaI']]= -(self.xa.imag+self.Kp_GCC*self.ua.imag)/(4*self.C)
			J[varInd['Vdc'],varInd['xaR']]= -self.ia.real/(4*self.C)
			J[varInd['Vdc'],varInd['xaI']]= -self.ia.imag/(4*self.C)
			J[varInd['Vdc'],varInd['uaR']]= -(self.Kp_GCC*self.ia.real)/(4*self.C)
			J[varInd['Vdc'],varInd['uaI']]= -(self.Kp_GCC*self.ia.imag)/(4*self.C)

			J[varInd['Vdc'],varInd['ibR']]= -(self.xb.real+self.Kp_GCC*self.ub.real)/(4*self.C)
			J[varInd['Vdc'],varInd['ibI']]= -(self.xb.imag+self.Kp_GCC*self.ub.imag)/(4*self.C)
			J[varInd['Vdc'],varInd['xbR']]= -self.ib.real/(4*self.C)
			J[varInd['Vdc'],varInd['xbI']]= -self.ib.imag/(4*self.C)
			J[varInd['Vdc'],varInd['ubR']]= -(self.Kp_GCC*self.ib.real)/(4*self.C)
			J[varInd['Vdc'],varInd['ubI']]= -(self.Kp_GCC*self.ib.imag)/(4*self.C)

			J[varInd['Vdc'],varInd['icR']]= -(self.xc.real+self.Kp_GCC*self.uc.real)/(4*self.C)
			J[varInd['Vdc'],varInd['icI']]= -(self.xc.imag+self.Kp_GCC*self.uc.imag)/(4*self.C)
			J[varInd['Vdc'],varInd['xcR']]= -self.ic.real/(4*self.C)
			J[varInd['Vdc'],varInd['xcI']]= -self.ic.imag/(4*self.C)
			J[varInd['Vdc'],varInd['ucR']]= -(self.Kp_GCC*self.ic.real)/(4*self.C)
			J[varInd['Vdc'],varInd['ucI']]= -(self.Kp_GCC*self.ic.imag)/(4*self.C)				
			
			J[varInd['Vdc'],varInd['Vdc']]= (-(self.q*self.Np*self.Irs*(self.Vdcbase**2))/(self.C*self.k*self.A*self.Ns*self.Tactual*self.Sbase))*math.exp((self.q*self.Vdc*self.Vdcbase)/(self.k*self.A*self.Ns*self.Tactual))
			#DC link voltage controller dynamics
			if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ	- self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
				if np.sign(self.Ki_DC*(self.Vdc_ref - self.Vdc)) == np.sign(self.xDC):
					#dxDC = 0.0
					J[varInd['xDC'],varInd['Vdc']]= 0.0
				else:
					#dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
					J[varInd['xDC'],varInd['Vdc']]=-self.Ki_DC
			else:
				#dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
				J[varInd['xDC'],varInd['Vdc']]=-self.Ki_DC
			
			# Reactive power controller dynamics
			if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ	- self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
			
				if np.sign(-self.Ki_Q*(self.Q_ref - self.S_PCC.imag)) == np.sign(self.xQ):
					#dxQ = 0.0
					J[varInd['xQ'],varInd['iaR']]= 0.0
					J[varInd['xQ'],varInd['iaI']]= 0.0
					J[varInd['xQ'],varInd['ibR']]= 0.0
					J[varInd['xQ'],varInd['ibI']]= 0.0
					J[varInd['xQ'],varInd['icR']]= 0.0
					J[varInd['xQ'],varInd['icI']]= 0.0
				else:
					#dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
					J[varInd['xQ'],varInd['iaR']]= (self.Ki_Q*self.va.imag/2)
					J[varInd['xQ'],varInd['iaI']]= -(self.Ki_Q*self.va.real/2)
					J[varInd['xQ'],varInd['ibR']]= (self.Ki_Q*self.vb.imag/2)
					J[varInd['xQ'],varInd['ibI']]= -(self.Ki_Q*self.vb.real/2)
					J[varInd['xQ'],varInd['icR']]= (self.Ki_Q*self.vc.imag/2)
					J[varInd['xQ'],varInd['icI']]= -(self.Ki_Q*self.vc.real/2)
		
			else:
				#dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
				J[varInd['xQ'],varInd['iaR']]= (self.Ki_Q*self.va.imag/2)
				J[varInd['xQ'],varInd['iaI']]= -(self.Ki_Q*self.va.real/2)
				J[varInd['xQ'],varInd['ibR']]= (self.Ki_Q*self.vb.imag/2)
				J[varInd['xQ'],varInd['ibI']]= -(self.Ki_Q*self.vb.real/2)
				J[varInd['xQ'],varInd['icR']]= (self.Ki_Q*self.vc.imag/2)
				J[varInd['xQ'],varInd['icI']]= -(self.Ki_Q*self.vc.real/2)
		
			#SRF-PLL dynamics
			#dxPLL = self.Ki_PLL*(self.vd)
			J[varInd['xPLL'],varInd['wte']] = (2/3)*self.Ki_PLL*(-ra*math.cos(theta_a)*math.sin(self.wte) +
												 0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)+
												 0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))
		
			#Frequency integration to get angle
			#dwte = self.we
			J[varInd['wte'],varInd['xPLL']]= 1
			J[varInd['wte'],varInd['wte']] = (2/3)*self.Kp_PLL*(-ra*math.cos(theta_a)*math.sin(self.wte) +
												0.5*rb*math.cos(theta_b)*math.sin(self.wte) + 0.8660254037*rb*math.cos(theta_b)*math.cos(self.wte)+
												0.5*rc*math.cos(theta_c)*math.sin(self.wte) - 0.8660254037*rc*math.cos(theta_c)*math.cos(self.wte))

			return J
		except:
			LogUtil.exception_handler()


