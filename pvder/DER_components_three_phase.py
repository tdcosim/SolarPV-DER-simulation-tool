"""Three phase PV-DER components."""

from __future__ import division
import numpy as np
import math
import cmath
import scipy
import six
import pdb
import warnings

from pvder.DER_components import SolarPVDER,PVModule
from pvder.DER_check_and_initialize import PVDER_SetupUtilities
from pvder.DER_features import PVDER_SmartFeatures
from pvder.DER_utilities import PVDER_ModelUtilities
from pvder.grid_components import BaseValues
from pvder import utility_functions
from pvder import defaults,templates
from pvder.logutil import LogUtil


class SolarPVDERThreePhase(PVModule,SolarPVDER):
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
			SolarPVDERThreePhase.count = SolarPVDERThreePhase.count+1 #Increment count to keep track of number of PV-DER model instances
			DER_arguments = self.setup_DER(events,configFile,**kwargs)		 
		
			if six.PY3:
				super().__init__(self.DER_config['basic_options']['Sinsol'])	#Initialize PV module class (base class)
			elif six.PY2:
				super(SolarPVDERThreePhase,self).__init__(self.DER_config['basic_options']['Sinsol'])
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
	
	#Apparent power output at inverter terminal
	def S_calc(self):
		"""Inverter apparent power output"""
		try:
			return (1/2)*(self.vta*self.ia.conjugate() + self.vtb*self.ib.conjugate() + self.vtc*self.ic.conjugate())*1.0
			#return utility_functions.S_calc(self.vta,self.vtb,self.vtc,self.ia,self.ib,self.ic)
		except:
			LogUtil.exception_handler()


	#Apparent power output at PCC - LV side
	def S_PCC_calc(self):
		"""Power output at PCC LV side"""
		try:
			return (1/2)*(self.va*self.ia.conjugate() + self.vb*self.ib.conjugate() + self.vc*self.ic.conjugate())*1.0
			#return utility_functions.S_calc(self.va,self.vb,self.vc,self.ia,self.ib,self.ic)
		except:
			LogUtil.exception_handler()


	def S_load1_calc(self):
		"""Power absorbed by load at PCC LV side."""
		try:
			return (1/2)*(self.va*(-(self.va/self.Zload1)).conjugate() + self.vb*(-(self.vb/self.Zload1)).conjugate() + self.vc*(-(self.vc/self.Zload1)).conjugate())
		except:
			LogUtil.exception_handler()


	def S_G_calc(self):
		"""Power absorbed/produced by grid voltage source."""
		try:
			return (1/2)*((-(self.ia-(self.va/self.Zload1))/self.a).conjugate()*self.grid_model.vag+(-(self.ib-(self.vb/self.Zload1))/self.a).conjugate()*self.grid_model.vbg+(-(self.ic-(self.vc/self.Zload1))/self.a).conjugate()*self.grid_model.vcg)
		except:
			LogUtil.exception_handler()


	#@property
	def Vtrms_calc(self):
		"""Inverter terminal voltage -	RMS"""
		try:
			return utility_functions.Urms_calc(self.vta,self.vtb,self.vtc)
		except:
			LogUtil.exception_handler()


	def Vrms_calc(self):
		"""PCC LV side voltage - RMS"""
		try:
			return utility_functions.Urms_calc(self.va,self.vb,self.vc)
		except:
			LogUtil.exception_handler()

	def Vrms_min_calc(self):
		"""PCC LV side voltage - RMS"""
		
		return utility_functions.Urms_min_calc(self.va,self.vb,self.vc)
		
	def Irms_calc(self):
		"""Inverter current - RMS"""
		try:
			return utility_functions.Urms_calc(self.ia,self.ib,self.ic)
		except:
			LogUtil.exception_handler()

	def Vtabrms_calc(self):
		"""Inverter terminal voltage - line to line	RMS"""
		
		return abs(self.vta-self.vtb)/math.sqrt(2)
	
	def Vabrms_calc(self):
		"""PCC LV side voltage - line to line	RMS"""
		try:
			return abs(self.va-self.vb)/math.sqrt(2)
		except:
			LogUtil.exception_handler()

	def update_inverter_states(self,ia,xa,ua,ib,xb,ub,ic,xc,uc,Vdc,xDC,xQ,xPLL,wte):
		"""Update inverter states"""
		try:
			self.ia = ia
			self.xa = xa
			self.ua = ua
		
			self.ib = ib
			self.xb = xb
			self.ub = ub
		
			self.ic = ic
			self.xc = xc
			self.uc = uc
		
			self.Vdc = Vdc
			self.xDC = xDC
			self.xQ = xQ
		
			self.xPLL = xPLL
			self.wte = wte
		except:
			LogUtil.exception_handler()


	def update_voltages(self):
		"""Update voltages."""
		try:
			#Update inverter terminal voltage
			self.vta = self.vta_calc()
			self.vtb = self.vtb_calc()
			self.vtc = self.vtc_calc()
				
			#Update PCC LV side voltage
			self.va = self.va_calc()
			self.vb = self.vb_calc()
			self.vc = self.vc_calc()
		
			#self.vtunbalance = utility_functions.Uunbalance_calc(self.vta,self.vtb,self.vtc)
			#self.vunbalance = utility_functions.Uunbalance_calc(self.va,self.vb,self.vc)
		except:
			LogUtil.exception_handler()


	def update_RMS(self):
		"""Update RMS voltages."""
		try:
			self.Vtrms = self.Vtrms_calc()
			self.Vrms = self.Vrms_calc()
			self.Vrms_min = self.Vrms_min_calc()
			self.Irms = self.Irms_calc()

			#Update RMS values
			if self.DO_EXTRA_CALCULATIONS:
				self.Vtabrms = self.Vtabrms_calc()
				self.Vabrms = self.Vabrms_calc()
		except:
			LogUtil.exception_handler()


	def update_power(self):
		"""Update RMS voltages."""
		try:
			#Update power output
			self.S = self.S_calc()
			self.S_PCC = self.S_PCC_calc()
		
			if self.DO_EXTRA_CALCULATIONS:
				self.S_PCCa = self.S_PCCph_calc(self.va,self.ia)
				self.S_PCCb = self.S_PCCph_calc(self.vb,self.ib)
				self.S_PCCc = self.S_PCCph_calc(self.vc,self.ic)
		
			if self.standAlone:	#Update load current and grid voltage source power only in stand alone mode
				self.iaload1 = self.iphload1_calc(self.va)
				self.ibload1 = self.iphload1_calc(self.vb)
				self.icload1 = self.iphload1_calc(self.vc)
				self.S_G = self.S_G_calc()
				self.S_load1 = self.S_load1_calc()		 
		except:
			LogUtil.exception_handler()


	def update_iref(self):
		"""Update reference reference current."""
		try:
			#Get current controller setpoint
			self.ia_ref = self.ia_ref_calc()
			self.ib_ref = self.ib_ref_calc()
			self.ic_ref = self.ic_ref_calc()
		except:
			LogUtil.exception_handler()


	def update_inverter_frequency(self,t):
		"""Update d-q quantities."""
		try:
			self.wgrid_measured = self.wgrid_calc(t) #Update grid frequency
			#d-q transformation		
			#Convert PCC LV side voltage from phasor to time domain
			self.vat,self.vbt,self.vct = utility_functions.phasor_to_time(upha = self.va,uphb = self.vb,uphc = self.vc,w=self.wgrid_measured,t=t)
		
			#Convert from 3ph time domain to d-q using Parks transformation
			self.vd,self.vq,self.v0 = utility_functions.abc_to_dq0(self.vat,self.vbt,self.vct,self.wte) #PCC LV side voltage
			self.we = self.we_calc() #Calculate inverter frequency from PLL equation
			self.winv = self.we#*self.wbase
		except:
			LogUtil.exception_handler()


	def ODE_model(self,y,t):
		"""Derivatives for the equation."""
		try:
			iaR, iaI, xaR, xaI, uaR, uaI,\
			ibR, ibI, xbR, xbI, ubR, ubI,\
			icR, icI, xcR, xcI, ucR, ucI,\
			Vdc, xDC, xQ, xPLL, wte = y	 # unpack current values of y
		
			self.update_inverter_states(iaR + 1j*iaI, xaR + 1j*xaI,uaR + 1j*uaI,\
										ibR + 1j*ibI, xbR + 1j*xbI,ubR + 1j*ubI,\
										icR + 1j*icI, xcR + 1j*xcI,ucR + 1j*ucI,\
										Vdc,xDC,xQ,\
										xPLL,wte)
		
			self.update_Ppv(t)
			self.update_Zload1(t)		
		
			self.update_voltages()
			self.update_power()
			self.update_RMS()
		
			self.update_Qref(t)
			self.update_Vdc_ref(t)	
			self.update_iref()
		
			self.update_inverter_frequency(t)
		
			self.update_ridethrough_flags(t)
			self.disconnect_or_reconnect(t)
		
			#Phase a inverter output current
			diaR = (1/self.Lf)*(-self.Rf*self.ia.real - self.va.real + self.vta.real) + (self.winv/self.wbase)*self.ia.imag 
			diaI = (1/self.Lf)*(-self.Rf*self.ia.imag - self.va.imag + self.vta.imag) - (self.winv/self.wbase)*self.ia.real	
			 
			#Current controller dynamics
			if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit:
				if np.sign(self.Ki_GCC*self.ua.real) == np.sign(self.xa.real):
					dxaR = 0.0
				else:
					dxaR = self.Ki_GCC*self.ua.real
				if np.sign(self.Ki_GCC*self.ua.imag) == np.sign(self.xa.imag):
					dxaI = 0.0
				else:
					dxaI = self.Ki_GCC*self.ua.imag
					#six.print_(dxaR+1j*dxaI,np.sign(self.Ki_GCC*self.ua))
			else:
				dxaR = self.Ki_GCC*self.ua.real
				dxaI = self.Ki_GCC*self.ua.imag
			
			if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit:
				if np.sign( (self.wp)*(-self.ua.real +	self.ia_ref.real - self.ia.real)) == np.sign(self.ua.real):
					duaR = 0.0
				else:
					duaR = (self.wp)*(-self.ua.real +	self.ia_ref.real - self.ia.real)
				
				if np.sign((self.wp)*(-self.ua.imag +	self.ia_ref.imag - self.ia.imag)) == np.sign(self.ua.imag):
					duaI = 0.0
				else:
					duaI = (self.wp)*(-self.ua.imag +	self.ia_ref.imag - self.ia.imag)
				
			else:
				duaR = (self.wp)*(-self.ua.real +	self.ia_ref.real - self.ia.real)
				duaI = (self.wp)*(-self.ua.imag +	self.ia_ref.imag - self.ia.imag)
		
			#Phase b inverter output current
			dibR = (1/self.Lf)*(-self.Rf*self.ib.real - self.vb.real + self.vtb.real) + (self.winv/self.wbase)*self.ib.imag 
			dibI = (1/self.Lf)*(-self.Rf*self.ib.imag - self.vb.imag + self.vtb.imag) - (self.winv/self.wbase)*self.ib.real	
		
			#Current controller dynamics - Phase b
			#dxbR = self.Ki_GCC*self.ub.real
			#dxbI = self.Ki_GCC*self.ub.imag
			#dubR = (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)
			#dubI = (self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)
		
			if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
				if np.sign(self.Ki_GCC*self.ub.real) == np.sign(self.xb.real):
					dxbR = 0.0
				else:
					dxbR = self.Ki_GCC*self.ub.real
				if np.sign(self.Ki_GCC*self.ub.imag) == np.sign(self.xb.imag):
					dxbI = 0.0
				else:
					dxbI = self.Ki_GCC*self.ub.imag
				
			else: 
				dxbR = self.Ki_GCC*self.ub.real
				dxbI = self.Ki_GCC*self.ub.imag
		
			if abs(self.Kp_GCC*self.ub + self.xb)>self.m_limit:
				if np.sign( (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)) == np.sign(self.ub.real):
					dubR = 0.0
				else:
					dubR = (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)
				
				if np.sign((self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)) == np.sign(self.ub.imag):
					dubI = 0.0
				else:
					dubI = (self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)
			else:
				dubR = (self.wp)*(-self.ub.real +	self.ib_ref.real - self.ib.real)
				dubI = (self.wp)*(-self.ub.imag +	self.ib_ref.imag - self.ib.imag)
		
			#Phase c inverter output current
			dicR = (1/self.Lf)*(-self.Rf*self.ic.real - self.vc.real + self.vtc.real) + (self.winv/self.wbase)*self.ic.imag 
			dicI = (1/self.Lf)*(-self.Rf*self.ic.imag - self.vc.imag + self.vtc.imag) - (self.winv/self.wbase)*self.ic.real 
		
			#Current controller dynamics - Phase c
			#dxcR = self.Ki_GCC*self.uc.real
			#dxcI = self.Ki_GCC*self.uc.imag
			#ducR = (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)
			#ducI = (self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)
		
			if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
				if np.sign(self.Ki_GCC*self.uc.real) == np.sign(self.xc.real):
					dxcR = 0.0
				else:
					dxcR = self.Ki_GCC*self.uc.real
				
				if np.sign(self.Ki_GCC*self.uc.imag) == np.sign(self.xc.imag):
					dxcI = 0.0
				else:
					dxcI = self.Ki_GCC*self.uc.imag
				
			else: 
				dxcR = self.Ki_GCC*self.uc.real
				dxcI = self.Ki_GCC*self.uc.imag
		
			if abs(self.Kp_GCC*self.uc + self.xc)>self.m_limit:
				if np.sign( (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)) == np.sign(self.uc.real):
					ducR = 0.0
				else:
					ducR = (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)
			
				if np.sign((self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)) == np.sign(self.uc.imag):
					ducI = 0.0
				else:
					ducI = (self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)
			else:
				ducR = (self.wp)*(-self.uc.real +	self.ic_ref.real - self.ic.real)
				ducI = (self.wp)*(-self.uc.imag +	self.ic_ref.imag - self.ic.imag)
		
			#DC link voltage dynamics
			dVdc = (self.Ppv - self.S.real)/(self.Vdc*self.C)
		
			if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ	- self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
				if np.sign(self.Ki_DC*(self.Vdc_ref - self.Vdc)) == np.sign(self.xDC):
					dxDC = 0.0
				else:
					dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
			else:
				dxDC = self.Ki_DC*(self.Vdc_ref - self.Vdc)
			
			# Reactive power controller dynamics
			if abs(self.xDC + self.Kp_DC*(self.Vdc_ref - self.Vdc) + 1j*(self.xQ	- self.Kp_Q*(self.Q_ref - self.S_PCC.imag)))>self.iref_limit:
			
				if np.sign(-self.Ki_Q*(self.Q_ref - self.S_PCC.imag)) == np.sign(self.xQ):
					dxQ = 0.0
				else:
					dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
			else:
				dxQ = -self.Ki_Q*(self.Q_ref - self.S_PCC.imag)
		
			#SRF-PLL dynamics
			dxPLL = self.Ki_PLL*(self.vd)
		
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
			iaR, iaI, xaR, xaI, uaR, uaI,\
			ibR, ibI, xbR, xbI, ubR, ubI,\
			icR, icI, xcR, xcI, ucR, ucI,\
			Vdc, xDC, xQ, xPLL, wte = y	 # unpack current values of y
		
			self.update_inverter_states(iaR + 1j*iaI, xaR + 1j*xaI,uaR + 1j*uaI,\
										ibR + 1j*ibI, xbR + 1j*xbI,ubR + 1j*ubI,\
										icR + 1j*icI, xcR + 1j*xcI,ucR + 1j*ucI,\
										Vdc,xDC,xQ,\
										xPLL,wte)

			J = self.J
			varInd = self.varInd 
			self.update_Ppv(t)
		
			self.update_voltages()
			self.update_power()
			self.update_RMS()
		
			self.update_Qref(t)
			self.update_Vdc_ref(t)	
			self.update_iref()
		
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
			if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit*1e1:
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
		
			if abs(self.Kp_GCC*self.ua + self.xa)>self.m_limit*1e1:
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


