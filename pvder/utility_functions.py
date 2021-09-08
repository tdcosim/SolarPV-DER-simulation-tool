"""Commonly used calculations on electrical quantities."""

from __future__ import division
import numpy as np
import math
import cmath
import sys
import time
import six
import json
import logging
import scipy.io as sio
from pvder.logutil import LogUtil


def Urms_time_series(ua,ub,uc):
	"""Function to calculate rms value of phasor quantities."""
	try:
		assert len(ua) == len(ub) == len(uc),  " The number of phasor quantities should be equal"
		return np.sqrt((np.square(np.abs(ua))+np.square(np.abs(ub))+np.square(np.abs(uc)))/3.0)/math.sqrt(2)
	except:
		LogUtil.exception_handler()


def Uphrms_time_series(uph):
	"""Function to calculate rms value of single phasor quantity."""
	try:
		return np.abs(uph)/math.sqrt(2)
	except:
		LogUtil.exception_handler()


def Uabsolute_time_series(u):
	"""Function to calculate rms value of phasor quantities."""
	try:
		return np.abs(u)
	except:
		LogUtil.exception_handler()


def Uunbalance_calc(ua,ub,uc):
	"""Calculate voltage/current unbalance."""
	try:
		uavg = (ua + ub + uc)/3
		return (max(ua,ub,uc) - min(ua,ub,uc))/uavg
	except:
		LogUtil.exception_handler()


#@jit(nopython=True)
def Ppv_calc(Iph,Np,Ns,Vdc_actual,Tactual,Sbase):
	"""Function to calculate PV module power output."""
	try:
		Irs = 1.2e-7		 #Cell reverse saturation current
		q = 1.602e-19		#Charge of an electron
		k = 1.38e-23		#Boltzmann's constant
		A = 1.92	  #p-n junction ideality factor
		Ipv = (Np*Iph)-(Np*Irs*(math.exp((q*Vdc_actual)/(k*Tactual*A*Ns))-1))   #Faster  with Pure Python functions
		return max(0,(Ipv*Vdc_actual))/Sbase	
	except:
		LogUtil.exception_handler()


#@jit(nopython=True)
def S_calc(va,vb,vc,ia,ib,ic):
	"""Function to calculate apparent power."""
	try:
		return (1/2)*(va*ia.conjugate() + vb*ib.conjugate() + vc*ic.conjugate())*1.0
	except:
		LogUtil.exception_handler()


#Average duty cycle - Phase A
#@jit(nopython=True)
def m_calc(Kp_GCC,u,x):
	"""Duty cycle for a single phase."""
	try:
		return Kp_GCC*u + x #PI controller equation  
	except:
		LogUtil.exception_handler()


#@jit(nopython=True)
def Urms_calc(ua,ub,uc):
	"""Function to calculate rms value of scalar phasor quantities."""
	try:
		return math.sqrt((pow(abs(ua),2)+pow(abs(ub),2)+pow(abs(uc),2))/3.0)/math.sqrt(2)  #Pure python implementation is faster	  
	except:
		LogUtil.exception_handler()


def Urms_min_calc(ua,ub,uc):
	"""Function to calculate minimum of rms value of scalar phasor quantities."""
	try:
		return min(abs(ua),abs(ub),abs(uc))/math.sqrt(2)  #Pure python implementation is faster 
	except:
		LogUtil.exception_handler()


def Urms_calc_1phase(ua):
	"""Function to calculate rms value of scalar phasor quantities for single phase."""
	try:
		return abs(ua)/math.sqrt(2)  #Pure python implementation is faster	  
	except:
		LogUtil.exception_handler()


#@jit(nopython=True)	
def Ub_calc(Ua):
	"""Convert phase A quantity to Phase B."""
	try:
		return Ua*np.power(np.e,1j*(-(2/3)*np.pi))
		#return Ua*pow(math.e,1j*(-(2/3)*math.pi))  #Shift by -120 degrees
	#@jit(nopython=True)
	except:
		LogUtil.exception_handler()


def Uc_calc(Ua):
	"""Convert phase A quantity to Phase C."""
	try:
		return Ua*np.power(np.e,1j*((2/3)*np.pi))
		#return Ua*pow(math.e,1j*((2/3)*math.pi))  #Shift by -120 degrees
	except:
		LogUtil.exception_handler()


def relative_phase_calc(Uph1,Uph2,DEGREES=False):
	"""Calculate relative phase between phasors between 0 to 2pi or 0 to 360 degrees."""
	try:
		if DEGREES:
			del_phase = math.degrees(cmath.phase(Uph1)-cmath.phase(Uph2)) % 360
		else:
			del_phase = math.radians(math.degrees(cmath.phase(Uph1)-cmath.phase(Uph2)) % 360)
		return del_phase
	except:
		LogUtil.exception_handler()


#@jit(nopython=True)	
def phasor_to_time(upha = 1+1j*0.0,uphb = -0.5-1j*0.867,uphc = -0.5+1j*0.867,w=2.0*math.pi*60.0,t=0.0):
	"""Convert a,b,c quantities from phasor domain to time domain."""
	try:
		ra,pha = cmath.polar(upha)
		rb,phb = cmath.polar(uphb)
		rc,phc = cmath.polar(uphc)
		#ua = (ra*np.exp(1j*(w*t+pha-(math.pi/2)))).real
		#ub = (rb*np.exp(1j*(w*t+phb-(math.pi/2)))).real
		#uc = (rc*np.exp(1j*(w*t+phc-(math.pi/2)))).real
		ua = ra*pow(math.e,1j*(w*t+pha-(math.pi/2))).real
		ub = rb*pow(math.e,1j*(w*t+phb-(math.pi/2))).real
		uc = rc*pow(math.e,1j*(w*t+phc-(math.pi/2))).real
		return ua,ub,uc   
	except:
		LogUtil.exception_handler()


def phasor_to_time_1phase(uph,w,t):
	"""Convert a,b,c quantities (time series) from phasor domain to time domain."""
	try:
		r,ph = cmath.polar(uph)
		return r*pow(math.e,1j*(w*t+ph-(math.pi/2))).real
	except:
		LogUtil.exception_handler()


#@jit(nopython=True)
def abc_to_dq0(ua,ub,uc,wt=2*math.pi):
	"""Convert to d-q."""
	try:
		Us = (2/3)*(ua + ub*pow(math.e,1j*((2/3)*math.pi)) + uc*pow(math.e,1j*(-(2/3)*math.pi)))*pow(math.e,1j*(-wt))
		ud = Us.real
		uq = Us.imag
		u0 = (1/3)*(ua+ub+uc)
		return ud,uq,u0
	except:
		LogUtil.exception_handler()


def dq0_to_abc(ud,uq,u0,wt=2*math.pi):
	"""Convert to abc."""
	try:
		ua = ud*math.cos(wt) - uq*math.sin(wt) + u0
		ub = ud*math.cos(wt-(2/3)*math.pi) - uq*math.sin(wt-(2/3)*math.pi) + u0
		uc = ud*math.cos(wt+(2/3)*math.pi) - uq*math.sin(wt+(2/3)*math.pi) + u0
		return ua,ub,uc
	except:
		LogUtil.exception_handler()


def alpha_beta_to_d_q(ualpha,ubeta,wt):
	"""Convert alpha-beta to d-q."""
	try:
		Us = (ualpha + 1j*ubeta)*pow(math.e,-1j*(wt))
		#Us = (ualpha + 1j*ubeta)*pow(math.e,-1j*(wt-(math.pi/2)))
	
		ud = Us.real
		uq = Us.imag
	
		#print(ud,uq)
		#ud = ualpha*math.sin(wt) - ubeta*math.cos(wt)
		#uq = ualpha*math.cos(wt) + ubeta*math.sin(wt)
	
		#ud = ualpha*math.cos(wt) + ubeta*math.sin(wt)
		#uq = -ualpha*math.sin(wt) + ubeta*math.cos(wt)
	
		return ud,uq
	except:
		LogUtil.exception_handler()


def phasor_to_symmetrical(upha,uphb,uphc):
	"""Convert to zero sequence."""
	try:
		a = pow(math.e,1j*((2/3)*math.pi))
		aa = pow(math.e,1j*((4/3)*math.pi))
		u0 = (1/3)*(upha + uphb + uphc)
		u1 = (1/3)*(upha + a*uphb + (aa)*uphc)  #Positive sequence
		u2 = (1/3)*(upha + (aa)*uphb + a*uphc) #Negative sequence

		return u0,u1,u2
	except:
		LogUtil.exception_handler()


def phasor_to_zero_sequence(upha,uphb,uphc):
	"""Convert to zero sequence.""" 
	try:
		u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
		return u0,u0,u0
	except:
		LogUtil.exception_handler()


def phasor_to_positive_sequence(upha,uphb,uphc):
	"""Convert to positive sequence.""" 
	try:
		u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
		return u1,u1*pow(math.e,1j*((4/3)*math.pi)),u1*pow(math.e,1j*((2/3)*math.pi))
	except:
		LogUtil.exception_handler()


def phasor_to_negative_sequence(upha,uphb,uphc):
	"""Convert to negative sequence.""" 
	try:
		u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
		return u2,u2*pow(math.e,1j*((4/3)*math.pi)),u2*pow(math.e,1j*((2/3)*math.pi))
	except:
		LogUtil.exception_handler()


def Vinv_terminal_time_series(m_t,Vdc_t):
	"""Function to generate time series inverter terminal voltage."""
	try:
		assert len(m_t) == len(Vdc_t) != None
		return m_t*(Vdc_t/2)
	except:
		LogUtil.exception_handler()

def m_time_series(u_t,x_t,Kp_GCC):
	"""Function to generate time series inverter terminal voltage."""
	try:
		assert len(u_t) == len(x_t) != None		
		return Kp_GCC*u_t + x_t
	except:
		LogUtil.exception_handler()


def extract_matlab_file(file_name,series_label):
	"""Program to extract contents of .mat file having structure with time format."""
	try:
		matlab_file = sio.loadmat(file_name)
		content_list = []
		sim_time = matlab_file[series_label][0,0][0]
		n_labels = len(matlab_file[series_label][0,0][1][0])
		for i in range(n_labels):
			content_list.append(matlab_file[series_label][0,0][1][0,i][0])
		return sim_time,content_list   
	except:
		LogUtil.exception_handler()


def limit_complex(z,low_limit=-1.0,high_limit=1.0):
	"""Check if complex number is within limits."""
	try:
		assert  low_limit < high_limit,'low limit should be higher than high limit'
		r,phi = cmath.polar(z)
		r = min(max(r,low_limit),high_limit)
	
		#z_real = min(max(z.real,low_limit),high_limit)
		#z_imag = min(max(z.imag,low_limit),high_limit)
		#z_imag = z.imag
		#return z_real + z_imag*1j
		return  cmath.rect(r, phi)
	except:
		LogUtil.exception_handler()


def limit_complex_time_series(z,low_limit=-1.0,high_limit=1.0):
	"""Check if complex number is within limits."""
	try:
		assert  low_limit < high_limit,'low limit should be higher than high limit'
		#z_real = np.minimum(np.maximum(z.real,np.full(len(z),-1)),np.full(len(z),1))
		#z_imag = np.minimum(np.maximum(z.imag,np.full(len(z),-1)),np.full(len(z),1))
		r = np.absolute(z)
		phi = np.angle(z)
		r = np.minimum(np.maximum(r,np.full(len(r),low_limit)),np.full(len(r),high_limit))
	
		#z_imag = z.imag
		return r*np.cos(phi) + r*np.sin(phi)*1j
	except:
		LogUtil.exception_handler()


def print_to_terminal(text_string='Printing to terminal!'):
	"""Print to terminal."""
	try:
		sys.__stdout__.write(text_string+'\n')
		sys.__stdout__.flush()	  #Flush buffer to terminal
	except:
		LogUtil.exception_handler()


def print_dictionary_keys(dictionary,dictionary_name):
	"""Print dictionary."""
	try:
		print(dictionary_name,':',','.join(dictionary))
	except:
		LogUtil.exception_handler()


def read_json(file_name):
	"""Load json file and return dictionary."""
	try:
		if six.PY3:
			string_type = (str)
		elif six.PY2:
			string_type = (str,unicode)		   
		assert isinstance(file_name,string_type),'File name should be a string!'
		assert file_name[-4:]=='json','File should be a JSON file!'
		with open(file_name, "r") as config_file:
			#logger.debug('Reading configuration file:{}'.format(config_file.name))
			confDict = json.load(config_file)

		return confDict
	except:
		LogUtil.exception_handler()


def get_logger(logging_level):
	"""Get logger."""
	try:
		logger=logging.getLogger()  #Creating an object 
		logger.setLevel(eval('logging.'+logging_level)) #Setting the threshold of logger to DEBUG 

		return logger
	except:
		LogUtil.exception_handler()


