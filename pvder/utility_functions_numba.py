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
import os
import glob
import scipy.io as sio

from pvder.logutil import LogUtil
import numba
from numba import float32, float64, complex64,complex128

debug_flag =  False
cache_flag =  True

def check_numba_cache(cache_folder):
    """Check Numba cache folder"""
    
    files = glob.glob(os.path.join(cache_folder,'*.nbi')) + glob.glob(os.path.join(cache_folder,'*.nbc'))
    print("Numba cache folder contains {} *.nbi,*.nbc files".format(len(files)))
    
def clear_numba_cache(cache_folder):
    """Clear Numba cache folder"""
    
    files = glob.glob(os.path.join(cache_folder,'*.nbi')) + glob.glob(os.path.join(cache_folder,'*.nbc'))
    print("Clearing {} *.nbi,*.nbc files from Cache".format(len(files)))
    for f in files:
          os.remove(f)
    files = glob.glob(os.path.join(cache_folder,'*.nbi')) + glob.glob(os.path.join(cache_folder,'*.nbc'))
    
    if len(files) >0:
        print("Following files could not be cleared:{}".format(files))

@numba.njit(cache=cache_flag)
def Uunbalance_calc(ua,ub,uc):
	"""Calculate voltage/current unbalance."""
	uavg = (ua + ub + uc)/3
	return (max(ua,ub,uc) - min(ua,ub,uc))/uavg

@numba.njit(cache=cache_flag)
def Ppv_calc(Iph,Np,Ns,Vdc_actual,Tactual,Sbase):
	"""Function to calculate PV module power output."""
	Irs = 1.2e-7		 #Cell reverse saturation current
	q = 1.602e-19		#Charge of an electron
	k = 1.38e-23		#Boltzmann's constant
	A = 1.92	  #p-n junction ideality factor
	Ipv = (Np*Iph)-(Np*Irs*(math.exp((q*Vdc_actual)/(k*Tactual*A*Ns))-1))   #Faster  with Pure Python functions
	return max(0,(Ipv*Vdc_actual))/Sbase	

@numba.njit(cache=cache_flag)
def we_calc(xPLL,vd,Kp_PLL):
	"""Calculate inverter frequency from PLL."""
	return  (Kp_PLL*(vd) + xPLL + 2*math.pi*60.0)#/wbase

@numba.njit(cache=cache_flag)
def S_calc(va,vb,vc,ia,ib,ic):
	"""Function to calculate apparent power."""
	return (1/2)*(va*ia.conjugate() + vb*ib.conjugate() + vc*ic.conjugate())*1.0

@numba.njit(cache=cache_flag)
def S_load_calc(va,vb,vc,Zload):
	"""Power absorbed by load at PCC LV side."""
	return (1/2)*(va*(-(va/Zload)).conjugate() + vb*(-(vb/Zload)).conjugate() + vc*(-(vc/Zload)).conjugate())

@numba.njit(cache=cache_flag)
def S_G_calc(va,vb,vc,vag,vbg,vcg,ia,ib,ic,Zload,a):
	"""Power absorbed/produced by grid voltage source."""
	return (1/2)*((-(ia-(va/Zload))/a).conjugate()*vag+(-(ib-(vb/Zload))/a).conjugate()*vbg+(-(ic-(vc/Zload))/a).conjugate()*vcg)
	
#Average duty cycle - Phase A
@numba.njit(cache=cache_flag) #complex128(float32, complex128 ,complex128),
def m_calc(Kp_GCC,u,x):
	"""Duty cycle for a single phase."""
	return Kp_GCC*u + x #PI controller equation  

@numba.njit(cache=cache_flag) #float32(complex128,complex128,complex128),
def ia_ref_Vdc_Q_control_calc(xDC,xQ,Vdc,S,Vdc_ref,Q_ref,Kp_DC,Kp_Q):
	"""Phase A current reference"""
	return xDC + Kp_DC*(Vdc_ref - Vdc) + 1j*(xQ  - Kp_Q*(Q_ref - S.imag)) #PI controller equation

@numba.njit(cache=cache_flag) #float32(complex128,complex128,complex128),
def i_load_calc(v,Zload):
	"""Current counsumed by load connected at PCC LV side -  Phase A/B/C."""
	return v/Zload

@numba.njit(cache=cache_flag) #complex128(complex128,float32),
def vt_calc(m,Vdc):
	"""Inverter terminal voltage -  Phase A/B/C"""
	return m*(Vdc/2)

@numba.njit(cache=cache_flag) #float32(complex128,complex128,complex128),
def Vabrms_calc(va,vb):
	"""Inverter terminal voltage - line to line	RMS"""
	return abs(va-vb)/math.sqrt(2)

@numba.njit(cache=cache_flag) #float32(complex128,complex128,complex128),
def Urms_calc(ua,ub,uc):
	"""Function to calculate rms value of scalar phasor quantities."""
	return math.sqrt((pow(abs(ua),2)+pow(abs(ub),2)+pow(abs(uc),2))/3.0)/math.sqrt(2)  #Pure python implementation is faster	  

@numba.njit(cache=cache_flag)
def Urms_min_calc(ua,ub,uc):
	"""Function to calculate minimum of rms value of scalar phasor quantities."""
	return min(abs(ua),abs(ub),abs(uc))/math.sqrt(2)  #Pure python implementation is faster 

@numba.njit(cache=cache_flag)
def Urms_calc_1phase(ua):
	"""Function to calculate rms value of scalar phasor quantities for single phase."""
	return abs(ua)/math.sqrt(2)  #Pure python implementation is faster	  

@numba.njit(cache=cache_flag)
def Ub_calc(Ua):
	"""Convert phase A quantity to Phase B."""
	return Ua*np.power(np.e,1j*(-(2/3)*np.pi))
	#return Ua*pow(math.e,1j*(-(2/3)*math.pi))  #Shift by -120 degrees

@numba.njit(cache=cache_flag)
def Uc_calc(Ua):
	"""Convert phase A quantity to Phase C."""
	return Ua*np.power(np.e,1j*((2/3)*np.pi))
	#return Ua*pow(math.e,1j*((2/3)*math.pi))  #Shift by -120 degrees

@numba.njit(cache=cache_flag)
def relative_phase_calc(Uph1,Uph2,DEGREES=False):
	"""Calculate relative phase between phasors between 0 to 2pi or 0 to 360 degrees."""
	if DEGREES:
		del_phase = math.degrees(cmath.phase(Uph1)-cmath.phase(Uph2)) % 360
	else:
		del_phase = math.radians(math.degrees(cmath.phase(Uph1)-cmath.phase(Uph2)) % 360)
	return del_phase


@numba.njit(cache=cache_flag)
def phasor_to_time(upha = 1+1j*0.0,uphb = -0.5-1j*0.867,uphc = -0.5+1j*0.867,w=2.0*math.pi*60.0,t=0.0):
	"""Convert a,b,c quantities from phasor domain to time domain."""
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

@numba.njit(cache=cache_flag)
def phasor_to_time_1phase(uph,w,t):
	"""Convert a,b,c quantities (time series) from phasor domain to time domain."""
	r,ph = cmath.polar(uph)
	return r*pow(math.e,1j*(w*t+ph-(math.pi/2))).real

@numba.njit(cache=cache_flag)
def abc_to_dq0(ua,ub,uc,wt=2*math.pi):
	"""Convert to d-q."""
	Us = (2/3)*(ua + ub*pow(math.e,1j*((2/3)*math.pi)) + uc*pow(math.e,1j*(-(2/3)*math.pi)))*pow(math.e,1j*(-wt))
	ud = Us.real
	uq = Us.imag
	u0 = (1/3)*(ua+ub+uc)
	return ud,uq,u0

@numba.njit(cache=cache_flag)
def dq0_to_abc(ud,uq,u0,wt=2*math.pi):
	"""Convert to abc."""
	ua = ud*math.cos(wt) - uq*math.sin(wt) + u0
	ub = ud*math.cos(wt-(2/3)*math.pi) - uq*math.sin(wt-(2/3)*math.pi) + u0
	uc = ud*math.cos(wt+(2/3)*math.pi) - uq*math.sin(wt+(2/3)*math.pi) + u0
	return ua,ub,uc

@numba.njit(cache=cache_flag)
def alpha_beta_to_d_q(ualpha,ubeta,wt):
	"""Convert alpha-beta to d-q."""
	Us = (ualpha + 1j*ubeta)*pow(math.e,-1j*(wt))
	#Us = (ualpha + 1j*ubeta)*pow(math.e,-1j*(wt-(math.pi/2)))
	
	ud = Us.real
	uq = Us.imag
	
	#ud = ualpha*math.sin(wt) - ubeta*math.cos(wt)
	#uq = ualpha*math.cos(wt) + ubeta*math.sin(wt)
	
	#ud = ualpha*math.cos(wt) + ubeta*math.sin(wt)
	#uq = -ualpha*math.sin(wt) + ubeta*math.cos(wt)
	
	return ud,uq

@numba.njit(cache=cache_flag)
def phasor_to_symmetrical(upha,uphb,uphc):
	"""Convert to zero sequence."""

	a = pow(math.e,1j*((2/3)*math.pi))
	aa = pow(math.e,1j*((4/3)*math.pi))
	u0 = (1/3)*(upha + uphb + uphc)
	u1 = (1/3)*(upha + a*uphb + (aa)*uphc)  #Positive sequence
	u2 = (1/3)*(upha + (aa)*uphb + a*uphc) #Negative sequence

	return u0,u1,u2

@numba.njit(cache=cache_flag)
def phasor_to_zero_sequence(upha,uphb,uphc):
	"""Convert to zero sequence.""" 
	u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
	return u0,u0,u0

@numba.njit(cache=cache_flag)
def phasor_to_positive_sequence(upha,uphb,uphc):
	"""Convert to positive sequence.""" 
	u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
	return u1,u1*pow(math.e,1j*((4/3)*math.pi)),u1*pow(math.e,1j*((2/3)*math.pi))

@numba.njit(cache=cache_flag)
def phasor_to_negative_sequence(upha,uphb,uphc):
	"""Convert to negative sequence.""" 
	u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
	return u2,u2*pow(math.e,1j*((4/3)*math.pi)),u2*pow(math.e,1j*((2/3)*math.pi))

@numba.njit(cache=cache_flag)
def limit_complex(z,low_limit=-1.0,high_limit=1.0):
	"""Check if complex number is within limits."""
	#assert  low_limit < high_limit,'low limit should be higher than high limit'
	r,phi = cmath.polar(z)
	r = min(max(r,low_limit),high_limit)
	
	#z_real = min(max(z.real,low_limit),high_limit)
	#z_imag = min(max(z.imag,low_limit),high_limit)
	#z_imag = z.imag
	#return z_real + z_imag*1j
	return  cmath.rect(r, phi)


