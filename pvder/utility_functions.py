"""Commonly used calculations on electrical quantities."""

from __future__ import division
import numpy as np
import math
import cmath
import sys
import time
import six
import scipy.io as sio
#from numba import jit

def Urms_time_series(ua,ub,uc):
    """Function to calculate rms value of phasor quantities."""
    
    assert len(ua) == len(ub) == len(uc),  " The number of phasor quantities should be equal"
    
    return np.sqrt((np.square(np.abs(ua))+np.square(np.abs(ub))+np.square(np.abs(uc)))/3.0)/math.sqrt(2)

def Uabsolute_time_series(u):
    """Function to calculate rms value of phasor quantities."""
    
    return np.abs(u)

#@jit(nopython=True)
def Ppv_calc(Iph,Np,Ns,Vdc_actual,Tactual,Sbase):
    """Function to calculate PV module power output."""
    
    Irs = 1.2e-7         #Cell reverse saturation current
    q = 1.602e-19        #Charge of an electron
    k = 1.38e-23        #Boltzmann's constant
    A = 1.92      #p-n junction ideality factor
    
    Ipv = (Np*Iph)-(Np*Irs*(math.exp((q*Vdc_actual)/(k*Tactual*A*Ns))-1))   #Faster  with Pure Python functions
    return max(0,(Ipv*Vdc_actual))/Sbase    

#@jit(nopython=True)
def S_calc(va,vb,vc,ia,ib,ic):
    """Function to calculate apparent power."""
    
    return (1/2)*(va*ia.conjugate() + vb*ib.conjugate() + vc*ic.conjugate())*1.0

#Average duty cycle - Phase A
#@jit(nopython=True)
def m_calc(Kp_GCC,u,x):
    """Duty cycle for a single phase."""
        
    return Kp_GCC*u + x #PI controller equation  

#@jit(nopython=True)
def Urms_calc(ua,ub,uc):
    """Function to calculate rms value of scalar phasor quantities."""
    
    return math.sqrt((pow(abs(ua),2)+pow(abs(ub),2)+pow(abs(uc),2))/3.0)/math.sqrt(2)  #Pure python implementation is faster      
#@jit(nopython=True)    
def Ub_calc(Ua):
    """Convert phase A quantity to Phase B."""
    
    return Ua*pow(math.e,1j*(-(2/3)*math.pi))  #Shift by -120 degrees
#@jit(nopython=True)    
def Uc_calc(Ua):
    """Convert phase A quantity to Phase C."""
    
    return Ua*pow(math.e,1j*((2/3)*math.pi))  #Shift by -120 degrees
#@jit(nopython=True)    
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

def phasor_to_time_1phase(uph,w,t):
    """Convert a,b,c quantities (time series) from phasor domain to time domain."""
    
    r,ph = cmath.polar(uph)
    
    return r*pow(math.e,1j*(w*t+ph-(math.pi/2))).real
#@jit(nopython=True)
def abc_to_dq0(ua,ub,uc,wt=2*math.pi):
    """Convert to d-q."""
    Us = (2/3)*(ua + ub*pow(math.e,1j*((2/3)*math.pi)) + uc*pow(math.e,1j*(-(2/3)*math.pi)))*pow(math.e,1j*(-wt))
        
    ud = Us.real
    uq = Us.imag
    u0 = (1/3)*(ua+ub+uc)
    return ud,uq,u0

def dq0_to_abc(ud,uq,u0,wt=2*math.pi):
    """Convert to abc."""
    
    ua = ud*math.cos(wt) - uq*math.sin(wt) + u0
    ub = ud*math.cos(wt-(2/3)*math.pi) - uq*math.sin(wt-(2/3)*math.pi) + u0
    uc = ud*math.cos(wt+(2/3)*math.pi) - uq*math.sin(wt+(2/3)*math.pi) + u0

    return ua,ub,uc

def alpha_beta_to_d_q(ualpha,ubeta,wt):
    """Convert alpha-beta to d-q."""
    
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

def phasor_to_symmetrical(upha,uphb,uphc):
    """Convert to zero sequence."""     
    a = pow(math.e,1j*((2/3)*math.pi))
    aa = pow(math.e,1j*((4/3)*math.pi))
    u0 = (1/3)*(upha + uphb + uphc)
    u1 = (1/3)*(upha + a*uphb + (aa)*uphc)  #Positive sequence
    u2 = (1/3)*(upha + (aa)*uphb + a*uphc) #Negative sequence
    
    return u0,u1,u2

def phasor_to_zero_sequence(upha,uphb,uphc):
    """Convert to zero sequence.""" 
    
    u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
    
    return u0,u0,u0	

def phasor_to_positive_sequence(upha,uphb,uphc):
    """Convert to positive sequence.""" 
    
    u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
    
    return u1,u1*pow(math.e,1j*((4/3)*math.pi)),u1*pow(math.e,1j*((2/3)*math.pi))

def phasor_to_negative_sequence(upha,uphb,uphc):
    """Convert to negative sequence.""" 
    
    u0,u1,u2 = phasor_to_symmetrical(upha,uphb,uphc)
    
    return u2,u2*pow(math.e,1j*((4/3)*math.pi)),u2*pow(math.e,1j*((2/3)*math.pi))
    
def Vinv_terminal_time_series(m_t,Vdc_t):
    """Function to generate time series inverter terminal voltage."""
    
    assert len(m_t) == len(Vdc_t) != None        
    
    return m_t*(Vdc_t/2)

def m_time_series(u_t,x_t,Kp_GCC):
    """Function to generate time series inverter terminal voltage."""
    
    assert len(u_t) == len(x_t) != None        
        
    return Kp_GCC*u_t + x_t

def extract_matlab_file(file_name,series_label):
    """Program to extract contents of .mat file having structure with time format."""
    
    matlab_file = sio.loadmat(file_name)
    content_list = []
    sim_time = matlab_file[series_label][0,0][0]
    n_labels = len(matlab_file[series_label][0,0][1][0])
    for i in range(n_labels):
        content_list.append(matlab_file[series_label][0,0][1][0,i][0])
    return sim_time,content_list   
 
def limit_complex(z,low_limit=-1.0,high_limit=1.0):
    """Check if complex number is within limits."""
    
    assert  low_limit < high_limit,'low limit should be higher than high limit'
    r,phi = cmath.polar(z)
    r = min(max(r,low_limit),high_limit)
    
    #z_real = min(max(z.real,low_limit),high_limit)
    #z_imag = min(max(z.imag,low_limit),high_limit)
    #z_imag = z.imag
    #return z_real + z_imag*1j
    return  cmath.rect(r, phi)
    
def limit_complex_time_series(z,low_limit=-1.0,high_limit=1.0):
    """Check if complex number is within limits."""
    
    assert  low_limit < high_limit,'low limit should be higher than high limit'
    #z_real = np.minimum(np.maximum(z.real,np.full(len(z),-1)),np.full(len(z),1))
    #z_imag = np.minimum(np.maximum(z.imag,np.full(len(z),-1)),np.full(len(z),1))
    
    r = np.absolute(z)
    phi = np.angle(z)
    r = np.minimum(np.maximum(r,np.full(len(r),low_limit)),np.full(len(r),high_limit))
    
    #z_imag = z.imag
    return r*np.cos(phi) + r*np.sin(phi)*1j

def print_to_terminal(text_string='Printing to terminal!'):
    """Print to terminal."""
    
    sys.__stdout__.write(text_string+'\n')
    sys.__stdout__.flush()      #Flush buffer to terminal
"""
def print_LVRT_events(simulation_time,voltage,timer_start=0.0,event_name='',print_inline = False,verbose = False):
        #Print logs for LVRT events.

        if event_name == 'LV1_start':
            text_string = '{}:{time_stamp:.4f}:LV1 zone entered at {timer_start:.4f}s for {voltage:.3f} V'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage)

        elif event_name == 'LV2_start':
            text_string = '{}:{time_stamp:.4f}:LV2 zone entered at {timer_start:.4f}s for {voltage:.3f} V'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage)

        elif event_name == 'LV1_reset':
            text_string = '{}:{time_stamp:.4f}:LV1 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LV1 zone for {voltage:.3f} V'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage)

        elif event_name == 'LV2_reset':
            text_string = '{}:{time_stamp:.4f}:LV2 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LV2 zone for {voltage:.3f} V'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage)

        elif event_name == 'LV1_zone' and verbose == True:
            text_string = '{}:{time_stamp:.4f}:LV1 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'LV2_zone' and verbose == True:
            text_string = '{}:{time_stamp:.4f}:LV2 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'inverter_trip_LV1':
            text_string = '{}:{time_stamp:.4f}:LV1 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {voltage:.3f} V - Inverter will be tripped'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage)
            six.print_(text_string)

        elif event_name == 'inverter_trip_LV2':
            text_string = '{}:{time_stamp:.4f}:LV2 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {voltage:.3f} V - Inverter will be tripped'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage)    
            six.print_(text_string)

        elif event_name == 'reconnect_start':
            text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s for {voltage:.3f} V'\
                            .format(self.name,time_stamp=simulation_time,timer_start=timer_start,voltage=voltage)

        elif event_name == 'reconnect_reset':
            text_string = '{}:{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {voltage:.3f} V'\
                           .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage)

        elif event_name == 'reconnect_zone' and verbose == True:
            text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
                           .format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

        elif event_name == 'inverter_reconnection':
            text_string = '{}:{time_stamp:.4f}:Inverter reconnecting after LV trip at {time_stamp:.4f}s after {time_elasped:.4f}s for {voltage:.3f} V'\
                            .format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage)
            six.print_(text_string)

        elif event_name == 'inverter_tripped' and verbose == True: 
            text_string = '{}:{time_stamp:.4f}:Inverter in tripped condition for {voltage:.3f} V'.format(self.name,time_stamp=simulation_time,voltage=voltage)

        else:
            text_string =''

        if text_string != '':
            if print_inline == True:  #Print in notebook window
                six.print_(text_string)

            else: #Print in console window
                print_to_terminal(text_string)
        else:
            pass
"""
def print_LFRT_events(simulation_time,frequency,timer_start=0.0,event_name='',print_inline = False,verbose = False):
    """Print LFRT events."""    
    
    if event_name == 'LF1_start':
        text_string = '{time_stamp:.4f}:LF1 zone entered at {timer_start:.4f}s for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)
    
    elif event_name == 'LF2_start':
        text_string = '{time_stamp:.4f}:LF2 zone entered at {timer_start:.4f}s for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)
    
    elif event_name == 'LF3_start':
        text_string = '{time_stamp:.4f}:LF3 zone entered at {timer_start:.4f}s for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)
       
    elif event_name == 'LF1_reset':
        text_string = '{time_stamp:.4f}:LF1 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LF1 zone for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
    
    elif event_name == 'LF2_reset':
        text_string = '{time_stamp:.4f}:LF2 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LF2 zone for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
    
    elif event_name == 'LF3_reset':
        text_string = '{time_stamp:.4f}:LF3 flag reset at {time_stamp:.4f}s after {time_elasped:.4f} s in LF3 zone for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
        
    elif event_name == 'LF1_zone' and verbose == True:
        text_string = '{time_stamp:.4f}:LF1 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                        .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
    
    elif event_name == 'LF2_zone' and verbose == True:
        text_string = '{time_stamp:.4f}:LF2 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                        .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
    
    elif event_name == 'LF3_zone' and verbose == True:
        text_string = '{time_stamp:.4f}:LF3 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
                        .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
    
    
    elif event_name == 'inverter_trip_LF1':
        text_string = '{time_stamp:.4f}:LF1 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
                        .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
        six.print_(text_string)
    
    elif event_name == 'inverter_trip_LF2':
        text_string = '{time_stamp:.4f}:LF2 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
                        .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)    
        six.print_(text_string)
    
    elif event_name == 'inverter_trip_LF3':
        text_string = '{time_stamp:.4f}:LF3 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
                        .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)    
        six.print_(text_string)
    
    
    elif event_name == 'reconnect_start':
        text_string = '{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)
            
    elif event_name == 'reconnect_reset':
        text_string = '{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {frequency:.3f} Hz'\
                       .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
    
    elif event_name == 'reconnect_zone' and verbose == True:
        text_string = '{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
                       .format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
           
    elif event_name == 'inverter_reconnection':
        text_string = '{time_stamp:.4f}:Inverter reconnecting after LF trip at {time_stamp:.4f}s after {time_elasped:.4f}s for {frequency:.3f} Hz'\
                        .format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
        six.print_(text_string)
    
    elif event_name == 'inverter_tripped' and verbose == True: 
        text_string = '{time_stamp:.4f}:Inverter in tripped condition for {frequency:.3f} Hz'.format(time_stamp=simulation_time,frequency=frequency)
    
    else:
        text_string =''
    
    if text_string != '':
        if print_inline == True:  #Print in notebook window
            logging.info(text_string)
        
        else: #Print in console window
            sys.__stdout__.write(text_string+'\n')
            sys.__stdout__.flush()      #Flush buffer to terminal
    else:
        pass
