"""Code for features inside PV-DER model instances."""

from __future__ import division
import six
import copy

import math

from pvder import utility_functions
from pvder import defaults,templates,specifications
from pvder.logutil import LogUtil


class PVDER_SmartFeatures():
	"""Class for describing smart inverter inverter features of PV-DER."""
	#Flags
	VOLT_VAR_ENABLE = False
	VOLT_VAR_FLAG = False
	
	VOLT_WATT_ENABLE = False
	VOLT_WATT_FLAG = False	
	
	V_threshold_high_limit = defaults.Vthreshold_high_limit #percentage
	t_threshold_high_limit = defaults.tthreshold_high_limit #seconds
	t_disconnect_low_limit = defaults.tdisconnect_low_limit #seconds #Minimum time to initiate DER disconnection (output cessation)
	t_reconnect_low_limit = defaults.treconnect_low_limit #seconds #Minimum time to initiate DER reconnection (output restoration)
	
	f_ref = 60.0
	DER_CONNECTED = True
	DER_MOMENTARY_CESSATION = False
	DER_TRIP	 = False		
	
	def initialize_Volt_VAR(self):
		"""Initialize the Volt-VAR controller settings."""
		try:
			self.Volt_VAR_dict = {'1':{'V':0.92,
										'Q':0.44},
								  '2':{'V':0.98,
										'Q':0.0},
								  '3':{'V':1.02,
										'Q':0.0},
								  '4':{'V':1.08,
										'Q':-0.44}
								  }  #From IEEE 1547-2018 Catergy B (Table 8 - page 39)		

			self.Volt_VAR_dict['m_inject'] = (self.Volt_VAR_dict['2']['Q'] - self.Volt_VAR_dict['1']['Q'])/(self.Volt_VAR_dict['2']['V'] - self.Volt_VAR_dict['1']['V'])
		
			self.Volt_VAR_dict['m_absorb'] = (self.Volt_VAR_dict['4']['Q'] - self.Volt_VAR_dict['3']['Q'])/(self.Volt_VAR_dict['4']['V'] - self.Volt_VAR_dict['3']['V'])
		
			self.Volt_VAR_dict['c_inject'] = self.Volt_VAR_dict['2']['Q'] - self.Volt_VAR_dict['2']['V']*self.Volt_VAR_dict['m_inject']
			self.Volt_VAR_dict['c_absorb'] = self.Volt_VAR_dict['4']['Q'] - self.Volt_VAR_dict['4']['V']*self.Volt_VAR_dict['m_absorb']
		
			self.VOLT_VAR_ACTIVE = False
		
			self.Qlimit = self.Qlimit_calc()
			self.Volt_VAR_logic_t = 0.0
		except:
			LogUtil.exception_handler()


	def Volt_VAR_logic(self,t):
		""" Volt-VAR."""
		try:
			#Select RMS voltage source
			Vrms_measured = self.Vrms
			self.Qlimit = self.Qlimit_calc()
		
			deadband_V = 0.01
		
			if t>self.Volt_VAR_logic_t: #Update volt-var logic only if solver time is greater that current time
				self.Volt_VAR_logic_t = t  
			
				if self.VOLT_VAR_ACTIVE:
				
					if self.Volt_VAR_inject_range(Vrms_measured,deadband_V):
						Qref = self.Q_Volt_VAR_inject(Vrms_measured)					
					elif self.Volt_VAR_absorb_range(Vrms_measured,deadband_V):
						Qref = self.Q_Volt_VAR_absorb(Vrms_measured)	   
					else:
						utility_functions.print_to_terminal("Volt-VAR control with reference voltage {:.3f} is deactivated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,Vrms_measured*self.Vbase))
						Qref = 0.0
						self.VOLT_VAR_ACTIVE = False					
				else:
					if self.Volt_VAR_inject_range(Vrms_measured,deadband_V):
						Qref = self.Q_Volt_VAR_inject(Vrms_measured)					
						utility_functions.print_to_terminal("Volt-VAR control (injection) with reference voltage {:.3f} is activated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,Vrms_measured*self.Vbase))
						self.VOLT_VAR_ACTIVE = True
					elif self.Volt_VAR_absorb_range(Vrms_measured,deadband_V):
						Qref = self.Q_Volt_VAR_absorb(Vrms_measured)
						utility_functions.print_to_terminal("Volt-VAR control (absorbption) with reference voltage {:.3f} is activated at {:.3f} s for {:.3f} V".format(self.Vrms_ref*self.Vbase,t,Vrms_measured*self.Vbase))			
						self.VOLT_VAR_ACTIVE = True

					else: #Don't supply/absorb reactive power outside operating range
						Qref = 0.0
			else:
				self.Volt_VAR_logic_t = self.Volt_VAR_logic_t
				Qref = self.Q_ref			
					
			return Qref
		except:
			LogUtil.exception_handler()


	def Qlimit_calc(self):
		"""Calculate maximum Q reference."""		
		try:
			Qmax = math.sqrt(math.pow(self.Sinverter_nominal,2)-math.pow(max(-self.Sinverter_nominal,min(self.S.real,self.Sinverter_nominal)),2))	#Qmax = sqrt(S^2 - P^2)
			Qlimit = Qmax - self.n_phases*self.Xf*(abs(self.ia)/math.sqrt(2))**2  #Qlimit = Qmax - Qfilter
			Qlimit = max(0.0,Qlimit)

			return Qlimit
		except:
			LogUtil.exception_handler()


	def Volt_VAR_inject_range(self,Vrms_measured,del_V=0.0):
		"""Check if voltage in volt-VAR operating range."""
		try:
			return Vrms_measured >= (self.Volt_VAR_dict['1']['V'] - del_V)*self.Vrms_ref and \
			Vrms_measured <= (self.Volt_VAR_dict['2']['V'] + del_V)*self.Vrms_ref
		except:
			LogUtil.exception_handler()


	def Volt_VAR_absorb_range(self,Vrms_measured,del_V=0.0):
		"""Check if voltage in volt-VAR operating range."""
		try:
			return Vrms_measured >= (self.Volt_VAR_dict['3']['V'] - del_V)*self.Vrms_ref and \
			Vrms_measured <= (self.Volt_VAR_dict['4']['V'] + del_V)*self.Vrms_ref	   
		except:
			LogUtil.exception_handler()


	def Q_Volt_VAR_inject(self,Vrms_measured):
		"""Calculate reactive power for Volt-VAR control."""
		try:
			if self.Volt_VAR_inject_range(Vrms_measured,del_V=0.0):
				Vrms_measured = (Vrms_measured/self.Vrms_ref)
				Q = max(0.0,((Vrms_measured)*self.Volt_VAR_dict['m_inject'] + self.Volt_VAR_dict['c_inject'])*self.Sinverter_nominal)
				Q = min(self.Qlimit,Q)
			else:
				Q = self.Q_ref

			return Q
		except:
			LogUtil.exception_handler()


	def Q_Volt_VAR_absorb(self,Vrms_measured):
		"""Calculate reactive power for Volt-VAR control."""
		try:
			if self.Volt_VAR_absorb_range(Vrms_measured,del_V=0.0):
				Vrms_measured = (Vrms_measured/self.Vrms_ref)
				Q = min(0.0,((Vrms_measured)*self.Volt_VAR_dict['m_absorb'] + self.Volt_VAR_dict['c_absorb'])*self.Sinverter_nominal)
				Q = max(-self.Qlimit,Q)
			else:
				Q = self.Q_ref

			return Q
		except:
			LogUtil.exception_handler()


	def update_ridethrough_flags(self,t):
		"""Check VRT and FRT logic."""
		try:
			if t < self.RT_logic_t_lock:
				LogUtil.logger.debug('{}:{:.4f}:Time going backwards in VRT/FRT logic - previous time was {:.4f} s.'.format(self.name,t,self.RT_logic_t_lock))
		
			elif t > self.RT_logic_t_lock: #Go through logic only if time is greater than previous time
				self.RT_logic_t_lock = t
				#LVRT logic		
				if self.LVRT_ENABLE:
					self.LVRT(t)

				#HVRT logic
				if self.HVRT_ENABLE:
					self.HVRT(t)

				#LFRT trip logic
				if self.LFRT_ENABLE:
					self.FRT(t)

				self.DER_TRIP = self.LVRT_TRIP or self.HVRT_TRIP# or self.LFRT_TRIP
				self.DER_MOMENTARY_CESSATION = self.LVRT_MOMENTARY_CESSATION or self.HVRT_MOMENTARY_CESSATION
		except:
			LogUtil.exception_handler()


	def disconnect_or_reconnect(self,t):
		"""Check flags and ither disconnect or reconnect DER."""
		try:
			if t < self.connect_logic_t_lock:
				LogUtil.logger.debug('{}:{:.4f}:Time going backwards in disconnect/reconnect logic - previous time was {:.4f} s .'.format(self.name,t,self.connect_logic_t_lock))
			
			elif t > self.connect_logic_t_lock: #Go through logic only if time is greater than previous time
				self.connect_logic_t_lock = t
				if self.DER_CONNECTED:
					if self.DER_TRIP or self.DER_MOMENTARY_CESSATION:
						self.DER_disconnect_logic(t)	   #Update status of DER_CONNECTED
				else:
					if not self.DER_TRIP:
						self.DER_reconnect_logic(t)
		
			if not self.DER_CONNECTED: #Keep DER connected only if DER_CONNECTED is True
				self.DER_disconnect()
		except:
			LogUtil.exception_handler()


	def DER_disconnect_logic(self,t):
		"""Logic for disconnecting/deenergizing DER from grid."""
		try:
			assert self.DER_CONNECTED, 'Disconnection logic can only be used if DER is already connected.'
		
			if self.t_disconnect_start == 0.0: #Start disconnect timer
				text_string = '{}:{:.4f}:DER disconnect timer started.'.format(self.name,t)
			
				self.print_event(text_string,True) 
				self.t_disconnect_start = t
			elif t-self.t_disconnect_start < self.t_disconnect_delay: #Disconnect DER only after disconnect time delay has elapsed
				text_string = '{}:{:.4f}:DER is in disconnect timer zone.'.format(self.name,t)
				LogUtil.logger.debug(text_string)
			elif t-self.t_disconnect_start >= self.t_disconnect_delay: #Disconnect DER only after disconnect time delay has elapsed
				text_string = '{}:{:.4f}:DER will be disconnected (t_disconnect_timer:{:.4f} s,Vrms measured = {:.4f} V,LVRT Momentary cessation:{},LVRT Trip:{},HVRT Momentary cessation:{},HVRT Trip:{})'.format(self.name,t,
																					t-self.t_disconnect_start,self.get_Vrms_measured()*self.Vbase,self.LVRT_MOMENTARY_CESSATION,self.LVRT_TRIP,self.HVRT_MOMENTARY_CESSATION,self.HVRT_TRIP)
				
				self.print_event(text_string,True) 
				self.t_disconnect_start = 0.0
				self.DER_CONNECTED = False
		except:
			LogUtil.exception_handler()


	def DER_disconnect(self):
		"""Function to disconnect PV-DER from grid."""
		try:
			#LogUtil.logger.debug('{}:DER disconnected'.format(self.name))
			self.VOLT_VAR_ENABLE = False
			self.VOLT_WATT_ENABLE = False
		
			self.Q_ref = 0.0  #Set reactive power reference to zero
			self.Vdc_ref = self.Vdc  #Maintain DC link voltage
			self.Ppv = 0.0	 #Disconnect PV panel
		
			self.ia_ref = 0.0 + 1j*0.0
			self.ib_ref = 0.0 + 1j*0.0
			self.ic_ref = 0.0 + 1j*0.0
		except:
			LogUtil.exception_handler()


	def DER_reconnect_logic(self,t):
		"""Logic used to decide reconnection."""
		try:
			Vrms_measured = self.get_Vrms_measured()
			fgrid = self.we/(2.0*math.pi)  #Use grid frequency as estimated by PLL
			assert not self.DER_CONNECTED, 'Reconnection logic can only be used if DER is disconnected.'
		
			if self.DER_MOMENTARY_CESSATION:
				if self.t_reconnect_start > 0.0:
					self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='reconnect_reset')
					self.t_reconnect_start = 0.0
				elif self.t_reconnect_start == 0.0:
					self.print_reconnect_events(t, Vrms_measured,fgrid,event_name='DER_tripped')	#Inverter output remains zero
				
			elif t > 0.0:
				if self.t_reconnect_start > 0.0:
					if t-self.t_reconnect_start < self.t_reconnect_delay:
						self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='reconnect_zone')
			
					if t-self.t_reconnect_start > self.t_reconnect_delay:
						self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='DER_reconnection')
						self.DER_CONNECTED = True
						self.t_reconnect_start = 0.0
						if self.RESTORE_Vdc: #Whether to restore Vdc to pre-cessation setpoint
							self.Vdc_ref_ramp(t,self.set_Vdc_ref()*self.Vbase)
							LogUtil.logger.info('Ramping Vdc to pre-anomaly setpoint with Vdc reference list:{}'.format(self.Vdc_ref_list))
			
				elif self.t_reconnect_start == 0.0:
					self.t_reconnect_start = t
					self.print_reconnect_events(t, Vrms_measured,fgrid,self.t_reconnect_start,event_name='reconnect_start')
		except:
			LogUtil.exception_handler()


	def RT_initialize(self):
		"""Initialize VRT and FRT settings."""	  
		try:
			self.VRT_initialize()
			self.FRT_initialize()
		except:
			LogUtil.exception_handler()


	def VRT_initialize(self):
		"""Initialize LVRT and HVRT settings."""
		try:
			#Common ride through variables
			self.RT_logic_t_lock = 0.0
			self.connect_logic_t_lock = 0.0
			self.t_disconnect_start = 0.0
			self.t_reconnect_start = 0.0
		
			#LVRT flags	
			self.LVRT_ENABLE = True
			self.LVRT_TRIP = False
			self.LVRT_MOMENTARY_CESSATION = False
		
			#HVRT flags
			self.HVRT_ENABLE = True
			self.HVRT_TRIP = False
			self.HVRT_MOMENTARY_CESSATION = False
			 
			self.LVRT_dict = copy.deepcopy(self.RT_config['LVRT'])
			self.HVRT_dict = copy.deepcopy(self.RT_config['HVRT'])
			
			self.LVRT_levels = [int(level) for level in list(self.LVRT_dict)]
			self.HVRT_levels = [int(level) for level in list(self.HVRT_dict)]
			self.LVRT_levels.sort()
			self.HVRT_levels.sort()
			
			self.t_disconnect_delay = self.RT_config['VRT_delays']['output_cessation_delay']#(1/120.0)
			self.t_reconnect_delay = self.RT_config['VRT_delays']['output_restore_delay'] 
			self.RESTORE_Vdc = self.RT_config['VRT_delays']['restore_Vdc'] #Restore Vdc to nominal reference by using a ramp
			
			self.check_VRT_settings()
		except:
			LogUtil.exception_handler()
	def LVRT(self,t):
		"""Function to implement LVRT ridethrough and trip logic."""
		try:
			Vrms_measured = self.get_Vrms_measured()
			
			if t > self.t_stable: #Go through logic only after a short time delay
				for LVRT_key,LVRT_values in self.LVRT_dict.items():
					zone_name = 'LV'+str(LVRT_key)
					
					V_threshold = LVRT_values['V_threshold']*self.Vrms_ref #Convert % thresholds into p.u.
					t_threshold = LVRT_values['t_threshold']
					t_min_ridethrough = LVRT_values['t_min_ridethrough']
					LVRT_mode = LVRT_values['mode']
					
					if Vrms_measured < V_threshold and not LVRT_values['threshold_breach']: #Check if voltage below threshold
						if LVRT_values['t_start'] == 0.0: #Start timer if voltage goes below threshold
							LVRT_values['t_start']  = t
							self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,LVRT_mode,event_name='zone_entered')
						
						elif t-LVRT_values['t_start'] <= t_threshold: #Remain in LV zone and monitor
							if LVRT_mode == 'momentary_cessation' and not self.LVRT_MOMENTARY_CESSATION: #Go into momentary cessation
								if t-LVRT_values['t_start'] >= t_min_ridethrough:
									self.LVRT_MOMENTARY_CESSATION = True
									self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,LVRT_mode,event_name='momentary_cessation')
							else:
								self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],V_threshold,t_threshold,LVRT_mode,event_name='zone_continue')
						
						elif t-LVRT_values['t_start'] >= t_threshold: #Trip DER if timer exceeds threshold
							self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,LVRT_mode,event_name='trip')
							LVRT_values['threshold_breach'] = True
							self.LVRT_TRIP = True
							LVRT_values['t_start'] = 0.0
					
					elif  Vrms_measured > V_threshold: #Check if voltage above threshold
						if LVRT_values['t_start'] > 0.0: #Reset timer if voltage goes above  threshold
							self.print_VRT_events(t,Vrms_measured,zone_name,LVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,LVRT_mode,event_name='zone_reset')
							LVRT_values['t_start']  = 0.0 
							self.LVRT_MOMENTARY_CESSATION = False #Reset momentary cessation flags
						else: #Do nothing
							pass
		except:
			LogUtil.exception_handler()
	def HVRT(self,t):
		"""Function to implement HVRT ridethrough and trip logic."""
		try:
			Vrms_measured = self.get_Vrms_measured()
			if t > self.t_stable: #Go through logic only after a short time delay
				for HVRT_key,HVRT_values in self.HVRT_dict.items():
					zone_name = 'HV'+str(HVRT_key)
					V_threshold = HVRT_values['V_threshold']*self.Vrms_ref #Convert % thresholds into p.u.
					t_threshold = HVRT_values['t_threshold']
					t_min_ridethrough = HVRT_values['t_min_ridethrough']
					HVRT_mode = HVRT_values['mode']
					
					if Vrms_measured > V_threshold and not HVRT_values['threshold_breach']: #Check if voltage above threshold
						if HVRT_values['t_start'] == 0.0: #Start timer if voltage goes above threshold
							HVRT_values['t_start']  = t
							self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,HVRT_mode,event_name='zone_entered')
					
						elif t-HVRT_values['t_start'] <= t_threshold: #Remain in LV zone and monitor
							if HVRT_mode == 'momentary_cessation' and not self.HVRT_MOMENTARY_CESSATION: #Go into momentary cessation
								if t-HVRT_values['t_start'] >= t_min_ridethrough:
									self.HVRT_MOMENTARY_CESSATION = True
									self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,HVRT_mode,event_name='momentary_cessation')
							else:
								self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],V_threshold,t_threshold,HVRT_mode,event_name='zone_continue')
						
						elif t-HVRT_values['t_start'] >= t_threshold: #Trip DER if timer exceeds threshold
							self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,HVRT_mode,event_name='trip')
							HVRT_values['threshold_breach'] = True
							self.HVRT_TRIP = True
							HVRT_values['t_start'] = 0.0
					
					elif  Vrms_measured < V_threshold: #Check if voltage below threshold
						if HVRT_values['t_start'] > 0.0: #Reset timer if voltage goes below threshold
							self.print_VRT_events(t,Vrms_measured,zone_name,HVRT_values['t_start'],V_threshold,t_threshold,t_min_ridethrough,HVRT_mode,event_name='zone_reset')
							HVRT_values['t_start']  = 0.0 
							self.HVRT_MOMENTARY_CESSATION = False #Reset momentary cessation flags
						else: #Do nothing
							pass
		except:
			LogUtil.exception_handler()

	def get_Vrms_measured(self):
		"""Get Vrms measurement"""
		try:
			#Select RMS voltage source
			if specifications.RT_measurement_type == 'minimum':
				Vrms_measured  = self.Vrms_min  #Select minimum value of PCC-Voltages as per IEEE 1547-2018 reccomendations
			elif specifications.RT_measurement_type == 'average':
				Vrms_measured = self.Vrms   #Select PCC - LV side voltage  

			return Vrms_measured
		except:
			LogUtil.exception_handler()


	def show_RT_flags(self):
		"""Show all RT flags."""
		try:
			print('DER_MOMENTARY_CESSATION:{},LVRT_MOMENTARY_CESSATION:{},HVRT_MOMENTARY_CESSATION:{}'.format(self.DER_MOMENTARY_CESSATION,self.LVRT_MOMENTARY_CESSATION,self.HVRT_MOMENTARY_CESSATION))
			print('DER_CONNECTED:{},DER_TRIP:{}'.format(self.DER_CONNECTED,self.DER_TRIP))
		except:
			LogUtil.exception_handler()


	def check_VRT_settings(self):
		"""Sanity check for VRT settings."""
		try:
			for RT_type in ['LV','HV']:
				V_threshold_list = []
				t_threshold_list = []
			
				if RT_type == 'LV':
					VRT_dict = self.LVRT_dict
				elif RT_type == 'HV':
					VRT_dict = self.HVRT_dict
			
				for VRT_key,VRT_values in VRT_dict.items():
					zone_name = RT_type+str(VRT_key)
					V_threshold_list.append(VRT_values['V_threshold'])
					t_threshold_list.append(VRT_values['t_threshold'])
					voltage_error_text = '{} is not a valid voltage threshold for zone {}!'.format(VRT_values['V_threshold'],zone_name) 
					time_error_text = '{} is not a valid time threshold for zone {}!'.format(VRT_values['t_threshold'],zone_name) 
				
					if RT_type == 'LV':
						if VRT_values['V_threshold'] >= 1.0:
							raise ValueError(voltage_error_text) 
						if not six.PY2: #Since dictionaries don't preserve insertion order in Python 2
							if not V_threshold_list == sorted(V_threshold_list,reverse=True): #Check if voltage threshold values are in descending order
								raise ValueError('LVRT voltage limits {} are not equal to sorted limits {}!'.format(V_threshold_list,sorted(V_threshold_list,reverse=True)))

					elif RT_type == 'HV':
						if VRT_values['V_threshold'] <= 1.0:
							raise ValueError(voltage_error_text) 
						if not six.PY2: #Since dictionaries don't preserve insertion order in Python 2
							if not V_threshold_list == sorted(V_threshold_list,reverse=False): #Check if voltage values are in ascending order
								raise ValueError('HVRT voltage limits {} are not equal to sorted limits {}!'.format(V_threshold_list,sorted(V_threshold_list,reverse=False)))
					
					if VRT_values['V_threshold'] <= 0.0:
						raise ValueError(voltage_error_text) 
					if VRT_values['V_threshold'] >=self.V_threshold_high_limit:
						raise ValueError(voltage_error_text) 
				
					if VRT_values['t_threshold'] <= 0.0:
						raise ValueError(time_error_text) 
					if VRT_values['t_threshold'] >= self.t_threshold_high_limit:
						raise ValueError(time_error_text) 

					if VRT_values['mode'] not in {'mandatory_operation','momentary_cessation'}:
						raise ValueError('{} is not a valid mode for zone {}'.format(VRT_values['mode'],zone_name))	
			
				if self.t_disconnect_delay < self.t_disconnect_low_limit:
					raise ValueError('DER output cessation time delay {} s is infeasible!'.format(self.t_disconnect_delay))   
			
				if self.t_reconnect_delay < self.t_reconnect_low_limit:
					raise ValueError('DER output restore time delay {} s after momentary cessation is infeasible!'.format(self.t_reconnect_delay))   
		except:
			LogUtil.exception_handler()


	def print_VRT_events(self,simulation_time,voltage,zone_name,timer_start=0.0,V_threshold=None,t_threshold=None,t_min_ridethrough=None,LVRT_mode=None,event_name='',print_inline = True,verbose = False):
		"""Print logs for VRT events."""
		try:
			voltage_actual = voltage*self.Vbase#)/(self.Vrms_ref*self.Vbase)#175
			Vrms_ref_actual = self.Vrms_ref*self.Vbase
			V_threshold_actual = V_threshold*self.Vbase
			V_threshold_frac = V_threshold/self.Vrms_ref
			if event_name == 'zone_entered':
				text_string = '{}:{:.4f}:{} zone entered at {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V,V_thresh:{:.2f} V({:.2f}),t_thresh:{:.2f} s,t_min_ridethrough:{:.2f} s,mode:{})'.format(self.name,simulation_time,zone_name,timer_start,voltage_actual,Vrms_ref_actual,V_threshold_actual,V_threshold_frac,t_threshold,t_min_ridethrough,LVRT_mode)
		
			elif event_name == 'zone_reset':
				text_string = '{}:{:.4f}:{} flag reset at {:.4f} s after {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V,V_thresh:{:.2f} V({:.2f}),t_thresh:{:.2f} s,t_min_ridethrough:{:.2f} s,mode:{})'.format(self.name,simulation_time,zone_name,simulation_time,simulation_time-timer_start,voltage_actual,Vrms_ref_actual,V_threshold_actual,V_threshold_frac,t_threshold,t_min_ridethrough,LVRT_mode)
		
			elif event_name == 'zone_continue' and verbose:
				text_string = '{}:{:.4f}:{} zone entered at:{:.4f} s and continuing for {:.4f} s'\
								.format(self.name,simulation_time,zone_name,timer_start,simulation_time-timer_start)
			elif event_name == 'momentary_cessation':
				text_string = '{}:{:.4f}:{} zone - momentary cessation after {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V,V_thresh:{:.2f} V({:.2f}),t_thresh:{:.2f} s,t_min_ridethrough:{:.2f} s,mode:{})'.format(self.name,simulation_time,zone_name,simulation_time-timer_start,voltage_actual,Vrms_ref_actual,V_threshold_actual,V_threshold_frac,t_threshold,t_min_ridethrough,LVRT_mode)
			elif event_name == 'trip':
				text_string = '{}:{:.4f}:{} violation at {:.4f}s after {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V,V_thresh:{:.2f} V({:.2f}),t_thresh:{:.2f} s,t_min_ridethrough:{:.2f} s,mode:{}) - DER will be tripped'\
								.format(self.name,simulation_time,zone_name,simulation_time,simulation_time-timer_start,voltage_actual,Vrms_ref_actual,V_threshold_actual,V_threshold_frac,t_threshold,t_min_ridethrough,LVRT_mode)
			else:
				text_string =''

			self.print_event(text_string,print_inline)
		except:
			LogUtil.exception_handler()


	def print_reconnect_events(self,simulation_time,voltage,frequency,timer_start=0.0,event_name='',print_inline = True,verbose = False):
		"""Print logs for VRT events."""
		try:
			voltage_actual = voltage*self.Vbase#)/(self.Vrms_ref*self.Vbase)#175#175
			Vrms_ref_actual = self.Vrms_ref*self.Vbase
			if event_name == 'reconnect_start':
				text_string = '{}:{:.4f}:Reconnect timer started at {:.4f} s for {:.3f} V p.u. (Vref:{:.2f} V)'\
								.format(self.name,simulation_time,timer_start,voltage_actual,Vrms_ref_actual)

			elif event_name == 'reconnect_reset':
				text_string = '{}:{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
								.format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage_actual,Vref=Vrms_ref_actual)

			elif event_name == 'reconnect_zone' and verbose:
				text_string = '{}:{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
								.format(self.name,time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)

			elif event_name == 'DER_reconnection':
				text_string = '{}:{time_stamp:.4f}:DER reconnecting after momentary cessation at {time_stamp:.4f}s after {time_elasped:.4f}s for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'\
								.format(self.name,time_stamp=simulation_time,time_elasped=simulation_time-timer_start,voltage=voltage_actual,Vref=Vrms_ref_actual)
			

			elif event_name == 'DER_tripped' and verbose: 
				text_string = '{}:{time_stamp:.4f}:Inverter in tripped condition for {voltage:.3f} V p.u. (Vref:{Vref:.2f} V)'.format(self.name,time_stamp=simulation_time,voltage=voltage_actual,Vref=Vrms_ref_actual)

			else:
				text_string =''
		
			self.print_event(text_string,print_inline)
		except:
			LogUtil.exception_handler()


	def check_anomaly(self):
		"""Check if voltage anomaly was detected."""
		try:
			for LVRT_key,LVRT_values in self.LVRT_dict.items():
				if LVRT_values['t_start'] > 0.0: #Check if any timer started
					print('Voltage anomaly started at {:.4f} due to zone {} with threshold {:.2f} V p.u.'.format(LVRT_values['t_start'],LVRT_key,LVRT_values['V_threshold']))
		except:
			LogUtil.exception_handler()


	def update_RT_config_old(self,derConfig):
		"""Check whether the config file is good."""		
		try:
			for item in templates.RT_config_template.keys():
				if item in derConfig:
					self.RT_config[item] = derConfig[item]
				elif item in self.DER_config:
					self.RT_config[item] = self.DER_config[item]
					LogUtil.logger.debug('{}:{} updated with alue from config file {}.'.format(self.name,item,self.DER_config[item]))
				elif item in self.default_RT_config:
					LogUtil.logger.debug('{}:{} updated with default value {}.'.format(self.name,item,self.default_RT_config[item]))	
					self.RT_config[item] = self.default_RT_config[item]
				else:
					raise KeyError('{}:Ridethrough setting {} could not be found!'.format(self.name,item))
		except:
			LogUtil.exception_handler()


	def update_RT_config(self,config_dict):
		"""Check whether the config file is good."""
		try:			
			for RT in list(templates.VRT_config_template.keys()) +  list(templates.FRT_config_template.keys()):
				if 'config' in self.DER_config[RT]: #Check if an RT settings are provided
					self.RT_config[RT] = self.DER_config[RT]['config'] 
					LogUtil.logger.debug('{}:{} updated with {} from DER config file.'.format(self.name,RT,self.DER_config[RT]['config'] ))
				elif 'config_id' in self.DER_config[RT]: #Check if an RT config id is provided
					self.RT_config[RT] = config_dict[self.DER_config[RT]['config_id']]['config']
					LogUtil.logger.debug('{}:{} updated with {} from config id {}.'.format(self.name,RT,config_dict[self.DER_config[RT]['config_id']]['config'],self.DER_config[RT]['config_id']))
				
			
			for RT in ['LVRT','HVRT']:
				for RT_level,RT_config in self.RT_config[RT].items():
					
					assert isinstance(RT_config,dict), 'Ridethrough level must be specified as a dictionary!'
					if 't_start' not in RT_config.keys():
						RT_config.update({'t_start':0.0})
					if 'threshold_breach' not in RT_config.keys():
						RT_config.update({'threshold_breach':False})

			for config_parameter in templates.VRT_config_template['VRT_delays']['config']:
				if config_parameter not in self.RT_config['VRT_delays']:
					self.RT_config['VRT_delays'][config_parameter] = templates.VRT_config_template['VRT_delays']['config'][config_parameter]
					LogUtil.logger.debug('{}: Updated {} from template with value {}'.format(self.name,config_parameter,templates.VRT_config_template['VRT_delays']['config'][config_parameter]))
		except:
			LogUtil.exception_handler()


	def check_RT_config(self):
		"""Check whether the config file is good."""
		try:
			for RT in list(templates.VRT_config_template.keys()):
				if RT != 'VRT_delays':
					for RT_setting in templates.VRT_config_template[RT]['config']['1']:
						for level,setting in self.RT_config[RT].items():
							if RT_setting not in setting:
								raise ValueError('{} not found for level {} in {}!'.format(RT_setting,level,RT))
			
			for RT in list(templates.FRT_config_template.keys()):
				if RT != 'FRT_delays':
					for RT_setting in templates.FRT_config_template[RT]['config']['1']:
						for level,setting in self.RT_config[RT].items():
							if RT_setting not in setting:
								raise ValueError('{} not found for level {} in {}!'.format(RT_setting,level,RT))
		except:
			LogUtil.exception_handler()


	def show_RT_settings(self,settings_type='LVRT',PER_UNIT=True):
		"""Method to show LVRT settings."""
		try:
			if settings_type not in {'LVRT','HVRT','LFRT'}:#
				raise ValueError('Unknown quantity: ' + str(settings_type))
		
			if PER_UNIT:
				V_multiplier = (1/self.Vrms_ref)
			else:
				V_multiplier = self.Vbase
		
			print('\n______{} - {}_____'.format(self.name,settings_type))
		
			if settings_type ==  'LVRT':
				print('______Flags______')
				print('LVRT_ENABLE:{}\nLVRT_TRIP:{},LVRT_MOMENTARY_CESSATION:{} '.format(self.LVRT_ENABLE,self.LVRT_TRIP,self.LVRT_MOMENTARY_CESSATION))	
				print('______Thresholds______')
				print('Vrms_ref:{:.2f} V'.format(self.Vrms_ref*self.Vbase))
				for LVRT_key,LVRT_values in self.LVRT_dict.items():
					print('Zone:{},Vthreshold:{:.2f},tthreshold:{:.2f},tminridethrough:{:.3f},mode:{}'.format(LVRT_key,LVRT_values['V_threshold'],LVRT_values['t_threshold'],
                                                                                                                       LVRT_values['t_min_ridethrough'],LVRT_values['mode']))
			
			if settings_type ==  'HVRT':
				print('______Flags______')
				print('HVRT_ENABLE:{}\nHVRT_TRIP:{},HVRT_MOMENTARY_CESSATION:{}'.format(self.HVRT_ENABLE,self.HVRT_TRIP,self.HVRT_MOMENTARY_CESSATION))
				print('______Thresholds______')
				print('Vrms_ref:{:.2f} V'.format(self.Vrms_ref*self.Vbase))
				for HVRT_key,HVRT_values in self.HVRT_dict.items():
					print('Zone:{},Vthreshold:{:.2f},tthreshold:{:.2f},tminridethrough:{:.3f},mode:{}'.format(HVRT_key,HVRT_values['V_threshold'],HVRT_values['t_threshold'],
                                                                                                HVRT_values['t_min_ridethrough'],HVRT_values['mode']))
			
		
			if settings_type ==  'LFRT':
				print('f_ref:{:.2f} Hz\nF_LF1:{:.2f} Hz\nF_LF2:{:.2f} Hz'.format(self.f_ref,self.LFRT_dict['1']['F_LF'],self.LFRT_dict['2']['F_LF']))
				print('t_LF1:{:.2f} s\nt_LF2:{:.2f} s'.format(self.LFRT_dict['1']['t_LF_limit'],self.LFRT_dict['2']['t_LF_limit']))
				print('______Flags______')
				print('LFRT_ENABLE:{}\nLFRT_TRIP:{} '.format(self.LFRT_ENABLE,self.LFRT_TRIP))		
		
			if settings_type ==  'LFRT' or settings_type ==  'HFRT':
				print('FRT_INSTANTANEOUS_TRIP:{}'.format(self.FRT_INSTANTANEOUS_TRIP))		
			print('OUTPUT_CESSATION_DELAY:{},OUTPUT_RESTORE_DELAY:{}'.format(self.t_disconnect_delay,self.t_reconnect_delay))
		except:
			LogUtil.exception_handler()


	def FRT_initialize(self):
		"""Initialize LFRT and HFRT settings."""
		try:
			#LVRT flags	
			self.LFRT_ENABLE = False
			self.LFRT_TRIP = False		
		
			#HVRT flags
			self.HFRT_ENABLE = False
			self.HFRT_TRIP = False
		
			self.del_f =0.02

			self.LFRT_dict = self.RT_config['LFRT']
			self.HFRT_dict = self.RT_config['HFRT']
			self.FRT_INSTANTANEOUS_TRIP = self.RT_config['FRT_delays']['FRT_INSTANTANEOUS_TRIP'] #Disconnects PV-DER within one cycle for frequency anomaly
			self.check_LFRT_settings()	
			self.check_HFRT_settings()
		except:
			LogUtil.exception_handler()


	def check_LFRT_settings(self):
		"""Sanity check for LFRT settings."""		
		try:
			if not self.LFRT_dict['2']['F_LF'] > self.LFRT_dict['1']['F_LF']:
				raise ValueError('LFRT frequency settings - F_LF1:{:.2f},F_LF2:{:.2f} are infeasible!'.format(self.LFRT_dict['1']['F_LF'],self.LFRT_dict['2']['F_LF']))
		except:
			LogUtil.exception_handler()


	def check_HFRT_settings(self):
		"""Sanity check for HFRT settings."""
		try:
			if not self.HFRT_dict['2']['F_HF'] > self.HFRT_dict['1']['F_HF']:
				raise ValueError('HFRT frequency settings - F_HF1:{:.2f},F_HF2:{:.2f} are infeasible!'.format(self.HFRT_dict['1']['F_HF'],self.HFRT_dict['2']['F_HF']))
		except:
			LogUtil.exception_handler()


	def FRT(self,t):
		"""Frequency ride through and trip logic. """
		try:
			#Select frequency source
			fgrid = self.we/(2.0*math.pi)  #Use grid frequency as estimated by PLL
			if t > self.t_stable: #Go through logic only after a short time delay
			
				if not self.LFRT_TRIP: #Logic before tripping/momentary cessation
				
					if any(self.LFRT_dict[key]['t_LFstart']==0.0 for key in self.LFRT_dict.keys()):
					
						for LFRT_key,LFRT_values in self.LFRT_dict.items():
							if LFRT_values['t_LFstart'] == 0.0 and fgrid < LFRT_values['F_LF']:
								LFRT_values['t_LFstart']  = t
								self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='LF'+LFRT_key+'_start')
								#print(fgrid, LFRT_values['F_LF'],LFRT_values['t_LFstart'])

					if any(self.LFRT_dict[key]['t_LFstart'] > 0.0 for key in self.LFRT_dict.keys()):	
	  
						for LFRT_key,LFRT_values in self.LFRT_dict.items():
							if fgrid > LFRT_values['F_LF']  + self.del_f and LFRT_values['t_LFstart'] > 0.0: #Reset timer if freqency goes above
								self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='LF'+LFRT_key+'_reset')
								LFRT_values['t_LFstart']  = 0.0

						for LFRT_key,LFRT_values in self.LFRT_dict.items():
							if fgrid <= LFRT_values['F_LF'] and t-LFRT_values['t_LFstart'] >= LFRT_values['t_LF_limit'] and  LFRT_values['t_LFstart'] > 0.0: #Set HVRT_TRIP flag if timer exeeds limit
								self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='inverter_trip_LF'+LFRT_key)
								self.LFRT_TRIP = True
								LFRT_values['t_LFstart'] = 0.0
								self.t_LF_reconnect = 0.0

							elif fgrid <= LFRT_values['F_LF'] and t-LFRT_values['t_LFstart'] < LFRT_values['t_LF_limit'] and  LFRT_values['t_LFstart'] > 0.0: #Remain in LF zone and monitor
								self.print_LFRT_events(t,fgrid,LFRT_values['t_LFstart'],event_name='LF'+LFRT_key+'_zone')

				elif self.LFRT_TRIP: #Logic after tripping/momentary cessation
					self.DER_reconnect_logic(t)
	 
			else:
				self.LFRT_TRIP = False
				for LFRT_key,LFRT_values in self.LFRT_dict.items():
					 LFRT_values['t_LFstart'] = 0.0
		except:
			LogUtil.exception_handler()


	@property
	def FRT_INSTANTANEOUS_TRIP(self):
		try:
			return self.__FRT_INSTANTANEOUS_TRIP
		except:
			LogUtil.exception_handler()


	@FRT_INSTANTANEOUS_TRIP.setter
	def FRT_INSTANTANEOUS_TRIP(self,FRT_INSTANTANEOUS_TRIP):
		try:
			self.__FRT_INSTANTANEOUS_TRIP = FRT_INSTANTANEOUS_TRIP
			if FRT_INSTANTANEOUS_TRIP:			
				self.LFRT_dict['1']['t_LF_limit'] = self.LFRT_dict['2']['t_LF_limit']  = 30*(1/60) #Disconnect within n cycles					 
				self.HFRT_dict['1']['t_LF_limit'] = self.HFRT_dict['2']['t_LF_limit']  = 30*(1/60) #Disconnect within n cycles  
			
			else:
				self.LFRT_dict['1']['t_LF_limit'] = self.RT_config['LFRT']['1']['t_LF_limit'] 
				self.LFRT_dict['2']['t_LF_limit'] = self.RT_config['LFRT']['2']['t_LF_limit'] 
				self.HFRT_dict['1']['t_HF_limit']  = self.RT_config['HFRT']['1']['t_HF_limit'] 
				self.HFRT_dict['2']['t_HF_limit']  = self.RT_config['HFRT']['2']['t_HF_limit'] 
					
			return self.__FRT_INSTANTANEOUS_TRIP
		except:
			LogUtil.exception_handler()


	def print_LFRT_events(self,simulation_time,frequency,timer_start=0.0,event_name='',print_inline = False,verbose = False):
		"""Print LFRT events."""	
		try:
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
			elif event_name == 'LF1_zone' and verbose:
				text_string = '{time_stamp:.4f}:LF1 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
								.format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
			elif event_name == 'LF2_zone' and verbose:
				text_string = '{time_stamp:.4f}:LF2 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
								.format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
			elif event_name == 'LF3_zone' and verbose:
				text_string = '{time_stamp:.4f}:LF3 zone entered at:{timer_start:.4f}s and continuing for {time_elasped:.4f}s'\
								.format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
			elif event_name == 'inverter_trip_LF1':
				text_string = '{time_stamp:.4f}:LF1 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
								.format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
			elif event_name == 'inverter_trip_LF2':
				text_string = '{time_stamp:.4f}:LF2 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
								.format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)	
			elif event_name == 'inverter_trip_LF3':
				text_string = '{time_stamp:.4f}:LF3 violation at {time_stamp:.4f}s after {time_elasped:.4f} s for {frequency:.3f} Hz - Inverter will be tripped'\
								.format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)	
			elif event_name == 'reconnect_start':
				text_string = '{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s for {frequency:.3f} Hz'\
								.format(time_stamp=simulation_time,timer_start=timer_start,frequency=frequency)
			elif event_name == 'reconnect_reset':
				text_string = '{time_stamp:.4f}:Reconnect timer reset after {time_elasped:.4f} s for {frequency:.3f} Hz'\
								.format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
			elif event_name == 'reconnect_zone' and verbose:
				text_string = '{time_stamp:.4f}:Reconnect timer started at {timer_start:.4f} s and continuing for {time_elasped:.4f} s'\
								.format(time_stamp=simulation_time,timer_start=timer_start,time_elasped=simulation_time-timer_start)
			elif event_name == 'inverter_reconnection':
				text_string = '{time_stamp:.4f}:Inverter reconnecting after LF trip at {time_stamp:.4f}s after {time_elasped:.4f}s for {frequency:.3f} Hz'\
								.format(time_stamp=simulation_time,time_elasped=simulation_time-timer_start,frequency=frequency)
			elif event_name == 'inverter_tripped' and verbose: 
				text_string = '{time_stamp:.4f}:Inverter in tripped condition for {frequency:.3f} Hz'.format(time_stamp=simulation_time,frequency=frequency)
			else:
				text_string =''
			self.print_event(text_string,print_inline)
		except:
			LogUtil.exception_handler()


	def print_event(self,text_string,print_inline):
		"""Print information about ride through events."""
		try:
			if text_string != '':
				if print_inline:  #Print log in notebook window
					LogUtil.logger.info(text_string)

				else: #Print in console window
					utility_functions.print_to_terminal(text_string)
			else:
				pass
		except:
			LogUtil.exception_handler()


