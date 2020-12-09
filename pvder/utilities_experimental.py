# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:28:24 2020

@author: splathottam
"""

from graphviz import Digraph

from pvder.logutil import LogUtil


class ModelUtilities():
	"""Class for model wide utilities."""

	def __init__(self,PV_model,simulation,grid_model=None):
		try:
			self.PV_model = PV_model
			self.simulation = simulation
			if self.PV_model.standAlone and grid_model is not None:
				self.grid_model = grid_model
			elif self.PV_model.standAlone and grid_model is None:
				raise ValueError('`Grid` instance need to provided in stand alone mode for creating `GridSimulation` instance!')
		except:
			LogUtil.exception_handler()


	def draw_model(self,display_value_type='per_unit'):
		"""Draw and render the model using graphs."""
		try:
			assert display_value_type in ['per_unit','actual'],'Display can only be in per unit or actual values'
			self.display_value_type =display_value_type
			if self.display_value_type  == 'per_unit': 
				_C = self.PV_model.C
				_Lf = self.PV_model.Lf
				_Rf = self.PV_model.Rf
				_Z1 = self.PV_model.Z1
				_Z2 = self.grid_model.Z2
			else: 
				_C = self.PV_model.C_actual
				_Lf = self.PV_model.Lf_actual
				_Rf = self.PV_model.Rf_actual
				_Z1 = self.PV_model.Z1_actual
				_Z2 = self.grid_model.Z2_actual

			dot = Digraph(comment='PV_DER and grid_model.')
			dot.node('Base','Vbase={},Sbase={},Zbase={:.4f},Lbase={:.4f},Cbase={:.4f}'.format(self.grid_model.Vbase,self.grid_model.Sbase,self.grid_model.Zbase,self.grid_model.Lbase,self.grid_model.Cbase),shape='rectangle')
			dot.node('Value_type','{}'.format(self.display_value_type))
			dot.node('DER','{}\nC={:.5f},Lf={:.5f},Rf= {:.4f}'.format(self.PV_model.name,_C,_Lf,_Rf),shape='box')
			dot.node('Transformer','{}\nZl = {:.5f}'.format(self.PV_model.transformer_name,_Z1),shape='rectangle')
			dot.node('Transmission_Line','{}\nZt = {:.5f}'.format(self.grid_model.transmission_name,_Z2),shape='rectangle')
			dot.node('Grid',self.grid_model.name)
			
			dot.edge('DER', 'Transformer')
			dot.edge('Transformer', 'Transmission_Line')
			dot.edge('Transmission_Line', 'Grid')
			dot.render('model_graphs/model.gv', view=True)
		except:
			LogUtil.exception_handler()


