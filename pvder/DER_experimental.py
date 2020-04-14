# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:25:23 2020

@author: splathottam
"""

class DERExperimental():
    """Class for experimental methods and features."""
  
    def update_inverter_states(self,iaR,iaI,xaR,xaI,uaI,uaR,xP,xQ,xPLL,wte):
        """Update inverter states"""
        
        Args:
             ia (complex): Inverter phase a current.
             xa (complex): Inverter controller state.
             ua (complex): Inverter controller state.
             Vdc (float): DC link voltage.             
        
        
        self.iaR = iaR
        self.iaI = iaI
        self.xaR = xaR
        self.xaI = xaI
        self.uaR = uaR
        self.uaI = uaI
        
        self.xP = xP
        self.xQ = xQ
        
        self.xPLL = xPLL
        self.wte = wte
        
        self.ia = self.convert2complex(iaR,iaI)
        self.xa = self.convert2complex(xaR,xaI)
        self.ua = self.convert2complex(uaR,uaI)
    """