#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 23:19:16 2020

@author: insauer
"""
import numpy as np
import matplotlib.pyplot as plt
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD

class Agent():
    
    def __init__(self, agent_type):
        
        self.agent_type = agent_type
        
        self._damage = []
        # share of income_0 coming from social transfers
        self._d_k_eff_t = 0.
        self._d_inc_t = 0.
        self._d_inc_sp_t = 0.
        self._d_con_t = 0.
        self._d_wb_t = 0.
        self._dt = 0.
        self._c_shock = 0
        
        self._aff = []


    @property
    def dt(self):
        return self._dt
    
    def t(self):
        return self._t

    @property
    def vul(self):
        return self._vul

    @property
    def tau(self):
        return self._tau

    @property
    def d_k_eff_t(self):
        return self._d_k_eff_t

    @property
    def d_inc_t(self):
        return self._d_inc_t

    @property
    def d_inc_sp_t(self):
        return self._d_inc_sp_t

    @property
    def d_con_t(self):
        return self._d_con_t

    @property
    def d_wb_t(self):
        return self._d_wb_t

    @property
    def affected(self):
        return self.__aff

    def set_vul(self, vul):
        self._vul = vul

    def shock(self, aff_flag=False,
              L=None, K=None, dt=None):
        """This function causes the household to be shocked. The recovery track is set up and the
           post disaster state is generated for all indicators.
        Parameters
        ----------
        aff_flag : bool, optional
            Bool indicating whether the household is affected by the current shock.
            The default is False.
        L : float, optional
            Total national damage.
        K : float, optional
            Total national capital stock.
        """
        self._dt = dt
        #print(self._c_shock)
        #print(aff_flag)
        if aff_flag:
            self._c_shock += 1
        self._aff.append(aff_flag)
        #  set the affected state
        self._set_shock_state(L, K, aff_flag)

        return

    def init_life(self):
        """This initialises the life.
        """
        self._t = 0.
        self._damage.append(0.)
        self._d_k_eff_t = 0.
        self._d_inc_sp_t = 0.
        self._d_inc_t = 0.
        self._d_con_t = 0.

        return

    
    # def _get_reco_fee(self):
    #     """Helperfunction.
    #     Parameters
    #     ----------
    #     t : float, optional
    #         DESCRIPTION. The default is 0.

    #     Returns
    #     -------
    #     TYPE
    #         recovery fee.

    #     """
    #     return self._damage[self._c_shock] * np.e**(-self._t*self._lmbda)
