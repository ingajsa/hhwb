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
        self._d_k_eff_t = None
        self._d_inc_t = None
        self._d_inc_sp_t = None
        self._d_con_t = None
        self._d_wb_t = 0.
        self._dt = 0.
        self._c_shock = 0
        
        self._aff = []

        self._k_eff_reco = np.empty(int(RECO_PERIOD*DT_STEP/TEMP_RES + 1))
        self._k_eff_reco[:] = np.nan
        self._inc_reco = np.empty(int(RECO_PERIOD*DT_STEP/TEMP_RES + 1))
        self._inc_reco[:] = np.nan
        self._inc_sp_reco = np.empty(int(RECO_PERIOD*DT_STEP/TEMP_RES + 1))
        self._inc_sp_reco[:] = np.nan
        self._cons_reco = np.empty(int(RECO_PERIOD*DT_STEP/TEMP_RES + 1))
        self._cons_reco[:] = np.nan
        self._wb_reco = np.empty(int(RECO_PERIOD*DT_STEP/TEMP_RES + 1))
        self._wb_reco[:] = np.nan

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

    @property
    def inc_reco(self):
        return self._inc_reco

    @property
    def inc_sp_reco(self):
        return self._inc_sp_reco

    @property
    def cons_reco(self):
        return self._cons_reco

    @property
    def wb_reco(self):
        return self._wb_reco

    @property
    def k_eff_reco(self):
        return self._k_eff_reco

    def update_reco(self, t_i=0., L_t=None, K=None):
        """
        Parameters
        ----------
        t : TYPE
            DESCRIPTION.
        t_i : TYPE
            DESCRIPTION.
        L_t : TYPE
            DESCRIPTION.
        K : TYPE
            DESCRIPTION.

        Returns
        -------
        None.
    
        """
        self._update_k_eff()
        self._update_income_sp(L_t, K)
        self._update_income()
        self._update_consum()
        self._update_wb()
        if t_i % TEMP_RES == 0:
            self._k_eff_reco[int(t_i/TEMP_RES)] = self._d_k_eff_t
            self._inc_reco[int(t_i/TEMP_RES)] = self._d_inc_t
            self._inc_sp_reco[int(t_i/TEMP_RES)] = self._d_inc_sp_t
            self._cons_reco[int(t_i/TEMP_RES)] = self._d_con_t
            self._wb_reco[int(t_i/TEMP_RES)] = self._d_wb_t
        self._t += self._dt
        return
    
    def shock(self, reco_period=RECO_PERIOD, temp_res=TEMP_RES, aff_flag=False,
              L=None, K=None, dt=None):
        """This function causes the household to be shocked. The recovery track is set up and the
           post disaster state is generated for all indicators.
        Parameters
        ----------
        reco_period : int, optional
            Time frame after first disaster in years where reconstruction is modeled (in years).
            The default is RECO_PERIOD.
        temp_res : TYPE, optional
            Temporal resolution after recovery in weeks. The default is TEMP_RES.
        aff_flag : bool, optional
            Bool indicating whether the household is affected by the current shock.
            The default is False.
        L : float, optional
            Total national damage.
        K : float, optional
            Total national capital stock.
        """
        self._dt = dt
        if aff_flag:
            self._t = 0
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

    
    def _get_reco_fee(self):
        """Helperfunction.
        Parameters
        ----------
        t : float, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        TYPE
            recovery fee.

        """
        return self._damage[self._c_shock] * np.e**(-self._t*self._lmbda)
