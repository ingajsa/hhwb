#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 21:41:02 2020

@author: insauer
"""

import numpy as np
from hhwb.agents.agent import Agent
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD

AGENT_TYPE = 'GV'


class Government(Agent):

    """Household definition. Computed from the FIES and interacts with the classes
    Government and Shock.

    Attributes:
        ***Attributes regarding the predisaster situation***
        tax_rate (float): flat income tax
        sp_cost (float): total spendings on social transfers
        K (float): Total national capital stock

        ***Attributes related to disaster recovery***
        L_0 (list): total national damage experienced at each shock
        vul (list): conceptional vulnerability at each shock as ratio L/K
        lmbda (float): optimal recovery rate
        tau (float): time at which 95% of capital stock is recovered
        L_t (float): loss at timestep t after disaster
    """

    def __init__(self):
        """Empty Definition"""
        Agent.__init__(self, AGENT_TYPE)
        self.__tax_rate = 0.
        self.__sp_cost = 0.
        self.__K = 0.

        self.__L_0 = 0
        self.__vul = []

        self.__lmbda = None
        self.__tau = None

        self.__L_t = 0.



    @property
    def tax_rate(self):
        return self.__tax_rate

    @property
    def sp_cost(self):
        return self.__sp_cost

    @property
    def K(self):
        return self.__K

    @property
    def L_t(self):
        return self._d_k_eff_t
    
    @property
    def con_smooth(self):
        return self.__con_smooth
    
    @property
    def wb_smooth(self):
        return self.__wb_smooth




    def set_tax_rate(self, reg_hh):
        """
        Sets the optimal tax rate basing on the total enrollment for social transfers and
        extracts capital stock from each hh to get total national capital stock
            Parameters:
                reg_hh (list): all registered households
        """

        tot_inc = 0.
        tot_sp = 0.
        # get total income and total social transfers from all households
        for hh in reg_hh:
            tot_inc += hh.weight*hh.income_0
            tot_sp += hh.weight*hh.income_sp

        # set governments total expenditure on social transfers and required tax rate
        self.__sp_cost = tot_sp
        self.__tax_rate = tot_sp/tot_inc
        # set tax rate for all households and get capital stock
        for hh in reg_hh:
            hh.set_tax_rate(self.__tax_rate)
            hh.init_life()
            self.__K += hh.weight*hh.k_eff_0
        return
    
    def update_reco(self, t_i=0., hh_reg=None):
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

        self.__update_all(hh_reg)

        return


    # def update(self, reg_hh):
    #     self.__L_t = 0.
    #     self.__d_cons = 0.
    #     for hh in reg_hh:
    #         self.__L_t += hh.d_k_eff_t
    #     return

    def __update_all(self, reg_hh):
        self._d_k_eff_t = 0.
        self._d_con_t = 0.
        self._d_inc_t = 0.
        self._d_inc_sp_t = 0.
        self._d_wb_t = 0.
        self.__con_smooth = 0.
        self.__wb_smooth = 0.
        for hh in reg_hh:
            self._d_k_eff_t += hh.weight* hh.d_k_eff_t
            self._d_con_t += hh.weight*hh.d_con_t
            self._d_inc_t += hh.weight*hh.d_inc_t
            self._d_inc_sp_t += hh.weight*hh.d_inc_sp_t
            self._d_wb_t += hh.weight*hh.d_wb_t
            self.__wb_smooth += hh.weight*hh.wb_smooth
            self.__con_smooth += hh.weight*hh.con_smooth
        return
    
    def _set_shock_state(self, L, K, aff_flag):
        """This function calculates the initial damage and the initial loss in income and
           consumption.
        Parameters
        ----------
        L : float, optional
            Total national damage. The default is 0.
        K : float, optional
            Total national capital stock. The default is 0.
        """
        #self._dt = dt

        self._vul = L/self.__K
        #self._optimize_reco()
        if aff_flag:
            self._d_k_eff_t = L
        else:
            self._d_k_eff_t = 0

        self._damage.append(self._d_k_eff_t)
        self._d_inc_sp_t = (L/self.__K) * self.__sp_cost
        self._d_inc_t = PI * L + self._d_inc_sp_t
        self._d_con_t = np.nan

        return

    # def _update_income_sp(self):
    #     self._d_inc_sp_t = (self._d_k_eff_t/self.__K) * self.__sp_cost
    #     return

    # def _update_income(self):
    #     self._d_inc_t = PI * self._d_k_eff_t
    #     return

    # def _update_consum(self, t):
    #     self._d_con_t = self._d_inc_t + self._lmbda * self._get_reco_fee(t=t)
    #     return

    # def _update_k_eff(self, t):

    #     self._d_k_eff_t = self._get_reco_fee(t=t)
    #     return
