#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 21:41:02 2020

@author: insauer
"""

import numpy as np
from hhwb.agents.agent import Agent
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP

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

        self.__L_0 = []
        self.__vul = []

        self.__lmbda = None
        self.__tau = None

        self.__L_t = np.array([])

    @property
    def tax_rate(self):
        return self.__tax_rate

    @property
    def sp_cost(self):
        return self.__sp_cost

    @property
    def K(self):
        return self.__K

    def update(self, hh):

        return

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
            hh.init_reco()
            self.__K += hh.weight*hh.k_eff_0
        return


    def start_hh_reco(self):
        self.__vul = self.__L_t/self.__K

        return

    def __update_L_t(self, reg_hh):
        for hh in reg_hh:
            self.__L_t += hh.d_k_eff_t
        return
    
    

    def __optimize_reco(self):
        """
        This is the core optimization function, that numerically optimizes
        the optimal reconstruction rate lmbda of the household and derives the
        reconstruction time tau
        TODO (- eventually implement a static jit version
              - this must be done multicore)
        """
        if self.__vul == 0:
            return
    
        last_integ = None
        last_lambda = None
    
        lmbda = 0.0
    
        while True:
    
            integ = 0.0
            for dt in np.linspace(0, T_RNG, DT_STEP*T_RNG):
                integ += np.e**(-dt * (RHO + lmbda)) * ((PI + lmbda) * dt - 1) * (PI - (PI + lmbda) * self.__vul * np.e**(-lmbda * dt))**(-ETA)
    
            if last_integ and ((last_integ < 0 and integ > 0) or
                               (last_integ > 0 and integ < 0)):
                print('\n Found the Minimum!\n lambda = ', last_lambda,
                      '--> integ = ', last_integ)
                print('lambda = ', lmbda, '--> integ = ', integ, '\n')
    
                out = (lmbda+last_lambda)/2
    
                self.__lmbda = out
                self.__tau = np.log(1./0.05) * (1./self.__lmbda)
                return
    
            last_integ = integ
            if last_integ is None:
                assert(False)
    
            last_lambda = lmbda
            lmbda += 0.01
