#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:48:03 2020

@author: insauer
"""
import numpy as np
from hhwb.agents.agent import Agent
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP

AGENT_TYPE = 'HH'


class Household(Agent):

    def __init__(self):

        Agent.__init__(self, AGENT_TYPE)
        self.vul = 0.2
        self.k_eff_0 = 100
        self.income = 0.
        self.consum = 0.
        self.wellfare = 0.
        """TODO: implement initial values for capital stock, consumption,
           savings coming from fies"""
        self.lmbda = None
        self.floor = None
        self.tf = None

    def update(self):
        self.__update_k_eff()
        self.__update_income()
        self.__update_consum()

        return

    def get_shock(self):
        self.__optimize_reco()
        self.__smooth_with_savings()
        return

    def __update_income(self):
        return

    def __update_consumption(self):
        return

    def __update_k_eff(self):
        return

    def __optimize_reco(self):
        """
        This is the core optimization function, that numerically optimizes
        the optimal reconstruction rate lmbda of the household and derives the
        reconstruction time tau
        TODO (- eventually implement a static jit version
              - this must be done multicore)
        """
        if self.vul == 0:
            return

        last_integ = None
        last_lambda = None

        lmbda = 0.0

        while True:

            integ = 0
            for dt in range(0, T_RNG, DT_STEP*T_RNG):
                integ += np.e**(-dt * (RHO + lmbda)) * ((PI + lmbda) * dt - 1)\
                         * (PI - (PI + lmbda) * self.vul *
                            np.e**(-lmbda * dt))**(-ETA)

            if last_integ and ((last_integ < 0 and integ > 0) or
                               (last_integ > 0 and integ < 0)):
                print('\n Found the Minimum!\n lambda = ', last_lambda,
                      '--> integ = ', last_integ)
                print('lambda = ', lmbda, '--> integ = ', integ, '\n')

                out = (lmbda+last_lambda)/2

                self.lmbda = out
                self.tau = np.log(1./0.05) * (1./self.lmbda)
                return

            last_integ = integ
            if last_integ is None:
                assert(False)

            last_lambda = lmbda
            lmbda += 0.01

        def __smooth_with_savings(self):
            """Sets the floor taken from savings to smoothen HH's consumption
            loss and the time tf when it runs out of savings
            """

            dc0 = self.k_eff_0 * self.vul*(PI+self.lmbda)

            if dc0 == 0:
                self.floor = 0
                self.tf = T_RNG
            if self.lmbda == 0:
                self.floor = int(round(min(dc0, max(dc0-(2/3)
                                                    * self.sav, 0.)), 0))
                self.tf = 1.

            gamma = dc0
            last_result = None

            while True:

                beta = gamma/dc0
                result = dc0 * (1-beta) + gamma * np.log(beta)\
                    - self.sav*self.lmbda

                if (last_result < 0 and result > 0) or\
                   (last_result > 0 and result < 0):

                    _t = -np.log(beta)/self.lmbda

                if _t < 0:
                    print('RESULT!:\ngamma = ', gamma, '& beta = ',
                          beta, ' & t = ', _t)
                    print('CHECK:', dc0 * np.e**(self.lmbda * _t),
                          ' gamma = ', gamma)

                if _t >= T_RNG:
                    self.floor = int(round(min(dc0, max(dc0-(2/3)
                                                        * self.sav, 0.)), 0.))
                    self.tf = 1.

                self.floor = int(round(gamma, 0))
                self.tf = round(_t, 3)

                last_result = result
                gamma -= 0.01 * dc0
                if gamma <= 0:
                    self.floor = 0
                    self.tf = T_RNG
