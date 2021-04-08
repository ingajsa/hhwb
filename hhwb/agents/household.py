#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:48:03 2020

@author: insauer
"""
import numpy as np
import matplotlib.pyplot as plt
from hhwb.agents.agent import Agent
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD, DT


AGENT_TYPE = 'HH'


class Household(Agent):
    """! Household definition. Computed from the FIES and interacts with the classes
    Government and Shock.
    Attributes:
        @param hhid (int): household id
        @param weight (float): household weight (summing up over all household weights returns the total
                        population of the administrative unit)
        @param vul (float): household vulnerability (independent of disaster magnitude)
        @param inc_0 (float): Total income in the predisaster situation
        @param inc_sp (float): Total income from social transfers in the predisaster situation
        @param con_0: predisaster consumption

        ***Attributes calculated for the predisaster situation (FIES) and national indicators***
        @param k_eff_0: predisaster effective capital stock
        @param lmbda: optimal reconstruction rate
        @param tau: time needed to achieve 95% recovery

        ***Attributes related to disaster recovery***
        @param ff (list): bools showing whether the household was affected by a shock or not
        @param damage (list): Experienced direct damage at each event
        @param d_k_eff_t (float): loss of effective capital stock at time step t after disaster
        @param d_inc_t (float): loss of income at timestemp t
        @param d_inc_sp_t (float): change in income from social transfers at time step t after disaster
        @param d_con_t (float): loss in consumption at time step t after disaster
        @param k_eff_reco (np.array): recovery track of capital stock
        @param inc_reco (np.array): recovery track of income
        @param inc_sp_reco (np.array): recovery track of income
        @param con_reco (np.array): recovery track of consumption
        @param wbls (float): aggregated wellbeing loss
    """

    def __init__(self, hhid=0, w=1., vul=0.2, i_0=1., i_sp=0.2, region=None,
                 savings=0., poverty_line=0., decile=None, isurban=1, ispoor=0):
        """! constructor"""

        Agent.__init__(self, AGENT_TYPE)
        #  Attributes set during initialisation
        self.__lmbda = [1.]
        self.__tau = []

        self.__hhid = hhid
        self.__weight = w
        self._vul = vul
        self.__inc_0 = i_0
        self.__inc_sp = i_sp
        self.__con_0 = i_0
        self.__sav_0 = savings
        self.__poverty_line = poverty_line
        self.__decile = decile
        self.__region = region
        self.__isurban = isurban
        self.__ispoor = ispoor

        self.__k_eff_0 = None

        self.__wbls = None

        #  self.__floor = None
        #  self.__tf = None

    @property
    def hhid(self):
        return self.__hhid

    @property
    def weight(self):
        return self.__weight

    @property
    def k_eff_0(self):
        return self.__k_eff_0

    @property
    def lmbda(self):
        return self.__lmbda

    @property
    def tau(self):
        return self.__tau

    @property
    def income_0(self):
        return self.__inc_0

    @property
    def income_sp(self):
        return self.__inc_sp

    @property
    def consum_0(self):
        return self.__con_0

    @property
    def decile(self):
        return self.__decile



    def set_tax_rate(self, tax_rate=0):
        """
        Prepares the recovery process by calculating the tax rate the household's
        effective capital stock .
        """
        self.__tax_rate = tax_rate
        self.__k_eff_0 = (self.__inc_0 - self.__inc_sp)/((1-self.__tax_rate)*PI)


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
        if aff_flag:
            
            self._d_k_eff_t += (self.__k_eff_0 - self._d_k_eff_t) * self._vul
            opt_vul = self._d_k_eff_t / self.__k_eff_0

            self._damage.append(self._d_k_eff_t)
        
            if len(self._damage) > 2:
                self.__optimize_reco(vul=opt_vul)
            else:
                self.__optimize_reco(self._vul)
        self._d_inc_sp_t = (L/K) * self.__inc_sp
        self._d_inc_t = (1-self.__tax_rate) * PI * self._d_k_eff_t + self._d_inc_sp_t
        self._d_con_t = self._d_inc_t + self.__lmbda[self._c_shock] * self._get_reco_fee()
        self._update_wb()

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
        return self._damage[self._c_shock] * np.e**(-self._t*self.__lmbda[self._c_shock])

    def set_k_eff_0(self):
        self._k_eff_0 = (self.__inc_0 - self.__inc_sp)/((1-self.__tax_rate)*PI)

    def _update_income_sp(self, L_t, K):
        self._d_inc_sp_t = (L_t/K) * self.__inc_sp
        return

    def _update_income(self):
        self._d_inc_t = (1-self.__tax_rate) * PI * self._d_k_eff_t + self._d_inc_sp_t
        return

    def _update_consum(self):
        self._d_con_t = self._d_inc_t + self.__lmbda[self._c_shock] * self._get_reco_fee()
        return

    def _update_k_eff(self):

        self._d_k_eff_t = self._get_reco_fee()
        return
    
    def _update_wb(self):

        self._d_wb_t += (self.__con_0**(1 - ETA))/(1 - ETA) *\
                      (1-((1 - (self._d_con_t / self.__con_0))**(1-ETA))) *\
                      np.e**(-RHO * self._t)


    def __smooth_with_savings(self):
        """Sets the floor taken from savings to smoothen HH's consumption
        loss and the time tf when it runs out of savings
        """

        dc0 = self.__k_eff_0 * self.__vul*(PI+self.__lmbda)

        if dc0 == 0:
            self.__floor = 0
            self.__tf = T_RNG
        if self.__lmbda == 0:
            self.__floor = int(round(min(dc0, max(dc0-(2/3)
                                                * self.sav, 0.)), 0))
            self.__tf = 1.

        gamma = dc0
        last_result = None

        while True:

            beta = gamma/dc0
            result = dc0 * (1-beta) + gamma * np.log(beta)\
                - self.sav*self.__lmbda

            if (last_result < 0 and result > 0) or\
               (last_result > 0 and result < 0):

                _t = -np.log(beta)/self.__lmbda

            if _t < 0:
                print('RESULT!:\ngamma = ', gamma, '& beta = ',
                      beta, ' & t = ', _t)
                print('CHECK:', dc0 * np.e**(self.__lmbda * _t),
                      ' gamma = ', gamma)

            if _t >= T_RNG:
                self.__floor = int(round(min(dc0, max(dc0-(2/3)
                                                      * self.sav, 0.)), 0.))
                self.__tf = 1.

            self.__floor = int(round(gamma, 0))
            self.__tf = round(_t, 3)

            last_result = result
            gamma -= 0.01 * dc0
            if gamma <= 0:
                self.__floor = 0
                self.__tf = T_RNG


    def __optimize_reco(self, vul=0.3):
        """
        This is the core optimization function, that numerically optimizes
        the optimal reconstruction rate lmbda of the household and derives the
        reconstruction time tau
        TODO (- eventually implement a static jit version
              - this must be done multicore)
        """
        if vul == 0:
            return

        last_integ = None
        last_lambda = None

        lmbda = 0.0

        print('vul='+str(vul))

        c=0

        while True:

            integ = 0.0
            for dt in np.linspace(0, T_RNG, DT_STEP*T_RNG):
                integ += np.e**(-dt * (RHO + lmbda)) * ((PI + lmbda) * dt - 1) * (PI - (PI + lmbda) * vul * np.e**(-lmbda * dt))**(-ETA)
                
            if last_integ and ((last_integ < 0 and integ > 0) or
                                (last_integ > 0 and integ < 0)):
                print('\n Found the Minimum!\n lambda = ', last_lambda,
                      '--> integ = ', last_integ)
                print('lambda = ', lmbda, '--> integ = ', integ, '\n')

                out = (lmbda+last_lambda)/2

                self.__lmbda.append(out)
                self.__tau = np.log(1./0.05) * (1./out)
                return

            last_integ = integ
            if last_integ is None:
                assert(False)

            last_lambda = lmbda
            lmbda += 0.01
            c+=1

    def plot_reco_trajec(self, timeframe=40, pred=5):

        pre_dis = np.linspace(-pred, 0, pred*DT_STEP)
        pre_cap_stock = np.ones(pred*DT_STEP)*100
        time = np.linspace(0, timeframe, timeframe*DT_STEP)
        cap_stock = (1-np.e**(-self.__lmbda*time))*100
        time = np.concatenate((pre_dis, time), axis=None)
        cap_stock = np.concatenate((pre_cap_stock, cap_stock), axis=None)
        plt.plot(time, cap_stock, label='HH'+str(int(self.__hhid)))
        plt.xlabel('years after shock')
        plt.ylabel('% recovered')
        plt.title('Optimal recovery track')
        plt.legend()
        return
    
    def plot_life_trajec(self, timeframe=40, pred=5):

        pre_dis = np.linspace(-pred, 0, pred*DT_STEP)
        pre_cap_stock = np.ones(pred*DT_STEP)*100
        time = np.linspace(0, timeframe, timeframe*DT_STEP)
        cap_stock = (1-np.e**(-self.__lmbda*time))*100
        time = np.concatenate((pre_dis, time), axis=None)
        cap_stock = np.concatenate((pre_cap_stock, cap_stock), axis=None)
        plt.plot(time, cap_stock, label='HH'+str(int(self.__hhid)))
        plt.xlabel('years after shock')
        plt.ylabel('% recovered')
        plt.title('Optimal recovery track')
        plt.legend()
        return
