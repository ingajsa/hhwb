#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:48:03 2020

@author: insauer
"""
import numpy as np
import matplotlib.pyplot as plt
from hhwb.agents.agent import Agent
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD


AGENT_TYPE = 'HH'


class Household(Agent):
    """Household definition. Computed from the FIES and interacts with the classes
    Government and Shock.

    Attributes:
        
        ***Attributes extracted from the FIES***
        hhid (int): household id
        weight (float): household weight (summing up over all household weights returns the total
                        population of the administrative unit)
        vul (float): household vulnerability (independent of disaster magnitude)
        inc_0 (float): Total income in the predisaster situation
        inc_sp (float): Total income from social transfers in the predisaster situation
        con_0: predisaster consumption

        ***Attributes calculated for the predisaster situation (FIES) and national indicators***
        k_eff_0: predisaster effective capital stock
        lmbda: optimal reconstruction rate
        tau: time needed to achieve 95% recovery

        ***Attributes related to disaster recovery***
        aff (list): bools showing whether the household was affected by a shock or not
        damage (list): Experienced direct damage at each event
        d_k_eff_t (float): loss of effective capital stock at time step t after disaster
        d_inc_t (float): loss of income at timestemp t
        d_inc_sp_t (float): change in income from social transfers at time step t after disaster
        d_con_t (float): loss in consumption at time step t after disaster
        k_eff_reco (np.array): recovery track of capital stock
        inc_reco (np.array): recovery track of income
        inc_sp_reco (np.array): recovery track of income
        con_reco (np.array): recovery track of consumption
        wbls (float): aggregated wellbeing loss
    """

    def __init__(self, hhid=0, w=1., vul=0.2, i_0=1., i_sp=0.2):
        

        Agent.__init__(self, AGENT_TYPE)
        #  Attributes set during initialisation
        self.__hhid = hhid
        self.__weight = w
        self.__vul = vul
        self.__inc_0 = i_0
        self.__inc_sp = i_sp
        self.__con_0 = i_0

        self.__k_eff_0 = None
        self.__lmbda = None
        self.__tau = None
        
        self.__aff = []
        self.__damage = []
        # share of income_0 coming from social transfers
        self.__d_k_eff_t = None
        self.__d_inc_t = None
        self.__d_inc_sp_t = None
        self.__d_con_t = None

        self.__k_eff_reco = np.array([])
        self.__inc_reco = np.array([])
        self.__inc_sp_reco = np.array([])
        self.__con_reco = np.array([])

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
    def affected(self):
        return self.__aff

    @property
    def k_eff_t(self):
        return self.__d_k_eff_t

    @property
    def lmbda(self):
        return self.__lmbda

    @property
    def vul(self):
        return self.__vul

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
    def d_inc_t(self):
        return self.__d_inc_t

    @property
    def d_inc_sp_t(self):
        return self.__d_inc_sp_t

    @property
    def d_con_t(self):
        return self.__d_con_t

    @property
    def k_eff_0(self):
        return self.__k_eff_0

    def set_tax_rate(self, tax_rate=0):
        self.__tax_rate = tax_rate

    def init_reco(self):
        """
        Prepares the recovery process by calculating the household's effective capital stock and
        the optimal recovery track.
        """
        self.__k_eff_0 = (self.__inc_0 - self.__inc_sp)/((1-self.__tax_rate)*PI)
        self.__optimize_reco()

    def update_reco(self, t, L_t, K):
        self.__update_k_eff(t)
        self.__update_income_sp(L_t, K)
        self.__update_income()
        self.__update_consum(t)
        if t % TEMP_RES == 0:
            self.__k_eff_reco[int(t/TEMP_RES)] == self.__d_k_eff_t
            self.__inc_reco[int(t/TEMP_RES)] == self.__d_inc_t
            self.__inc_sp_reco[int(t/TEMP_RES)] == self.__d_inc_sp_t
            self.__cons_reco[int(t/TEMP_RES)] == self.__d_con_t
        return

    def shock(self, reco_period=RECO_PERIOD, temp_res=TEMP_RES, aff_flag=False,
              L=0, K=0):
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
            Total national damage. The default is 0.
        K : float, optional
            Total national capital stock. The default is 0.
        """

        self.__aff.append(aff_flag)
        #  initialize recovery track
        self.__inc_reco = np.zeros(RECO_PERIOD*DT_STEP/TEMP_RES)
        self.__inc_sp_reco = np.zeros(RECO_PERIOD*DT_STEP/TEMP_RES)
        self.__cons_reco = np.zeros(RECO_PERIOD*DT_STEP/TEMP_RES)
        self.__k_eff_reco = np.zeros(RECO_PERIOD*DT_STEP/TEMP_RES)
        #  set the affected state
        self.__set_shock_state(L, K)

        return
    
    def __set_shock_state(self, L, K):

        self.__d_k_eff_t = self.__d_k_eff_t * self.__vul
        self.__damage.append(self.__d_k_eff_t)
        self.__d_inc_sp_t = (L/K) * self.__inc_sp
        self.__d_inc_t = (1-self.__tax_rate) * PI * self.__d_k_eff_t + self.__d_inc_sp_t
        self.__d_con_t = self.__lmbda * self.__get_reco_fee()

        return

    def __get_reco_fee(self, t=0.):

        return self.__damage * np.e**(-t*self.__lmbda)
    
    def set_k_eff_0(self):
        self.__k_eff_0 = (self.__inc_0 - self.__inc_sp)/((1-self.__tax_rate)*PI)
    
    def __update_income_sp(self, L_t, K):
        self.__d_inc_sp_t = (1-self.__tax_rate) * PI * self.__d_k_eff_t + self.__d_inc_sp_t
        return

    def __update_income(self):
        self.__d_inc_t = (1-self.__tax_rate) * PI * self.__d_k_eff_t + self.__d_inc_sp_t
        return

    def __update_consum(self, t):
        self.__d_con_t = self.__lmbda * self.__get_reco_fee(t=t)
        return

    def __update_k_eff(self, t):

        self.__d_k_eff_t = self.__get_reco_fee(t=t)
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
