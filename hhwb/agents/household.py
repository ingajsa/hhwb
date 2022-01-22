#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 21:48:03 2020

@author: insauer
"""
import numpy as np
import matplotlib.pyplot as plt
from hhwb.agents.agent import Agent
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, DT, SUBS_SAV_RATE, OPT_DYN


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

    def __init__(self, hhid=0, n_inds=1, w=1., vul=0.2, i_0=1., i_sp=0.2, region=None,
                 savings=0., subsistence_line=0., decile=None, isurban=1, ispoor=0):
        """! constructor"""

        Agent.__init__(self, AGENT_TYPE)
        #  Attributes set during initialisation
        self.__lmbda = [1.]
        self.__tau = []

        self.__hhid = hhid
        self.__n_inds = n_inds
        self.__weight = w
        self._vul = vul
        self.__inc_0 = i_0
        self.__inc_sp = i_sp
        self.__con_0 = i_0
        self.__sav_0 = savings
        self.__sav_t = savings
        self.__subsistence_line = subsistence_line
        self.__decile = decile
        self.__region = region
        self.__isurban = isurban
        self.__ispoor = ispoor
        self.__recovery_type = 0
        self.__poverty_trap = False
        self.__recovery_spending = 0.
        self.__con_smooth = 0
        
        self.__tf = -0.01
        self.__floor = 0.

        self.__k_eff_0 = None
        self.__twb = 0.
        self.__wb_0 = 0.
        self.__wb_t = 0.
        self.__con_smooth = 0.
        self.__wb_smooth = 0.
        self.__wb_0_sm = 0.
        self.__wb_t_sm = 0.
        self.__wbls = 0.
        self.__cnt_ts = 0
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
    def con_smooth(self):
        return self.__con_smooth
    
    @property
    def wb_smooth(self):
        return self.__wb_smooth

    @property
    def decile(self):
        return self.__decile

    @property
    def poverty_trap(self):
        return self.__poverty_trap
    
    @property
    def recovery_type(self):
        return self.__recovery_type

    @property
    def subsistence_line(self):
        return self.__subsistence_line

    def update_reco(self, L_t=None, K=None):
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
        self._update_wb_sav()
        self.__update_savings()
        # if self.__recovery_type == 1:
        #     self._update_wb_sav()
        return



    def set_tax_rate(self, tax_rate=0):
        """
        Prepares the recovery process by calculating the tax rate the household's
        effective capital stock.
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
        #print(self.__hhid)
        if aff_flag:
            
            if (self.__k_eff_0 - self._d_k_eff_t) > 0.0:
                self._t = 0.
                self._d_k_eff_t += (self.__k_eff_0 - self._d_k_eff_t) * self._vul
                opt_vul = self._d_k_eff_t / self.__k_eff_0
                self._damage.append(self._d_k_eff_t)
                self.__optimize_reco(vul=opt_vul)
                
            else:
                #self.__poverty_trap = True
                self._damage.append(self._d_k_eff_t)
                self.__lmbda.append(self.__lmbda[self._c_shock-1])
                self.update_reco(L_t=L, K=K)
                return
            #if not self.__poverty_trap:
            self._d_inc_sp_t = (L/K) * self.__inc_sp
            self._d_inc_t = ((1-self.__tax_rate) * PI * self._d_k_eff_t + self._d_inc_sp_t)/DT_STEP
    
            self._set_recovery_path()
        
            if self.__recovery_type == 1:
                self.__cnt_ts += 1
                self._d_con_t = self._d_inc_t + self.__recovery_spending
                self.__smooth_with_savings_1()
                
                self.__con_smooth = self.__floor
                
            elif self.__recovery_type == 2:
                self._d_con_t = self._d_inc_t + self._possible_reco()
                self.__smooth_with_savings_2()
                self.__con_smooth = self.__floor
            
            elif self.__recovery_type == 3:
                self._d_con_t = self._d_inc_t
                self.__smooth_with_savings_3()
                self.__con_smooth = self.__floor
            
            elif self.__recovery_type == 4:
                self._d_con_t = self._d_inc_t
                self.__smooth_with_savings_3()
                self.__con_smooth = self.__floor
            self._update_wb()
            self._update_wb_sav()

        else:
            
            self.update_reco(L_t=L, K=K)
            
        # if self.__recovery_type == 1:
        #     self._update_wb_sav()
        #self.__determine_sav_opt()
            # else:
            #     self._d_inc_sp_t = (L/K) * self.__inc_sp
            #     self._d_inc_t = self._d_inc_t - (self.__inc_sp-self._d_inc_sp_t)
            #     self._d_con_t = self.__con_0/100.#self._d_inc_t
            #     self._update_wb()
        return

    def _get_reco_fee(self, t1=None, t2=None):
        """
        Returns
        -------
        TYPE
            recovery fee.
        """
        if not t1:
            t1 = self._t

        if not t2:
            t2 = self._t + self._dt

        dam_0 = self._damage[self._c_shock] * np.e**(-t1*self.__lmbda[self._c_shock])
        dam_1 = self._damage[self._c_shock] * np.e**(-t2*self.__lmbda[self._c_shock])
        return dam_0-dam_1
    
    def _set_recovery_path(self):
        """Check for recovery below subsistence level.
        Parameters
        ----------
        t : float, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        TYPE
            recovery fee.
        """
        # hh not effected --> 0
        if self._d_k_eff_t <= 0:
            self.__recovery_type = 0
            self.__recovery_spending = 0.
            return
        # optimal_recovery_spending = self._get_reco_fee()
        # self.__recovery_type = 1
        # self.__recovery_spending = optimal_recovery_spending
        # return

        # optimal recovery --> 1
        # second best recovery --> 2
        # recovery generally below substistence --> 3
        # recovery starting below substistence --> 4

        optimal_recovery_spending = self._get_reco_fee()
        if self._check_subs(optimal_recovery_spending) > 0:
            self.__recovery_type = 1
            self.__recovery_spending = optimal_recovery_spending
            return
        elif self._possible_reco() > 0:
            self.__recovery_type = 2
            self.__recovery_spending = self._possible_reco() + SUBS_SAV_RATE

            return

        elif self.__con_0 < self.__subsistence_line:
            self.__recovery_type = 3
            self.__recovery_spending = SUBS_SAV_RATE
            return
        else:
            self.__recovery_type = 4
            self.__recovery_spending = SUBS_SAV_RATE
            return

    def _check_subs(self, opt_reco):

        cons_level = self.__con_0 - (self._d_inc_t + opt_reco)

        if cons_level >= self.__subsistence_line:

            return opt_reco

        else:
            return -1

    def _possible_reco(self):

        possible_reco = self.__con_0 - self._d_inc_t - self.__subsistence_line

        if possible_reco > self._d_k_eff_t:
            possible_reco = self._d_k_eff_t

        if possible_reco > 0:
            return possible_reco
        else:
            return -1

    def _update_reco_spend(self):


        if self.__recovery_type == 0:
            self.__recovery_spending = 0.
            return

        if self.__recovery_type == 1:
            if OPT_DYN:
                opt_vul = self._d_k_eff_t / self.__k_eff_0
                self._t = 0.
                del self.__lmbda[-1]
                self.__optimize_reco(vul=opt_vul)
                self.__cnt_ts += 1
                del self._damage[-1]
                self._damage.append(self._d_k_eff_t)
                self.__recovery_spending = self._get_reco_fee()
                self.__smooth_with_savings_1()
            else:
                self.__recovery_spending = self._get_reco_fee()
            return

        if self.__recovery_type == 2:

            opt_vul = self._d_k_eff_t / self.__k_eff_0
            del self.__lmbda[-1]
            self.__optimize_reco(vul=opt_vul)
            del self._damage[-1]
            self._damage.append(self._d_k_eff_t)
            self._t = 0.
            optimal_recovery_spending = self._get_reco_fee()
            if self._check_subs(optimal_recovery_spending) > 0:
                self.__recovery_type = 1
                self.__smooth_with_savings_1()
                self.__recovery_spending = optimal_recovery_spending
                self.__con_smooth = self.__floor
            else:
                self.__recovery_spending = self._possible_reco() + SUBS_SAV_RATE
                self.__smooth_with_savings_2()
            return

        if self.__recovery_type == 3:

            if self._d_k_eff_t <= 0:
                self.__recovery_spending = 0
                self.__recovery_type = 0
                return
            elif self._d_k_eff_t < SUBS_SAV_RATE:
                self.__recovery_spending = self._d_k_eff_t
            else:
                self.__recovery_spending = SUBS_SAV_RATE
            return

        if self.__recovery_type == 4:
            if self._possible_reco() > 0:
                self._t = 0.
                self.__recovery_type = 2
                self.__recovery_spending = self._possible_reco() + SUBS_SAV_RATE
                self.__smooth_with_savings_2()
            else:
                
                if self._d_k_eff_t <= 0:
                    self.__recovery_spending = 0
                    self.__recovery_type = 0
                    return

                if self._d_k_eff_t < SUBS_SAV_RATE:
                    self.__recovery_spending = self._d_k_eff_t
                
                    
                else:
                    self.__recovery_spending = SUBS_SAV_RATE
        return


    def set_k_eff_0(self):
        self._k_eff_0 = (self.__inc_0 - self.__inc_sp)/((1-self.__tax_rate)*PI)

        return

    def _update_income_sp(self, L_t, K):

        self._d_inc_sp_t = ((L_t/K) * self.__inc_sp)/DT_STEP
        
        if self._d_inc_sp_t <0:
            self._d_inc_sp_t=0
        return

    def _update_income(self):

        #if not self.__poverty_trap:
        self._d_inc_t = ((1-self.__tax_rate) * PI * self._d_k_eff_t + self._d_inc_sp_t)/DT_STEP
        if self._d_inc_t <0:
            self._d_inc_t=0
        #else:
            #self._d_inc_t = self.__inc_0 - (self.__inc_sp - self._d_inc_sp_t)
        return

    def _update_consum(self):
        #if not self.__poverty_trap:

        self._update_reco_spend()
        if self.__recovery_type < 2:
            self._d_con_t = self._d_inc_t + self.__recovery_spending
            if self.__recovery_type == 1:
                if self._t <= self.__tf:
                    self.__con_smooth = self.__floor
                else:
                    self.__con_smooth = self._d_con_t
            else:
                self.__con_smooth = self._d_con_t
        elif self.__recovery_type == 2:
            self._d_con_t = self._d_inc_t + self._possible_reco()
            if self._t <= self.__tf:
                self.__con_smooth = self.__floor
            else:
                self.__con_smooth = self._d_con_t
        else:
            self._d_con_t = self._d_inc_t
            if self._t <= self.__tf:
                self.__con_smooth = self.__con_smooth
            else:
                self.__con_smooth = self._d_con_t
    
        if self._d_con_t <0:
            self._d_con_t=0
            self.__recovery_spending = 0.
            self.__recovery_type = 0
        
        if self.__con_smooth <0:
            self.__con_smooth=0
        
        return

    def _update_k_eff(self):
        #if not self.__poverty_trap:

        self._d_k_eff_t -= self.__recovery_spending
        
        if self._d_k_eff_t <=0:
            self._d_k_eff_t=0
            self.__recovery_spending = 0.
            self.__recovery_type = 0
        return

    def _update_wb(self):

        self.__twb += self._dt
        con_0 = self.__con_0/self.__n_inds
        d_con_t = self._d_con_t/self.__n_inds
        self.__wb_0 += ((con_0**(1-ETA))/(1-ETA)) * self._dt * np.e**(-RHO * self.__twb)
        self.__wb_t += (1/(1-ETA)) * (con_0 - d_con_t)**(1-ETA) * self._dt * np.e**(-RHO * self.__twb) 
        self._d_wb_t = (self.__wb_0 - self.__wb_t)#/self.__wb_0

        return
    
    def _update_wb_sav(self):
        con_0 = self.__con_0/self.__n_inds
        d_con_t = self.__con_smooth/self.__n_inds
        self.__wb_0_sm += ((con_0**(1-ETA))/(1-ETA)) * self._dt * np.e**(-RHO * self.__twb)
        self.__wb_t_sm += (1/(1-ETA)) * (con_0 - d_con_t)**(1-ETA) * self._dt * np.e**(-RHO * self.__twb) 
        self.__wb_smooth = (self.__wb_0_sm - self.__wb_t_sm)#/self.__wb_0

        return
    
    def __update_savings(self):

        if (self._d_con_t < 0.05*self.__con_0/DT_STEP) & (self.__sav_t<self.__sav_0):
            self.__sav_t += self.__sav_0/DT_STEP
        elif (self._t <= self.__tf) & (self.__sav_t> 0.0):
            self.__sav_t -= self._d_con_t - self.__floor
        
        if self.__sav_t< 0.0:
            self.__sav_t=0.0
        return
    
    #def __determine_sav_opt(self):
        
        # recovery below subsistence 
        
        

    def __smooth_with_savings_0(self, vul=0.3):
        """Sets the floor taken from savings to smoothen HH's consumption
        loss and the time tf when it runs out of savings
        """

        dc0 = self.__k_eff_0 * vul * (PI+self.__lmbda[self._c_shock])

        if dc0 == 0:
            self.__floor = 0
            self.__tf = T_RNG
            return
        if self.__lmbda[self._c_shock] == 0:
            self.__floor = int(round(min(dc0, max(dc0-(2/3)
                                                * self.sav_t, 0.)), 0))
            self.__tf = 1.
            return

        gamma = dc0
        last_result = None

        while True:

            beta = gamma/dc0
            result = dc0 * (1-beta) + gamma * np.log(beta) - self.__sav_t * self.__lmbda[self._c_shock]
            try:
                if (last_result < 0 and result > 0) or\
                   (last_result > 0 and result < 0):
    
                    _t = -np.log(beta)/self.__lmbda[self._c_shock]
    
                    if _t < 0:
                        print('RESULT!:\ngamma = ', gamma, '& beta = ',
                              beta, ' & t = ', _t)
                        print('CHECK:', dc0 * np.e**(self.__lmbda[self._c_shock] * _t),
                              ' gamma = ', gamma)
        
                    if _t >= T_RNG:
                        self.__floor = int(round(min(dc0, max(dc0-(2/3)
                                                              * self.__sav_t, 0.)), 0.))
                        self.__tf = 1.
                        return

                    self.__floor = int(round(gamma, 0))
                    self.__tf = round(_t, 3)
                    return

            except: pass

            last_result = result
            gamma -= 0.01 * dc0
            if gamma <= 0:
                self.__floor = 0
                self.__tf = T_RNG
                return
        return

    def __smooth_with_savings_1(self):
        """Sets the floor taken from savings to smoothen HH's consumption
        loss and the time tf when it runs out of savings
        """

        dc0 = self._d_con_t
        
        if self.__sav_t <= 0:
            self.__floor = self._d_con_t
            self.__tf = 0.
            return
        
        if self.__sav_t > dc0/(self.__lmbda[self._c_shock]/DT_STEP):
            self.__floor = 0.
            self.__tf = 40.
            return
        
        f = 0.01
        last_result = (1/(self.__lmbda[self._c_shock]/DT_STEP)) * (dc0-f*(np.log(dc0) - np.log(f)+1)) - self.__sav_t
        
        f += 0.1

        while f<dc0:

            result = (1/(self.__lmbda[self._c_shock]/DT_STEP)) *(dc0-f*(np.log(dc0) -np.log(f)+1)) - self.__sav_t

            if (last_result < 0 and result > 0) or\
               (last_result > 0 and result < 0):
                       
                self.__floor = f
                self.__tf = -(1/(self.__lmbda[self._c_shock])) * np.log(f/dc0)
                return
            else:
                last_result = result
                f += 1

        self.__floor = 0.
        self.__tf = 0.

        return
    
    def __smooth_with_savings_2(self):
        """Sets the floor taken from savings to smoothen HH's consumption
        loss and the time tf when it runs out of savings
        """
        


        if self.__sav_t <= 0:
            self.__floor = self._d_con_t
            self.__tf = 0.0
            return

        dc0 = self._d_con_t

        f_1 = dc0 + np.sqrt(dc0**2 - (dc0**2 - 2 * self.__sav_t*self.__recovery_spending*PI/DT_STEP))

        f_2 = dc0 - np.sqrt(dc0**2 - (dc0**2 - 2 * self.__sav_t*self.__recovery_spending*PI/DT_STEP))

        if f_1 <=dc0:
            self.__floor = f_1

        elif f_2 <=dc0:
            self.__floor = f_2

        else:
            raise ValueError

        self.__tf = ((dc0 - self.__floor)/(self.__recovery_spending*PI/DT_STEP)/DT_STEP)

        if self.__floor < 0:
            self.__tf = ((dc0 - 0.)/(self.__recovery_spending*PI/DT_STEP)/DT_STEP)
            self.__floor = 0.

        return
    
    def __smooth_with_savings_3(self):
        """Sets the floor taken from savings to smoothen HH's consumption
        loss and the time tf when it runs out of savings
        """
        
        if self.__sav_t <= 0:
            self.__floor = self._d_con_t
            self.__tf = 0.0
            return

        dc0 = self._d_con_t
        
        f_1 = dc0 + np.sqrt(dc0**2 - (dc0**2 - 2 * self.__sav_t*SUBS_SAV_RATE*PI/DT_STEP))
        
        f_2 = dc0 - np.sqrt(dc0**2 - (dc0**2 - 2 * self.__sav_t*SUBS_SAV_RATE*PI/DT_STEP))
        
        if f_1 <=dc0:
            self.__floor = f_1
        
        elif f_2 <=dc0:
            self.__floor = f_2
        
        else:
            raise ValueError

        self.__tf = ((dc0 - self.__floor)/(SUBS_SAV_RATE*PI/DT_STEP)/DT_STEP)

        if self.__floor < 0:
            self.__tf = ((dc0 - 0.)/(SUBS_SAV_RATE*PI/DT_STEP)/DT_STEP)
            self.__floor = 0.

        return
    
    def __smooth_with_savings_4(self, vul=0.3):
        """Sets the floor taken from savings to smoothen HH's consumption
        loss and the time tf when it runs out of savings
        """

        dc0 = self._d_con_t


        f_1 = dc0 + np.sqrt(dc0**2 - (dc0**2 - 2 * self.__sav_t*SUBS_SAV_RATE*PI/DT_STEP))

        f_2 = dc0 - np.sqrt(dc0**2 - (dc0**2 - 2 * self.__sav_t*SUBS_SAV_RATE*PI/DT_STEP))
        if f_1 <=dc0:
            self.__floor = f_1

        elif f_2 <=dc0:
            self.__floor = f_2

        else:
            raise ValueError

        self.__tf = ((dc0 - self.__floor)/(SUBS_SAV_RATE*PI/DT_STEP)/DT_STEP)

        if self.__floor < 0:
            self.__tf = ((dc0 - 0.)/(SUBS_SAV_RATE*PI/DT_STEP)/DT_STEP)
            self.__floor = 0.

        return


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
        
        if vul == np.nan:
            
            print('invalid k_eff')
            print(self.__hhid)
            self.__k_eff_0 = (self.__inc_0 - self.__inc_sp)/((1-self.__tax_rate)*PI)
            print(self.__k_eff_0)
            vul = 0.3

        last_integ = None
        last_lambda = None

        lmbda = 0.0

        #print('vul='+str(vul))

        c=0

        while True:

            integ = 0.0
            
            # if self.__cnt_ts == 0:
            
            #     for dt in np.linspace(0, T_RNG, DT_STEP*T_RNG):
            #         integ += np.e**(-dt * (RHO + lmbda)) * ((PI + lmbda) * dt - 1) * (PI - (PI + lmbda) * vul * np.e**(-lmbda * dt))**(-ETA)
            # else:
                
            #     for dt in np.linspace(0, T_RNG, DT_STEP*T_RNG)[:-self.__cnt_ts]:
            #         integ += np.e**(-dt * (RHO + lmbda)) * ((PI + lmbda) * dt - 1) * (PI - (PI + lmbda) * vul * np.e**(-lmbda * dt))**(-ETA)
            
            if self.__cnt_ts == 0:
            
                for dt in np.linspace(0, T_RNG, DT_STEP*T_RNG):
                    integ += np.e**(-dt * lmbda) * ((PI + lmbda) * dt - 1) 
            else:
                
                for dt in np.linspace(0, T_RNG, DT_STEP*T_RNG)[:-self.__cnt_ts]:
                    integ += np.e**(-dt * lmbda) * ((PI + lmbda) * dt - 1) 
            #print(integ)
            
            if last_integ and ((last_integ < 0 and integ > 0) or
                                (last_integ > 0 and integ < 0)):
                # print('\n Found the Minimum!\n lambda = ', last_lambda,
                #       '--> integ = ', last_integ)
                # print('lambda = ', lmbda, '--> integ = ', integ, '\n')

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
