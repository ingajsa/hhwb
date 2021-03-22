#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 12:41:42 2020

@author: insauer
"""
import os
import random
import numpy as np
import pandas as pd
from hhwb.agents.agent import Agent
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))


AGENT_TYPE = 'SH'

class Shock(Agent):
    """Shock definition. This class builds the intersection with Climada and provides direct
       damage obtained from Climada and the affected households.

        Attributes:
            aff_hh (list): list with affected households
            unaff_hh (list): list with uneffected households
            aff_hh_id (list): list with IDs of affected households
            unaff_hh_id (list): list with IDs of unaffected households
            L (float): total damage
    """

    def __init__(self):

        self.__time_stemps = []
        self.__event_names = []
        self.__aff_ids = np.array([[]])

        self.__L = 0

    @property
    def aff_ids(self):
        return self.__aff_ids
    
    @property
    def time_stemps(self):
        return self.__time_stemps

    @property
    def unaff_hh(self):
        return self.__unaff_hh

    @property
    def L(self):
        return self.__L
    
    @property
    def dt(self):
        return self.__dt
    
    

    def set_random_shock(self, n_events=2, n_hhs=10):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade intersection.

            Parameters
            ----------
            hhs : list
                list with all households
        """

        self.__aff_ids = np.zeros((n_hhs, n_events))
        self.__set_time_stemps(n_events)

        for ev in range(n_events):
            sel_hhs = self.__set_aff_hhs(n_hhs)
            self.__aff_ids[sel_hhs, ev] = 1

        return
    
    def set_shock_from_csv(self, path='../../data/hazard/damage_ts_PHL_wfdei_dbh.csv',
                           sim_period=31, dist=52, hh_reg=None):

        shock_series = pd.read_csv(DATA_DIR + path)

        self.__set_time_stemps_disaster_set(shock_series, dist)
        
        self.__set_aff_hhs_disaster_set(shock_series, hh_reg)

        return

    def __set_time_stemps_disaster_set(self, df, dist):
        
        df['Total'] = df.iloc[:, 2:].sum(axis=1)

        event_years = np.array(df.loc[df['Total'] != 0, 'Year'])

        event_timestemps_wk = (np.array(df.loc[df['Total'] != 0, 'Year'])
                               - event_years[0]) * dist

        add_period_yr = RECO_PERIOD - (event_years[-1]-event_years[0])

        pre_period_wk = int(add_period_yr/2)*dist

        self.__time_stemps = event_timestemps_wk + pre_period_wk

        return
    
    def __set_aff_hhs_disaster_set(self, df, hh_reg):

        region_dict = {'PH150000000': 15,
                       'PH140000000': 14,
                       'PH130000000': 13,
                       'PH010000000': 1,
                       'PH020000000': 2,
                       'PH030000000': 3,
                       'PH040000000': 41,
                       'PH170000000': 42,
                       'PH090000000': 9,
                       'PH050000000': 5,
                       'PH060000000': 6,
                       'PH070000000': 7,
                       'PH080000000': 8,
                       'PH100000000': 10,
                       'PH110000000': 11,
                       'PH120000000': 12,
                       'PH160000000': 16}
        
        self.__get_k_eff_reg(hh_reg)
        

        return

    def __get_k_eff_reg(self, hh_reg):

        for reg in hh_reg.region_hhs:
            
            for hhid in hh_reg.region_hhs[reg]:
                

            #print('hä')


    def __set_aff_hhs(self, n_hhs):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade intersection.
    
            Parameters
            ----------
            hhs : list
                list with all households
        """
        n_aff_hh = 2000
        sel_hhs = []

        while len(sel_hhs) < n_aff_hh:
            hh = random.randint(0, n_hhs-1)
            if not np.isin(hh, sel_hhs):
                sel_hhs.append(hh)
                print('Household appended')
                
        return sel_hhs

        
    def __set_time_stemps(self, n_events=2):
        """The function selects households randomly hocks them. (Function is only a
           placeholder for a real Climade intersection.

            Parameters
            ----------
            n_events : int
                number events
        """

        self.__time_stemps = [50]
        self.__event_names = ['0']
        add = 0
        for run in range(n_events-1):
            self.__time_stemps.append(random.randint(200+add, 300+add))
            self.__event_names.append(str(self.__time_stemps[run]))
            add+=300

    def shock(self, event_index, gov, hhs, dt_reco):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade intersection.

            Parameters
            ----------
            gov : Government
                the government of all households
        """
        affected_hhs = np.where(self.__aff_ids[:, event_index]==1)[0]
        L_t = 0
        for h_ind, hh in enumerate(hhs):
            if h_ind in affected_hhs:
                L_t += hh.vul * (hh.k_eff_0 - hh.d_k_eff_t)
            else:
                L_t += hh.d_k_eff_t

        for h_ind, hh in enumerate(hhs):
            if h_ind in affected_hhs:
                hh.shock(aff_flag=True, L=L_t, K=gov.K, dt=dt_reco)
            else:
                hh.shock(aff_flag=False, L=L_t, K=gov.K, dt=dt_reco)

        gov.shock(aff_flag=True, L=L_t)

        return

    def __plot(self, ax):
        """The function selects households randomly and plots their recovery.

            Parameters
            ----------
            gov : Government
                the government of all households
            hhs : Households
                list with households
        """

        ax.set_xlim((0, RECO_PERIOD*DT_STEP+1))
        ax.set_xlabel('weeks after shock')
        ax.set_ylabel('')
        ax.set_title('Consumption_recovery')
        #ax.legend()

        return
