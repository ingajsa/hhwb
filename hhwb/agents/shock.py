#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 12:41:42 2020

@author: insauer
"""
import random
import numpy as np
from hhwb.agents.agent import Agent

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

        Agent.__init__(self, AGENT_TYPE)

        self.__aff_hh = []
        self.__unaff_hh = []
        self.__aff_hh_id = []
        self.__unaff_hh_id = []
        self.__L = None

    @property
    def aff_hh(self):
        return self.__aff_hh

    @property
    def unaff_hh(self):
        return self.__unaff_hh

    @property
    def L(self):
        return self.__L

    def set_random_shock(self, hhs):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade intersection.

            Parameters
            ----------
            hhs : list
                list with all households
        """
        n_aff_hh = 4

        while len(self.__aff_hh_id) < n_aff_hh:
            hh = random.randint(0, 9)
            if not np.isin(hh, self.__aff_hh_id):
                self.__aff_hh_id.append(hh)
                self.__aff_hh.append(hhs[hh])

        self.__unaff_hh_id = list(np.arange(len(hhs)))
        self.__unaff_hh_id = list(set(self.__unaff_hh_id).difference(set(self.__aff_hh_id)))
        self.__unaff_hh = [hhs[i] for i in self.__unaff_hh_id]

        return

    def shock(self, gov):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade intersection.

            Parameters
            ----------
            gov : Government
                the government of all households
        """

        for ah in self.__aff_hh:
            ah.shock(aff_flag=True)
            self.__L += ah.vul * ah.k_eff_0

        for uh in self.__unaff_hh:
            uh.shock(aff_flag=False)
        gov.shock()

        return
