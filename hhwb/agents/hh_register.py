#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 13:32:14 2020

@author: insauer
"""
import os
import numpy as np
import pandas as pd
from hhwb.agents.household import Household

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))


class HHRegister():
    """This class represents the intersection between the FIES and the model. All households
       listed in the FIES are initialised as an household object extracting information from the
       relevant columns. All required column names need to be provided during the initialisation
       of the HHRegister. The major attribute needed for further modeling is then a list of all
       registered households.

       Attributes:

           hh_list (list): list with the registered households
           n_hh (int): number of registered households
    """

    def __init__(self):

        self.__hh_list = []
        self.__n_hh = 0

    @property
    def hh_list(self):
        return self.__hh_list

    @property
    def n_hh(self):
        return self.__n_hh

    def set_from_csv(self, path='/data/test_data.csv', id_col='HHID', weight_col='weight',
                     vul_col='vul', income_col='income', income_sp='income_sp'):
        """
        This function reads the household information from the FIES if it is given in a
        csv file. It provides the full list of household agents.

        Parameters:
            path (str): path of the input data frame
            id_col (str): column name of the column containing the household ID
            income_col (str): column name of the column containing the household's total income
            income_sp (str): column name of the column containing the household's income from
                              social transfers
            weight_col (str): column name of the column containing the household's weight
            vul_col (str): column name of the column containing the household's vulnerability
        """
        data = pd.read_csv(DATA_DIR + path)

        self.__n_hh = data.shape[0]
        hh_list = []
        for hhid in range(self.__n_hh):

            hh_data = data.iloc[hhid, :]
            hh_id = hh_data[id_col]
            hh_w = hh_data[weight_col]
            hh_vul = hh_data[vul_col]
            hh_inc = hh_data[income_col]
            hh_inc_sp = hh_data[income_sp]

            hh = Household(hhid=hh_id, w=hh_w, vul=hh_vul, i_0=hh_inc, i_sp=hh_inc_sp)
            hh_list.append(hh)

        self.__hh_list = hh_list

        return
