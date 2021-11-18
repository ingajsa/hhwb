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
        self.__region_hhs = {}

    @property
    def hh_list(self):
        return self.__hh_list

    @property
    def n_hh(self):
        return self.__n_hh

    @property
    def region_hhs(self):
        return self.__region_hhs

    def set_from_csv(self, path='/data/test_data.csv', id_col='HHID', n_ind='n_individuals',
                     weight_col='weight', vul_col='vul', income_col='income',
                     income_sp='income_sp', region='region', decile='decile', savings='savings',
                     subsistence_line='subsistence_line', ispoor=None, isurban=None):
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
            region (str): column name of the column containing the region in which the household
                          is located
            decile (str): column name of the income decile to which the household belongs
            savings (str): column name of the household's savings'
            subsistence_line (str): name of the column containing the poverty line in the area of the
                         household
            ispoor (str): name of the column that indicates whether the household is poor
            isurban (str): name of the column that indicates whether the household is located in
                           an urban or rural area
        """

        data = pd.read_csv(DATA_DIR + path)
        self.__n_hh = data.shape[0]
        self.__extract_meta_info(data, id_col)

        hh_list = []
        for hh, hhid in enumerate(data[id_col]):

            hh_data = data.iloc[hh, :]
            hh_id = hh_data[id_col]
            hh_w = hh_data[weight_col]
            hh_vul = hh_data[vul_col]
            hh_inc = hh_data[income_col]
            hh_inc_sp = hh_data[income_sp]

            hh_reg = hh_data[region]
            hh_dec = hh_data[decile]
            hh_sav = hh_data[savings]
            hh_inds = hh_data[n_ind]
            hh_pov_line = hh_data[subsistence_line]
            if ispoor:
                hh_poor = hh_data[ispoor]
            else:
                hh_poor = None

            if isurban:
                hh_urban = hh_data[isurban]
            else:
                hh_urban = None

            hh = Household(hhid=hh_id, n_inds=hh_inds, w=hh_w, vul=hh_vul, i_0=hh_inc, i_sp=hh_inc_sp,
                           region=hh_reg, savings=hh_sav, subsistence_line=hh_pov_line, decile=hh_dec,
                           isurban=hh_urban, ispoor=hh_poor)
            print(hhid)
            hh_list.append(hh)

        self.__hh_list = hh_list

        return

    def __extract_meta_info(self, data, id_col):

        regions = list(set(data['region']))

        for region in regions:

            hhids = list(data.loc[data['region'] == region, id_col])

            self.__region_hhs.update({region: hhids})

        return
