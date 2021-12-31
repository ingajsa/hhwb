#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 11:34:05 2021

@author: insauer
"""
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd

#survey_data = pd.read_csv('/home/insauer/projects/WB_model/hhwb/data/FIES_prepared/survey_PHL.csv')

# region_dict = {312: 'MW312',
#                305: 'MW305',
#                315: 'MW315',
#                310: 'MW310',
#                304: 'MW304',
#                101: 'MW101',
#                208: 'MW208',
#                204: 'MW204',
#                102: 'MW102',
#                201: 'MW201',
#                #106: 'MW106',
#                206: 'MW206',
#                210: 'MW210',
#                302: 'MW302',
#                301: 'MW301',
#                207: 'MW207',
#                308: 'MW308',
#                306: 'MW306',
#                105: 'MW105',
#                107: 'MW107',
#                313: 'MW313',
#                103: 'MW103',
#                202: 'MW202',
#                311: 'MW311',
#                209: 'MW209',
#                203: 'MW203',
#                309: 'MW309',
#                104: 'MW104',
#                205: 'MW205',
#                307: 'MW307',
#                303: 'MW303',
#                314: 'MW314'}

REGION_DICT = {15: 'PH150000000',
                14: 'PH140000000',
                13: 'PH130000000',
                1: 'PH010000000',
                2: 'PH020000000',
                3: 'PH030000000',
                41: 'PH040000000',
                42: 'PH170000000',
                9: 'PH090000000',
                5: 'PH050000000',
                6: 'PH060000000',
                7: 'PH070000000',
                8: 'PH080000000',
                10: 'PH100000000',
                11: 'PH110000000',
                12: 'PH120000000',
                16: 'PH160000000'}

regions = list(REGION_DICT)
full_df = pd.DataFrame()

for reg in regions:
    
    new_df = pd.read_csv('/home/insauer/projects/WB_model/hhwb/data/surveys_prepared/PHL/region_hh_30as_{}_pack.csv'.format(str(reg)))
    
    full_df = full_df.append(new_df, ignore_index=True)

full_df = full_df.drop(columns=['Unnamed: 0'])
full_df['fhhid'] = np.arange(full_df.shape[0])
print(full_df.shape[0])
full_df.to_csv('/home/insauer/projects/WB_model/hhwb/data/surveys_prepared/PHL/region_hh_full_pack_PHL.csv')