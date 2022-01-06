#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:05:40 2020

@author: insauer
"""
import numpy as np
from hhwb.agents.household import Household
from hhwb.agents.hh_register import HHRegister
from hhwb.agents.government import Government
from hhwb.agents.shock import Shock
from hhwb.application.climate_life import ClimateLife
import psutil

CLUSTER=False

if CLUSTER==True:
    work_path='/p/projects/ebm/inga/hhwb'
    cores=16
else:
    work_path='/home/insauer/projects/WB_model/hhwb'
    cores=7
    
print(psutil.cpu_count(logical = True))
hh_reg = HHRegister()

# create HH agents in the register
hh_reg.set_from_csv(work_path=work_path, path='/data/survey_data/PHL/region_hh_full_pack_PHL.csv', id_col='fhhid', n_ind = 'n_individuals', weight_col='weight',
                      vul_col='vul', income_col='income', income_sp='income_sp', region='region',
                      decile='decile', savings='savings', subsistence_line='subsistence_line',
                      ispoor='ispoor', isurban='isurban')

# print('Households registered')
# ## get number of registered households

all_hhs = hh_reg.hh_list
# # ## set up a goverment
gov = Government()
gov.set_tax_rate(all_hhs)

fld = Shock()
#fld.set_shock_from_csv()
fld.read_shock(work_path=work_path, path='/data/output/shocks/shocks.csv', event_identifier='-', run='')

# print('Shocks prepared')
# # print(fld.aff_ids)
cl = ClimateLife(all_hhs, fld, gov)
cl.start(work_path=work_path, result_path='/data/output/', cores=cores)
# hh_reg.write_output_files(all_hhs,len(all_hhs), cl.dt_life, '/data/output/')
# print('test')