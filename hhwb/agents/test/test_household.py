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

hh_reg = HHRegister()

## create HH agents in the register
hh_reg.set_from_csv(path='/data/FIES_prepared/region_hh_full_pack.csv', id_col='fhhid', weight_col='weight',
                      vul_col='vul', income_col='income', income_sp='income_sp', region='region',
                      decile='decile', savings='savings', poverty_line='poverty_line',
                      ispoor='ispoor', isurban='isurban')

#hh_reg.set_from_csv()

# print('Households registered')
# ## get number of registered households

all_hhs = hh_reg.hh_list
# # ## set up a goverment
gov = Government()
gov.set_tax_rate(all_hhs)

fld = Shock()
fld.set_shock_from_csv(hh_reg=all_hhs, k_eff=gov.K-gov.L_t)
aff_ids = fld.aff_ids
print('Shocks prepared')
# print(fld.aff_ids)
cl = ClimateLife(all_hhs, fld, gov)
cl.start()
# hh_reg.write_output_files(all_hhs,len(all_hhs), cl.dt_life, '/data/output/')
# print('test')