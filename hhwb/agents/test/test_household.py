#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 13:05:40 2020

@author: insauer
"""
import numpy as np
import sys
sys.path.append('/home/insauer/projects/WB_model/hhwb')
from hhwb.agents.household import Household
from hhwb.agents.hh_register import HHRegister
from hhwb.agents.government import Government
from hhwb.agents.shock import Shock

hh_reg = HHRegister()

## create HH agents in the register
hh_reg.set_from_csv(path='/data/test_data.csv', id_col='HHID', weight_col='weight', vul_col='vul',
                 income_col='income', income_sp='income_sp')
## get number of registered households 
hh_reg.n_hh

all_hhs = hh_reg.hh_list
## set up a goverment

## select information
gov = Government()
gov.set_tax_rate(all_hhs)
print('Tax rate: ' + str(gov.tax_rate))
print('Total expenditure on social programs: ' + str(gov.sp_cost))
print('Total national capital stock: ' + str(gov.K))

for hh in all_hhs:
    print('Capital stock of HH ' + str(int(hh.hhid))+': '+str(hh.k_eff_0))
    
print('Lambda HH 1: '+ str(all_hhs[1].lmbda) + ' per year')
print('Tau HH 1: '+ str(all_hhs[1].tau) + ' years')
all_hhs[5].plot_reco_trajec()
all_hhs[1].plot_reco_trajec()

fld = Shock()
fld.affect(all_hhs)
fld.aff_hh

fld.unaff_hh

print('test')



