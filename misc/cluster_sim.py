#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 16:53:44 2022

@author: insauer
"""
import sys
import os
import inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import argparse

from hhwb.agents.household import Household
from hhwb.agents.hh_register import HHRegister
from hhwb.agents.government import Government
from hhwb.agents.shock import Shock
from hhwb.application.climate_life import ClimateLife
from hhwb.application.data_analysis import DataAnalysis
import psutil

parser = argparse.ArgumentParser(
    description='run hhwb for different shock series')
parser.add_argument(
    '--run_name', type=str, default='shocks',
    help='runoff model')

parser.add_argument(
    '--run_time', type=int, default=160,
    help='runoff model')

args = parser.parse_args()

CLUSTER=True

if CLUSTER==True:
    work_path='/p/projects/ebm/inga/hhwb'
    cores=psutil.cpu_count(logical = True)
else:
    work_path='/home/insauer/projects/WB_model/hhwb'
    cores=7
    
print('Number threads = ' + str(cores))
#hh_reg = HHRegister()

if args.run_name == 'shocks_syn':
    hh_path = '/data/survey_data/PHL/region_hh_full_pack_PHL_pop_syn.csv'
else:
    hh_path = '/data/survey_data/PHL/region_hh_full_pack_PHL_pop.csv'

# create HH agents in the register
# hh_reg.set_from_csv(work_path=work_path, path=hh_path, id_col='fhhid', n_ind = 'n_individuals', weight_col='weight',
#                       vul_col='vul', income_col='income', income_sp='income_sp', region='region',
#                       decile='decile', savings='savings', subsistence_line='subsistence_line',
#                       ispoor='ispoor', isurban='isurban')

# # # print('Households registered')
# # # ## get number of registered households

# all_hhs = hh_reg.hh_list
# # # ## set up a goverment
# gov = Government()
# gov.set_tax_rate(all_hhs)

# fld = Shock()
# #fld.set_shock_from_csv()
# fld.read_shock(work_path=work_path, path='/data/shock_data/'+args.run_name+'.csv',
#                event_identifier='-', run=args.run_name)

# # print('Shocks prepared')
# # # print(fld.aff_ids)
# cl = ClimateLife(all_hhs, fld, gov)
# # cl.start(work_path=work_path, result_path='/data/output_'+args.run_name+'/',
# #          cores=cores, reco_period=args.run_time)
# cl.start(work_path='', result_path='',
#           cores=cores, reco_period=args.run_time)

survey_data_path=work_path+ '/data/survey_data/PHL/survey_'+args.run_name+'analysed.csv'#hh_path
shock_data_path=work_path+'/data/shock_data/'+args.run_name+'.csv'
output_data_path='/p/projects/ebm/inga/hhwb/data/output_'+args.run_name+'/'

da=DataAnalysis(survey_data_path, shock_data_path, output_data_path, column_id='', run_name=args.run_name)

da.analyse_time(step=20000)
da.analyse_wb(step=20000)