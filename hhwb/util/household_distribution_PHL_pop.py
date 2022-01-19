#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 14:09:46 2022

@author: insauer
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 20:08:57 2021

@author: insauer
"""

import sys
sys.path.append('/home/insauer/projects/Climada/climada_python')
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry.multipolygon import MultiPolygon
from climada.entity.exposures.base import Exposures
from climada.engine import Impact


survey_data = pd.read_csv('/home/insauer/projects/WB_model/hhwb/data/surveys_prepared/PHL/survey_PHL.csv')

shape_path = '/home/insauer/projects/WB_model/hhwb/data/shapes/phl_adm_psa_namria_20200529_shp/phl_admbnda_adm1_psa_namria_20200529.shp'

prvS = gpd.GeoDataFrame()
prvS = gpd.read_file(shape_path)

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

region_dict = { 15:'PH150000000',
                14:'PH140000000',
                13:'PH130000000',
                1:'PH010000000',
                2:'PH020000000',
                3:'PH030000000',
                41:'PH040000000',
                42:'PH170000000',
                9:'PH090000000',
                5:'PH050000000',
                6:'PH060000000',
                7:'PH070000000',
                8:'PH080000000',
                10:'PH100000000',
                11:'PH110000000',
                12:'PH120000000',
                16:'PH160000000'}

regions = list(region_dict)
hh_pack_size = 1

hh_inst = 0
#Select only the first region
c_all=0
for reg in regions:
    print('Region' + str(reg))
    
    reg_imp_df = pd.read_csv('/home/insauer/projects/WB_model/hhwb/data/hazard_data/PHL_sat/haz_dat_30as_{}_vul.csv'.format(str(reg)))
    
    reg_data = survey_data.loc[survey_data['region'] == reg]
    
    reg_data = reg_data.reset_index().drop(columns=['index'])
    

    reg_imp_df=reg_imp_df.drop(columns=reg_imp_df.filter(regex='Unnamed').columns)

    #tot_hh = reg_data['weight'].sum().round()
    
    tot_pop = reg_data['weight'].sum()#.round()
    
    shape = prvS.loc[prvS['ADM1_PCODE']==region_dict[reg], 'geometry'].values[0]
    
    
    gdpa = Exposures()
    
    try:
        gdpa.set_from_raster(file_name='/home/insauer/mnt/bene_data/population/columbia_population_resampled/gpw_population_2016_30_sec.tif', band=1, 
                        geometry=shape)
    except TypeError:
        shape = MultiPolygon([shape])
        gdpa.set_from_raster(file_name='/home/insauer/mnt/bene_data/population/columbia_population_resampled/gpw_population_2016_30_sec.tif', band=1,
                        geometry=shape)
    gdpa = gdpa.gdf
    
    gdpa['centroid_id'] = np.arange(gdpa.shape[0])
    gdpa['people_share'] = gdpa['value']/gdpa['value'].sum()

    #gdpa['nhh'] = gdpa['people_share'] * reg_data['weight'].sum()
    gdpa['pop'] = gdpa['people_share'] * reg_data['weight'].sum()
    
    #gdpa['res'] = gdpa['nhh'] - gdpa['nhh'].astype(int)
    #gdpa['nhh'] = gdpa['nhh'].astype(int)
    
    #total_hh = gdpa['nhh'].sum()
    #rest_hh = int(tot_hh - total_hh)
    
    #gdpa = gdpa.sort_values(['res'], ascending=False)
    
    #gdpa.iloc[:rest_hh, 5] += 1
    
    #gdpa = gdpa.sort_values(['nhh'], ascending=True)
    
    # region_df = region_df.sort_values(['hhid', 'hh_instance'])

    reg_imp_df['n_events'] = reg_imp_df.iloc[:, :44].sum(axis=1)
    
    cells_aff = reg_imp_df.loc[reg_imp_df['n_events']>0.000000000]
    cells_unaff = reg_imp_df.loc[reg_imp_df['n_events']==0]

    reg_imp_df.to_csv('/home/insauer/projects/WB_model/hhwb/data/hazard_data/PHL_sat/haz_dat_30as_{}_vul.csv'.format(str(reg)))
    aff_ids = list(cells_aff)
    
    aff_ids = list(cells_aff['Centroid_ID'])
    
    aff_exp = gdpa.loc[gdpa['centroid_id'].isin(aff_ids)]
    
    aff_exp = aff_exp.sort_values(['pop'], ascending=False)
    
    n_aff_hh = aff_exp['pop'].sum()
    print(n_aff_hh)
    
    reg_data['centroid_id'] = -1
    reg_data['n_copied'] = 0
    reg_data['hh_instance'] = 0
    
    aff_hh_data = pd.DataFrame()
    
    index_counter = 0
    print('entries:')
    print(reg_data.shape[0])
    index_max = reg_data.shape[0]-1
    
    c_aff_hhs=0
    
    for centr in aff_ids:
        
        
        n_pop_cntr = aff_exp.loc[aff_exp['centroid_id']==centr,'pop'].sum()
        c_hh = 0
        hh_pack_size = 1
        while c_hh < n_pop_cntr:
            if n_pop_cntr-c_hh <hh_pack_size:
                hh_pack_size = n_pop_cntr-c_hh
                
            if reg_data.iloc[index_counter,:].loc['weight']<hh_pack_size:
                if index_counter == index_max:
                    index_counter = 0
                else:
                    index_counter += hh_pack_size
                continue
            
            aff_hh_data = aff_hh_data.append(reg_data.iloc[index_counter,:], ignore_index=True)
            aff_hh_data.iloc[c_aff_hhs,:].loc['centroid_id'] = centr
            aff_hh_data.iloc[c_aff_hhs,:].loc['weight'] =  hh_pack_size*reg_data.loc[:,'n_individuals'].iloc[index_counter]
            aff_hh_data.iloc[c_aff_hhs,:].loc['hh_instance'] = reg_data.iloc[index_counter,:].loc['n_copied'].sum()+1
            reg_data.loc[:,'n_copied'].iloc[index_counter] += 1
            reg_data.loc[:,'weight'].iloc[index_counter] -= hh_pack_size*reg_data.loc[:,'n_individuals'].iloc[index_counter]
            
            c_hh +=hh_pack_size *reg_data.loc[:,'n_individuals'].iloc[index_counter]
            if index_counter == index_max:
                index_counter = 0
            else:
                index_counter += 1
            c_aff_hhs += 1
            #print(c_aff_hhs)
    
    nhh=aff_hh_data['weight'].sum()
    reg_data = reg_data.append(aff_hh_data, ignore_index= True)
    reg_data = reg_data.reset_index().drop(columns=['index'])
    reg_data.to_csv('/home/insauer/projects/WB_model/hhwb/data/surveys_prepared/PHL/region_hh_30as_{}_pack_pop.csv'.format(str(reg)))
    c_all += reg_data.shape[0]
    print(reg_data.shape[0])
    print(c_all)

print('Total households:')
print(c_all)