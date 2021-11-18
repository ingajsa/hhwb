#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 10:06:38 2021

@author: insauer
"""
import sys
sys.path.append('/home/insauer/Climada/climada_python')
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
from climada.hazard.river_flood import RiverFlood
from climada.hazard.centroids import Centroids
from shapely.geometry.multipolygon import MultiPolygon
from climada.entity.exposures.gdp_asset import GDP2Asset
from climada.entity.exposures.base import Exposures
from climada.entity.impact_funcs.river_flood import flood_imp_func_set
from climada.engine import Impact

PHL_FRC = '/home/insauer/mnt/cama-flood/results/isimip2a/gev/fldfrc_150arcsec_clm40_princeton_flopros.nc'
PHL_DPH  = '/home/insauer/mnt/cama-flood/results/isimip2a/gev/flddph_150arcsec_clm40_princeton_flopros.nc'
gdp_path = '/home/insauer/mnt/ebm/inga/climada_exposures/asset/'

survey_data = pd.read_csv('/home/insauer/projects/WB_model/hhwb/data/surveys_prepared/MWI/final_survey_MWI.csv')

shape_path = '/home/insauer/projects/WB_model/hhwb/data/shapes/mwi_adm_nso_20181016_shp/mwi_admbnda_adm2_nso_20181016.shp'

prvS = gpd.GeoDataFrame()
prvS = gpd.read_file(shape_path)

region_dict = {312: 'MW312',
               305: 'MW305',
               315: 'MW315',
               310: 'MW310',
               304: 'MW304',
               101: 'MW101',
               208: 'MW208',
               204: 'MW204',
               102: 'MW102',
               201: 'MW201',
               #106: 'MW106',
               206: 'MW206',
               210: 'MW210',
               302: 'MW302',
               301: 'MW301',
               207: 'MW207',
               308: 'MW308',
               306: 'MW306',
               105: 'MW105',
               107: 'MW107',
               313: 'MW313',
               103: 'MW103',
               202: 'MW202',
               311: 'MW311',
               209: 'MW209',
               203: 'MW203',
               309: 'MW309',
               104: 'MW104',
               205: 'MW205',
               307: 'MW307',
               303: 'MW303',
               314: 'MW314'}

# region_dict = {15:'PH150000000',
#                14:'PH140000000',
#                13:'PH130000000',
#                1:'PH010000000',
#                2:'PH020000000',
#                3:'PH030000000',
#                41:'PH040000000',
#                42:'PH170000000',
#                9:'PH090000000',
#                5:'PH050000000',
#                6:'PH060000000',
#                7:'PH070000000',
#                8:'PH080000000',
#                10:'PH100000000',
#                11:'PH110000000',
#                12:'PH120000000',
#                16:'PH160000000'}

regions = list(region_dict)
hh_pack_size = 20

hh_inst = 0
#Select only the first region

for reg in regions[17:]:
    print('Region' + str(reg))
    reg_data = survey_data.loc[survey_data['region'] == reg]


    tot_hh = reg_data['weight'].sum().round()

    hh_pack = (reg_data['weight']/hh_pack_size).astype(int)
    n_hh_pack = ((reg_data['weight']/hh_pack_size).astype(int)*hh_pack_size).sum()

    reg_data['res_weight'] = reg_data['weight'] - (reg_data['weight']/hh_pack_size).astype(int)\
    * hh_pack_size

    reg_data['res_weight_dec'] = reg_data['res_weight'] - reg_data['res_weight'].astype(int)

    reg_data['new_weight'] = (reg_data['weight']/hh_pack_size).astype(int)

    c_weight = ((reg_data['weight']/hh_pack_size).astype(int) * hh_pack_size).sum() + reg_data['res_weight'].astype(int).sum()
    res_dec_boarder = tot_hh - c_weight

    reg_data = reg_data.sort_values(['res_weight_dec'], ascending=False)

    region_df = pd.DataFrame()
    c = 0
    hhids = list(reg_data['hhid'])
    print(len(hhids))
    print('packaging')
    for i, hhid in enumerate(hhids):

        #print(str(i) + '   '+str(hhid))
        cur_hh = reg_data.loc[reg_data['hhid'] == hhid]
        w = cur_hh['new_weight'].sum()
        res_w = cur_hh['res_weight'].sum()
        if c < res_dec_boarder:
            cur_hh_1 = cur_hh
            cur_hh_1.loc[:, 'weight'] = hh_pack_size
            region_df = region_df.append([cur_hh_1] * w, ignore_index=True)
            cur_hh_1.loc[:, 'weight'] = int(res_w)
            if int(res_w) > 0:
                hh_instance = np.arange(w+2)
                region_df = region_df.append([cur_hh_1], ignore_index=True)
            else:
                hh_instance = np.arange(w+1)
            cur_hh_1.loc[:, 'weight'] = 1
            region_df = region_df.append([cur_hh_1], ignore_index=True)
            region_df.loc[region_df['hhid'] == hhid, 'hh_instance'] = hh_instance
        else:
            cur_hh_1 = cur_hh
            cur_hh_1.loc[:, 'weight'] = hh_pack_size
            region_df = region_df.append([cur_hh_1]*w, ignore_index=True)
            cur_hh_1.loc[:, 'weight'] = int(res_w)
            if int(res_w) > 0:
                hh_instance = np.arange(w+1)
                region_df = region_df.append([cur_hh_1], ignore_index=True)
            else:
                hh_instance = np.arange(w)
            region_df.loc[region_df['hhid'] == hhid, 'hh_instance'] = hh_instance
        c += 1

    region_df = region_df.sort_values(['hh_instance'])
    region_df = region_df.drop(columns=['Unnamed: 0'])

    
    shape = prvS.loc[prvS['ADM2_PCODE']==region_dict[reg], 'geometry'].values[0]
    

    gdpa = GDP2Asset()
    
    try:
        gdpa.set_from_raster(file_name='/home/insauer/data/Tobias/hyde_ssp2_1980-2015_0150as_yearly_zip.nc4', band=16, src_crs=None, window=False,
                        geometry=shape)
    except TypeError:
        shape = MultiPolygon([shape])
        gdpa.set_from_raster(file_name='/home/insauer/data/Tobias/hyde_ssp2_1980-2015_0150as_yearly_zip.nc4', band=16, src_crs=None, window=False,
                        geometry=shape)

    gdpa['centroid_id'] = np.arange(gdpa.shape[0])
    gdpa['people_share'] = gdpa['value']/gdpa['value'].sum()

    gdpa['nhh'] = gdpa['people_share'] * region_df['weight'].sum()
    gdpa['res'] = gdpa['nhh'] - gdpa['nhh'].astype(int)
    gdpa['nhh'] = gdpa['nhh'].astype(int)
    
    total_hh = gdpa['nhh'].sum()
    rest_hh = int(tot_hh - total_hh)
    
    gdpa = gdpa.sort_values(['res'], ascending=False)
    
    gdpa.iloc[:rest_hh, 5] += 1
    
    gdpa = gdpa.sort_values(['nhh'], ascending=True)
    
    aff_gdpa = gdpa.loc[gdpa['nhh']>0]
    aff_gdpa = aff_gdpa.sort_values(['nhh'], ascending=False)
    
    region_df['Longitude'] = np.nan
    region_df['Latitude'] = np.nan
    region_df['Centroid_ID'] = np.nan
    
    region_df_copy = region_df
    
    count=0
    end = 0
    unfitting_centroids = []
    print('household distribution')
    for centr in aff_gdpa['centroid_id']:
        if centr == 18585:
            print('stop')
        nhh = aff_gdpa.loc[aff_gdpa['centroid_id'] == centr, 'nhh'].values[0]
        rest_hh = nhh
        while rest_hh > 0:
            possible_hh = region_df_copy.loc[region_df_copy['weight'] <= rest_hh]
            if possible_hh.shape[0] == 0:
                possible_hh = region_df_copy.loc[region_df_copy['weight'] == region_df_copy['weight'].min()]
                if region_df_copy.shape[0] == 0:
                    print('missing')
                    break
                print('unfitting')
                unfitting_centroids.append(centr)

            possible_hh = possible_hh.sort_values(['weight'], ascending=False)
            rest_hh = rest_hh-possible_hh['weight'].iloc[0]
            region_df_copy = region_df_copy.drop(labels=possible_hh.iloc[0].name)
            region_df['Longitude'].iloc[possible_hh.iloc[0].name] = gdpa.loc[gdpa['centroid_id'] == centr,'longitude'].values[0]
            region_df['Latitude'].iloc[possible_hh.iloc[0].name] = gdpa.loc[gdpa['centroid_id'] == centr,'latitude'].values[0]
            region_df['Centroid_ID'].iloc[possible_hh.iloc[0].name] = gdpa.loc[gdpa['centroid_id'] == centr,'centroid_id'].values[0]

    region_df = region_df.sort_values(['hhid', 'hh_instance'])

    data = pd.read_csv('/home/insauer/projects/NC_Submission/Data/check_data/assembled_data_regions.csv')

    data_phl = data.loc[(data['Country'] == 'MWI') & (data['Year'] > 1979) & (data['clim_forc'] != 'watch')]
    natcat = data_phl.loc[((data_phl['clim_forc'] == 'princeton') & (data_phl['GHM'] == 'dbh')),
                          ['Year', 'natcat_flood_damages_2005_CPI']]

    years = np.arange(1980, 2011)
    # calculate country level damage
    rf_reg = RiverFlood()
    rf_reg.set_from_nc(shape=shape, years=years, dph_path=PHL_DPH, frc_path=PHL_FRC)

    if_set = flood_imp_func_set()

    gdpa = GDP2Asset()

    gdpa.set_from_raster(file_name='/home/insauer/data/Tobias/gdp_1850_2100_150arcsec.nc', band=161
                         , src_crs=None, window=False, geometry=shape)

    gdpa['centroid_id'] = np.arange(gdpa.shape[0])

    gdpa['if_RF'] = 7
    gdpa['value'] = 1
    gdpa.value_unit = 'USD'

    imp_0 = Impact()
    imp_0.calc(gdpa, if_set, rf_reg, save_mat=True)

    rf_cnt = RiverFlood()
    rf_cnt.set_from_nc(countries=['MWI'], years=years, dph_path=PHL_DPH, frc_path=PHL_FRC)

    gdpa_1 = Exposures()
    gdpa_1.read_hdf5(gdp_path + 'asset_MWI_{}.h5'.format(str(2010)))
    imp_cnt = Impact()
    imp_cnt.calc(gdpa_1, if_set, rf_cnt, save_mat=True)

    imp_df = np.array(imp_0.imp_mat.todense()).T
    c_fac = np.array(natcat['natcat_flood_damages_2005_CPI'])/imp_cnt.at_event

    fac_bool = np.where(c_fac > 0)
    c_fac[fac_bool[0]] = 1
    c_fac = np.reshape(c_fac, (31))

    imp_df_vul = imp_df*c_fac

    column_names = np.arange(1980, 2011).astype(str)
    reg_imp_df = pd.DataFrame(data=imp_df_vul, columns=column_names)
    reg_imp_df['Longitude'] = imp_0.coord_exp[:, 1]
    reg_imp_df['Latitude'] = imp_0.coord_exp[:, 0]
    reg_imp_df['Centroid_ID'] = gdpa['centroid_id']

    reg_imp_df = reg_imp_df.drop(index=reg_imp_df[reg_imp_df.iloc[:, :31].sum(axis=1) == 0].index)

    non_affected_cells = list(set(gdpa['centroid_id']) ^ set(reg_imp_df['Centroid_ID']))
    affected_cells = list(set(gdpa['centroid_id']) & set(reg_imp_df['Centroid_ID']))

    reg_imp_df.to_csv('/home/insauer/projects/WB_model/hhwb/data/hazard_data/MWI/haz_data_reg{}.csv'.format(str(reg)))

    hhids = list(set(region_df.loc[region_df['Centroid_ID'].isin(non_affected_cells), 'hhid']))
    print('household aggregation')
    for hhid in hhids:
        region_df.loc[region_df['Centroid_ID'].isin(non_affected_cells) &
                      (region_df['hhid'] == hhid),
                      'weight'] = region_df.loc[region_df['Centroid_ID'].isin(non_affected_cells)
                                                & (region_df['hhid'] == hhid), 'weight'].sum()
        region_df = region_df.drop(index=region_df.loc[region_df['Centroid_ID'].isin(non_affected_cells)
                                                       & (region_df['hhid'] == hhid)].index[1:])
    print(region_df.shape[0])
    
    hh_inst += region_df.shape[0]
    
    region_df.to_csv('/home/insauer/projects/WB_model/hhwb/data/surveys_prepared/MWI/region_hh_{}_pack.csv'.format(str(reg)))

print('Total households:')
print(hh_inst)