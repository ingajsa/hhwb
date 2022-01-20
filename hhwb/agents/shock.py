#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 12:41:42 2020

@author: insauer
"""
import os
import random
import numpy as np
import pandas as pd
import multiprocessing as mp
from hhwb.agents.agent import Agent
from functools import partial
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD
from datetime import datetime, date


AGENT_TYPE = 'SH'

REGION_DICT = { 
                15: 'PH150000000',
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

# REGION_DICT = {312: 'MW312',
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

REGIONS = list(REGION_DICT)



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

        self.__time_stemps = []
        self.__event_names = []
        self.__aff_ids = np.array([[]])

        self.__L = 0

    @property
    def aff_ids(self):
        return self.__aff_ids

    @property
    def time_stemps(self):
        return self.__time_stemps

    @property
    def unaff_hh(self):
        return self.__unaff_hh

    @property
    def L(self):
        return self.__L

    @property
    def dt(self):
        return self.__dt
    
    @staticmethod
    def apply_shock(hh, L, K, dt_reco, affected_hhs):

        if hh.hhid in affected_hhs:
            
            hh.shock(aff_flag=True, L=L, K=K, dt=dt_reco)
        else:
            hh.shock(aff_flag=False, L=L, K=K, dt=dt_reco)
        return hh


    def set_random_shock(self, n_events=3, n_hhs=10):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade intersection.

            Parameters
            ----------
            hhs : list
                list with all households
        """

        self.__aff_ids = np.zeros((n_hhs, n_events))
        self.__set_time_stemps(n_events)

        for ev in range(n_events):
            sel_hhs = self.__set_aff_hhs(n_hhs)
            self.__aff_ids[sel_hhs, ev] = 1

        return
    
    def set_shock_from_csv(self, work_path='/home/insauer/projects/WB_model/hhwb',
                           path_haz='/data/hazard_data/PHL_sat/haz_dat_30as_{}_vul.csv',
                           path_hh='/data/surveys_prepared/PHL/region_hh_full_pack_PHL.csv',
                           hh_reg=None, k_eff=0):

        

        df_hh = pd.read_csv(work_path + path_hh)
        
        
        for r, reg in enumerate(REGIONS):

            event_names, weeks = self.__set_time_stemps_disaster_set(work_path, path_haz, reg, '-')
            
            if r==0:
                self.__aff_ids = np.zeros((len(df_hh), len(event_names)))
                self.__time_stemps = weeks
                

            self.__set_aff_hhs_disaster_set(work_path, path_haz, path_hh, reg, event_names)
        
        
        shock_df = pd.DataFrame(data=self.__aff_ids, columns=event_names)
        
        shock_df['region']=df_hh['region']
        
        shock_df.to_csv(work_path + '/data/output/shocks/shocks_99.csv')


        return
    
    def generate_single_shocks(self, work_path='/home/insauer/projects/WB_model/hhwb',
                           path_haz='/data/output/shocks/shocks_99_aggregated.csv',
                           path_hh='/data/surveys_prepared/PHL/region_hh_full_pack_PHL_pop.csv',
                           path_hh_orig='/data/surveys_prepared/PHL/survey_PHL.csv',
                           hh_reg=None, k_eff=0):


        df_hh = pd.read_csv(work_path + path_hh)
        df_hh_orig = pd.read_csv(work_path + path_hh_orig)
        df_shock = pd.read_csv(work_path + path_haz)
        
        self.__time_stemps = df_shock.columns[1:].astype(int)
        
        event_names = df_shock.columns[1:]
        
        df_shock['region']=df_hh['region']
        
        df_shock['fhhid']=df_shock.index
        
        new_survey_data=pd.DataFrame()
        
        for r, reg in enumerate(REGIONS):
            
            print(reg)
            
            df_hh_reg = df_hh.loc[df_hh['region']==reg]
            df_hh_orig_reg = df_hh_orig.loc[df_hh_orig['region']==reg]
            df_shock_reg = df_shock.loc[df_shock['region']==reg]
            
            aff_hh_data = self.__shock_hh_without_location(df_hh_reg, df_hh_orig_reg, df_shock_reg, reg)
            
            new_survey_data=new_survey_data.append(aff_hh_data, ignore_index=True)
        
        self.__aff_ids = np.zeros((len(new_survey_data), len(event_names)))
        
        for i, sh in enumerate(self.__time_stemps):
            
            self.__aff_ids[:,i]=(new_survey_data.loc[:, 'event']==sh).astype(int)
        
        shocks = pd.DataFrame(data=self.__aff_ids, columns=event_names)
        
        shocks['region']=new_survey_data['region']
        
        new_survey_data.to_csv(work_path + '/data/surveys_prepared/PHL/region_hh_full_pack_PHL_pop_syn.csv')
        
        shocks.to_csv('/home/insauer/projects/WB_model/hhwb/data/shocks/shocks_synthetic.csv')
        
        return


    def __shock_hh_without_location(self, df_hh_reg, df_hh_orig_reg, df_shock_reg, reg):
        
        c_aff_hhs = 0
        
        hhids=list(df_hh_orig_reg['hhid'])
        
        df_hh_orig_reg['n_copied']=0
        df_hh_orig_reg['event']=0
        df_hh_orig_reg['hh_instance']=0
        
        aff_hh_data = pd.DataFrame()
        
        for shock in self.__time_stemps:
            print(shock)
            aff_ids_shock = list(df_shock_reg.loc[df_shock_reg[str(shock)]==1,'fhhid'])
            n_aff_shock = df_hh_reg.loc[df_hh_reg['fhhid'].isin(aff_ids_shock), 'n_individuals'].sum()
            
            c_hh = 0
            print(n_aff_shock)
            
            while c_hh < n_aff_shock:
                
                hh_ind = random.choice(hhids)
                
                if df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind,'weight'].sum() <= \
                    df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind,'n_individuals'].sum():
                    hhids.remove(hh_ind)
                    print('stop')
                    print(len(hhids))
                    continue
                
                aff_hh_data = aff_hh_data.append(df_hh_orig_reg.loc[df_hh_orig_reg['hhid']==hh_ind], ignore_index=True)
                aff_hh_data.loc[c_aff_hhs,'weight'] =  df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind,'n_individuals'].sum()
                aff_hh_data.loc[c_aff_hhs,'hh_instance'] = df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind,'n_copied'].sum()+1
                aff_hh_data.loc[c_aff_hhs, 'event'] = shock
                df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind, 'n_copied'] += 1
                df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind, 'weight'] -= df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind,'n_individuals'].sum()
                
                c_hh +=df_hh_orig_reg.loc[df_hh_orig_reg['hhid']== hh_ind,'n_individuals'].sum()

                c_aff_hhs += 1
                
            aff_hh_data['region']=reg
        aff_hh_data = df_hh_orig_reg.append(aff_hh_data, ignore_index=True)
        
        return aff_hh_data

    def __set_time_stemps_disaster_set(self, work_path, path_haz, reg, event_identifier):
        
        start_date = date(2000,1,1)
        
        shock_series = pd.read_csv((work_path + path_haz).format(reg))
            
        event_names = [col for col in shock_series.columns if '-' in col]
        
        dates_list = [datetime.strptime(event, '%Y-%m-%d').date() for event in event_names]
        
        weeks = [np.round((event_date - start_date).days/7) for event_date in dates_list]
        
        return event_names, weeks
    
    def read_shock(self, work_path, path, event_identifier,run):
        
        shock_data = pd.read_csv(work_path + path)
        
        start_date = date(2002,1,1)
        
        event_names = [col for col in shock_data.columns if event_identifier in col]
        
        event_data=shock_data[event_names]
        
        dates_list = [datetime.strptime(event, '%Y-%m-%d').date() for event in event_names]
        
        month = [np.round((event_date - start_date).days/28) for event_date in dates_list]
        
        self.__aff_ids=np.zeros((event_data.shape[0],len(list(set(month)))))
        
        week_nums=list(set(month))
        
        week_nums.sort()
        for m, week in enumerate(week_nums):
            
            locat = np.where(np.array(month)==week)
            
            self.__aff_ids[:, m]=event_data.iloc[:, locat[0]].sum(axis=1).clip(upper=1)
        
        self.__time_stemps = week_nums
        
        shock_df = pd.DataFrame(data=self.__aff_ids, columns=np.array(week_nums).astype(int).astype(str))
        
        #shock_df.to_csv(work_path + '/data/output_{}/shocks_aggregated.csv'.format(run))
        #shock_df.to_csv(work_path + '/data/output/shocks/shocks_99_aggregated.csv')
        
        return
    
    # def __set_aff_hhs_disaster_set_damage(self, df_haz, df_hh, hh_reg, k_eff_reg):

    #     region_dict = {'PH150000000': 15,
    #                    'PH140000000': 14,
    #                    'PH130000000': 13,
    #                    'PH010000000': 1,
    #                    'PH020000000': 2,
    #                    'PH030000000': 3,
    #                    'PH040000000': 41,
    #                    'PH170000000': 42,
    #                    'PH090000000': 9,
    #                    'PH050000000': 5,
    #                    'PH060000000': 6,
    #                    'PH070000000': 7,
    #                    'PH080000000': 8,
    #                    'PH100000000': 10,
    #                    'PH110000000': 11,
    #                    'PH120000000': 12,
    #                    'PH160000000': 16}
        
    #     aff_ids = np.zeros((len(hh_reg), len(self.__time_stemps)))
        
    #     for ev_id, event in enumerate(self.__event_names):
    #         aff_hhs = []
    #         affected_cells = df_haz.loc[df_haz[event] !=0, 'Centroid_ID']
    #         for cell in affected_cells:
    #             hhids = np.array(df_hh.loc[df_hh['Centroid_ID']==cell, 'fhhid'])
    #             min_k_eff = 0
    #             sum_k_eff = k_eff_reg
    #             rel_dam = df_haz.loc[(df_haz['Centroid_ID']==cell),event].values[0]
    #             abs_dam = rel_dam * sum_k_eff
    #             dam = 0
    #             print(abs_dam)
    #             while dam<abs_dam:
    #                 print(dam)
    #                 print(ev_id)
    #                 print(cell)
    #                 hh_ind = random.randint(0, hhids.shape[0]-1)
    #                 if not np.isin(hhids[hh_ind], aff_hhs):
    #                     aff_hhs.append(hhids[hh_ind])
    #                     dam += hh_reg[hhids[hh_ind]].weight*hh_reg[hhids[hh_ind]].vul*(hh_reg[hhids[hh_ind]].k_eff_0 - hh_reg[hhids[hh_ind]].d_k_eff_t)
    #                 all_aff = np.isin(hhids, aff_hhs).sum()
    #                 if all_aff == hhids.shape[0]:
    #                     hh_ind = random.randint(0, hhids.shape[0]-1)
    #                     new_vul = hh_reg[hhids[hh_ind]].vul+0.1
    #                     if new_vul < 0.6:
    #                         hh_reg[hhids[hh_ind]].set_vul(new_vul)
    #                         dam+= 0.1* hh_reg[hhids[hh_ind]].weight*(hh_reg[hhids[hh_ind]].k_eff_0 - hh_reg[hhids[hh_ind]].d_k_eff_t)
                        
    #         aff_ids[aff_hhs, ev_id]=1
    #     return
    
    def __set_aff_hhs_disaster_set(self, work_path, path_haz, path_hh, reg, event_names):

        print('region')
        print(reg)
        df_haz = pd.read_csv((work_path + path_haz).format(reg))
        df_hh = pd.read_csv(work_path + path_hh)
        df_hh = df_hh[df_hh['region'] == reg]
        for ev_id, event in enumerate(event_names):
            print(event)
            aff_hhs = []
            affected_cells = df_haz.loc[df_haz[event] !=0, 'Centroid_ID']
            if affected_cells.shape[0]==0:
                continue
            for cell in affected_cells:
                hhids = list(df_hh.loc[df_hh['centroid_id'] == cell, 'fhhid'])
                
                n_people = df_hh.loc[df_hh['centroid_id'] == cell, 'weight'].sum()
                if n_people == 0.0:
                    continue
                frac = df_haz.loc[(df_haz['Centroid_ID'] == cell), event].values[0]
                n_aff_people = np.round(frac* n_people)
                c_n_aff_hhs = 0
                while (c_n_aff_hhs < n_aff_people) and (len(hhids)>0):
                    hh_ind = random.choice(hhids)
                    aff_hhs.append(hh_ind)
                    c_n_aff_hhs += df_hh.loc[df_hh['fhhid']==hh_ind].n_individuals.values[0]
                    hhids.remove(hh_ind)
                
            self.__aff_ids[aff_hhs, ev_id]=1
                
        return

    # def __get_k_eff_reg(self, hh_reg):

    #     for reg in hh_reg.region_hhs:
            
    #         for hhid in hh_reg.region_hhs[reg]:
                

    #         #print('hä')


    def __set_aff_hhs(self, n_hhs):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade intersection.
    
            Parameters
            ----------
            hhs : list
                list with all households
        """
        n_aff_hh = 6
        sel_hhs = []

        # while len(sel_hhs) < n_aff_hh:
        #     hh = random.randint(0, n_hhs-1)
        #     if not np.isin(hh, sel_hhs):
        #         sel_hhs.append(hh)
        #         print('Household appended')

        return [0,1,2,3,4,5]

        
    def __set_time_stemps(self, n_events=2):
        """The function selects households randomly hocks them. (Function is only a
           placeholder for a real Climade intersection.

            Parameters
            ----------
            n_events : int
                number events
        """

        self.__time_stemps = [50]
        self.__event_names = ['0']
        add = 0
        for run in range(n_events-1):
            self.__time_stemps.append(random.randint(200+add, 300+add))
            self.__event_names.append(str(self.__time_stemps[run]))
            add+=300
     

    def shock(self, event_index, gov, hhs, dt_reco, cores):
        """The function selects households randomly and shocks them. (Function is only a
           placeholder for a real Climade interface.

            Parameters
            ----------
            gov : Government
                the government of all households
        """
        
        affected_hhs = np.where(self.__aff_ids[:, event_index]==1)[0]
        L_t = 0
        p = mp.Pool(cores)
        for h_ind, hh in enumerate(hhs):
            if hh.hhid in affected_hhs:
                L_t += hh.vul * (hh.k_eff_0 - hh.d_k_eff_t)
            else:
                L_t += hh.d_k_eff_t
                
        prod_x = partial(Shock.apply_shock, L=L_t, K=gov.K, dt_reco=dt_reco, affected_hhs=affected_hhs)
        
        hhs = p.map(prod_x, hhs)
        p.close()
        p.join()
        # print('shock shocks')
        # for h_ind, hh in enumerate(hhs):
        #     if h_ind in affected_hhs:
        #         hh.shock(aff_flag=True, L=L_t, K=gov.K, dt=dt_reco)
        #     else:
        #         hh.shock(aff_flag=False, L=L_t, K=gov.K, dt=dt_reco)

        gov.shock(aff_flag=True, L=L_t)

        return hhs

    def __plot(self, ax):
        """The function selects households randomly and plots their recovery.

            Parameters
            ----------
            gov : Government
                the government of all households
            hhs : Households
                list with households
        """

        ax.set_xlim((0, RECO_PERIOD*DT_STEP+1))
        ax.set_xlabel('weeks after shock')
        ax.set_ylabel('')
        ax.set_title('Consumption_recovery')
        #ax.legend()

        return
