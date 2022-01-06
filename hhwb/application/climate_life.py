#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 21:19:13 2020

@author: insauer
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import multiprocessing as mp
from functools import partial
import time
from hhwb.agents.government import Government
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD
import pandas as pd
import csv



class ClimateLife():

    def __init__(self, hhs, shock, gov):
        
        self.__hhs = hhs
        self.__shock = shock
        self.__gov = gov
        self.__pt = np.array([])
        self.__dt_life = np.array([])
        
        #plt.figure(figsize=(10, 5))
        #plt.suptitle('Household Resilience Model', size=20)
        
        # self.__info_summary = gridspec.GridSpec(3, 3)
        # self.__info_gov = plt.subplot(self.__info_summary[0, :])
        # self.__abs_wb = plt.subplot(self.__info_summary[1, :-1])
        # self.__abs_cons = plt.subplot(self.__info_summary[1, -1])
        # self.__abs_k = plt.subplot(self.__info_summary[-1, -1])
        # self.__abs_inc = plt.subplot(self.__info_summary[-1, 0])
        # self.__abs_inc_sp = plt.subplot(self.__info_summary[-1, -2])
        # save with RECO_PERIOD*DT_STEP/TEMP_RES + 1
        # self.k_eff_reco = np.empty((int((RECO_PERIOD*DT_STEP)/TEMP_RES +1), len(self.__hhs)))
        # self.k_eff_reco[:] = np.nan
        # self.inc_reco = np.empty((int((RECO_PERIOD*DT_STEP)/TEMP_RES +1), len(self.__hhs)))
        # self.inc_reco[:] = np.nan
        # self.inc_sp_reco = np.empty((int((RECO_PERIOD*DT_STEP)/TEMP_RES +1), len(self.__hhs)))
        # self.inc_sp_reco[:] = np.nan
        # self.cons_reco = np.empty((int((RECO_PERIOD*DT_STEP)/TEMP_RES +1), len(self.__hhs)))
        # self.cons_reco[:] = np.nan
        # self.wb_reco = np.empty((int((RECO_PERIOD*DT_STEP)/TEMP_RES +1), len(self.__hhs)))
        # self.wb_reco[:] = np.nan
        
        # self.cons_reco_sm = np.empty((int((RECO_PERIOD*DT_STEP)/TEMP_RES +1), len(self.__hhs)))
        # self.cons_reco_sm[:] = np.nan
        
        
    @property
    def dt_life(self):
        
        return self.__dt_life
    @staticmethod
    def update_reco(hh, gov):
        hh.update_reco(gov.L_t, gov.K)
        
        return hh
    
    # def __update_reco(self, t_i):

    #     for hh in self.__hhs:
            
    #         if not t_i in self.__shock.time_stemps:
    #             hh.update_reco(t_i, self.__gov.L_t, self.__gov.K)

    #         if t_i % TEMP_RES == 0:
    #             self.k_eff_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_k_eff_t
    #             self.inc_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_inc_t
    #             self.inc_sp_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_inc_sp_t
    #             self.cons_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_con_t
    #             self.wb_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.wb_smooth
    #             self.cons_reco_sm[int(t_i/TEMP_RES), int(hh.hhid)] = hh.con_smooth
    #         hh._t += hh._dt

    #     return
    
    def __update_records(self, t_i):
        
        keff = np.zeros((len(self.__hhs)))
        inc = np.zeros((len(self.__hhs)))
        inc_sp = np.zeros((len(self.__hhs)))
        cons = np.zeros((len(self.__hhs)))
        cons_sm = np.zeros((len(self.__hhs)))
        wb = np.zeros((len(self.__hhs)))
        
        for hh in self.__hhs:

            if t_i % TEMP_RES == 0:
                keff[int(hh.hhid)] = hh.d_k_eff_t
                inc[int(hh.hhid)] = hh.d_inc_t
                inc_sp[int(hh.hhid)] = hh.d_inc_sp_t
                cons[int(hh.hhid)]  = hh.d_con_t
                wb[int(hh.hhid)] = hh.wb_smooth
                cons_sm[int(hh.hhid)] = hh.con_smooth
            hh._t += hh._dt

        return keff, inc, inc_sp, cons, cons_sm, wb
    
    

    def start(self, work_path='/home/insauer/projects/WB_model/hhwb',
              result_path='/data/output/', cores=1):
        
        print('Life started')
        
        
        
        pt = np.linspace(0, RECO_PERIOD, RECO_PERIOD*DT_STEP+1)
        dt_reco = np.diff(pt)[0]
        
        self.__pt = pt[::4]
        
        self.__dt_life = np.arange(0, RECO_PERIOD*DT_STEP+1)
        

        print('Tax rate: ' + str(self.__gov))
        print('Total expenditure on social programs: ' + str(self.__gov.sp_cost))
        print('Total national capital stock: ' + str(self.__gov.K))
        
        colnames = np.arange(len(self.__hhs)).astype(str)
        
        keff = pd.DataFrame(columns=colnames)
        inc =  pd.DataFrame(columns=colnames)
        inc_sp =  pd.DataFrame(columns=colnames)
        cons =  pd.DataFrame(columns=colnames)
        cons_sm =  pd.DataFrame(columns=colnames)
        wb =  pd.DataFrame(columns=colnames)
        
        keff.to_csv(work_path+result_path+'keff.csv')
        inc.to_csv(work_path+result_path+'inc.csv')
        inc_sp.to_csv(work_path+result_path+'inc_sp.csv')
        cons.to_csv(work_path+result_path+'cons.csv')
        cons_sm.to_csv(work_path+result_path+'cons_sm.csv')
        wb.to_csv(work_path+result_path+'wb.csv')
        
        for hh in self.__hhs:
            print('Capital stock of HH ' + str(int(hh.hhid))+': '+str(hh.k_eff_0))

        
        #plot_ids = self.__get_plot_hhs()

        #gov.shock()

        #dt_s = 0

        #plt.ion()
        n_shock = 0
        
        
        with open(work_path+result_path+'keff.csv', 'w') as f_keff,\
             open(work_path+result_path+'inc.csv', 'w') as f_inc,\
             open(work_path+result_path+'inc_sp.csv', 'w') as f_inc_sp,\
             open(work_path+result_path+'cons.csv', 'w') as f_cons,\
             open(work_path+result_path+'cons_sm.csv', 'w') as f_cons_sm,\
             open(work_path+result_path+'wb.csv', 'w') as f_wb:
                 
            writer_keff=csv.writer(f_keff, delimiter=',')
            writer_inc=csv.writer(f_inc, delimiter=',')
            write_incsp=csv.writer(f_inc_sp, delimiter=',')
            writer_cons=csv.writer(f_cons, delimiter=',')
            writer_conssm=csv.writer(f_cons_sm, delimiter=',')
            writer_wb=csv.writer(f_wb, delimiter=',')
            
            for t_i in self.__dt_life:
    
                print(t_i)
    
                if t_i in self.__shock.time_stemps:
                    print('shock start')
                    self.__hhs = self.__shock.shock(n_shock, self.__gov, self.__hhs, dt_reco, cores)
                    n_shock += 1
                    #dt_s = 0
                self.__gov.update_reco(t_i, self.__hhs)
                
                if not t_i in self.__shock.time_stemps:
                    p = mp.Pool(cores)
                    prod_x=partial(ClimateLife.update_reco, gov=self.__gov)
                    self.__hhs=p.map(prod_x, self.__hhs)
                    p.close()
                    p.join()
    
                keff, inc, inc_sp, cons, cons_sm, wb = self.__update_records(t_i)
    

                writer_keff.writerow(keff)
                writer_inc.writerow(inc)
                write_incsp.writerow(inc_sp)
                writer_cons.writerow(cons)
                writer_conssm.writerow(cons_sm)
                writer_wb.writerow(wb)
            
            f_keff.close()
            f_inc.close()
            f_inc_sp.close()
            f_cons.close()
            f_cons_sm.close()
            f_wb.close()
    
                # self.__plot_info(n_plot_hhs=5, plot_hhs=plot_ids)
    
                # if t_i%12 ==0:
                #     #plt.tight_layout()
                #     plt.show(block=False)
                #     plt.pause(0.01)
        return

    def _set_agents(self):

        return
    
    def __get_plot_hhs(self, n_plot_hhs=10):
        
        sort_aff = np.argsort(self.__shock.aff_ids.sum(axis=1))
        plot_hhs = sort_aff[-n_plot_hhs:]
        #print(plot_hhs)
        return [0, 1, 5, 7]
    
    def __plot_info(self, n_plot_hhs=4, plot_hhs=None):

        self.__plot_cons(n_plot_hhs=4, plot_hhs=plot_hhs)
        self.__plot_inc(n_plot_hhs=4, plot_hhs=plot_hhs)
        self.__plot_inc_sp(n_plot_hhs=4, plot_hhs=plot_hhs)
        self.__plot_k_eff(n_plot_hhs=4, plot_hhs=plot_hhs)
        self.__plot_wb(n_plot_hhs=4, plot_hhs=plot_hhs)
        self.__plot_info_gov()

        #self.__info_summary.suptitle('Household Resilience Model', va= 'top')

    def __plot_cons(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_cons.clear()
        self.__abs_cons.set_xlim((-1, RECO_PERIOD))
        #self.__abs_cons.set_ylim((0, 30000))
        self.__abs_cons.set_title('HH total consumption loss USD')
        self.__abs_cons.set_ylabel('Absolute loss in USD')
        #self.__abs_cons.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']

        for h, phh in enumerate(plot_hhs):
            self.__abs_cons.plot(self.__pt, self.cons_reco[:,phh], color = colours[h])
            # self.__abs_cons.plot(self.__pt, self.cons_reco_sm[:,phh], color = colours[h])
            self.__abs_cons.axhline(y=self.__hhs[phh].consum_0 - self.__hhs[phh].subsistence_line, color = colours[h])

    def __plot_inc(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_inc.clear()
        self.__abs_inc.set_xlim((-1, RECO_PERIOD))
        #self.__abs_inc.set_ylim((0, 20000))
        self.__abs_inc.set_title('HH total income loss USD')
        self.__abs_inc.set_ylabel('Absolute loss in USD')
        self.__abs_inc.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']
        
        for h, phh in enumerate(plot_hhs):
            self.__abs_inc.plot(self.__pt, self.inc_reco[:, phh], color = colours[h],
                               label='HH' + str(self.__hhs[phh].hhid))
        
        self.__abs_inc.legend()
        
    def __plot_inc_sp(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_inc_sp.clear()
        self.__abs_inc_sp.set_xlim((-1, RECO_PERIOD))
        #self.__abs_inc_sp.set_ylim((0, 3000))
        self.__abs_inc_sp.set_title('HH income loss from social programs USD')
        self.__abs_inc_sp.set_ylabel('Absolute loss in USD')
        self.__abs_inc_sp.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']
        
        for h, phh in enumerate(plot_hhs):
            self.__abs_inc_sp.plot(self.__pt, self.inc_sp_reco[:, phh], color = colours[h],
                               label='HH' + str(self.__hhs[phh].hhid))
        
        
    def __plot_k_eff(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_k.clear()
        self.__abs_k.set_xlim((-1, RECO_PERIOD))
        #self.__abs_k.set_ylim((0, 80000))
        self.__abs_k.set_title('HH capital stock damage')
        self.__abs_k.set_ylabel('Absolute damage in USD')
        self.__abs_k.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']
        
        for h, phh in enumerate(plot_hhs):
            self.__abs_k.plot(self.__pt, self.k_eff_reco[:, phh], color = colours[h],
                               label='HH' + str(self.__hhs[phh].hhid))
    
    def __plot_wb(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_wb.clear()
        self.__abs_wb.set_xlim((-1, RECO_PERIOD))
        #self.__abs_wb.set_ylim((0, 100000))
        self.__abs_wb.set_title('Accumlated HH WB loss')
        self.__abs_wb.set_ylabel('')
        #self.__abs_wb.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']
        for h, phh in enumerate(plot_hhs):
            self.__abs_wb.plot(self.__pt, self.wb_reco[:, phh], color = colours[h],
                                label='HH' + str(self.__hhs[phh].hhid) + ' dec: ' + str(self.__hhs[phh].decile))
            # self.__abs_wb.plot(self.__pt, self.wb_reco_sm[:, phh], color = colours[h],
            #                     label='HH' + str(self.__hhs[phh].hhid) + ' dec: ' + str(self.__hhs[phh].decile))
        self.__abs_wb.legend()

    def __plot_info_gov(self):
        self.__info_gov.clear()
        self.__info_gov.set_xlim((-1, RECO_PERIOD))
        #self.__info_gov.set_ylim((0, 250000))
        self.__info_gov.set_title('National losses')
        self.__info_gov.set_ylabel('Absolute national loss USD')
        self.__info_gov.set_xlabel('time in years')
        self.__info_gov.plot(self.__pt, self.__gov.k_eff_reco, color = 'blue',
                               label='national capital stock damage', alpha = 0.5)
        self.__info_gov.plot(self.__pt, self.__gov.inc_reco, color = 'red',
                               label='national income loss', alpha = 0.5)
        self.__info_gov.plot(self.__pt, self.__gov.inc_sp_reco, color = 'green',
                               label='national loss in social spendings', alpha = 0.5)
        self.__info_gov.plot(self.__pt, self.__gov.cons_reco, color = 'gold',
                               label='loss in national consumption', alpha = 0.5)
        self.__info_gov.legend()
    
    # def write_output_files(self, t_i, result_path):
        
    #     #path = '/home/insauer/projects/WB_model/hhwb/data/output/'

    #     k_eff = pd.DataFrame(data = self.k_eff_reco)
    #     k_eff.to_csv(result_path + 'k_eff.csv')
    #     del k_eff
    #     inc = pd.DataFrame(data = self.inc_reco)
    #     inc.to_csv(result_path + 'inc.csv')
    #     del inc
    #     inc_sp = pd.DataFrame(data = self.inc_sp_reco)
    #     inc_sp.to_csv(result_path + 'inc_sp.csv')
    #     del inc_sp
    #     hh_wb = pd.DataFrame(data = self.wb_reco)
    #     hh_wb.to_csv(result_path + 'hh_wb.csv')
    #     del hh_wb
    #     cons = pd.DataFrame(data = self.cons_reco)
    #     cons.to_csv(result_path + 'cons.csv')
    #     del cons
    #     cons_sav = pd.DataFrame(data = self.cons_reco_sm)
    #     cons_sav.to_csv(result_path + 'cons_sav.csv')
    #     del cons_sav

    #     return

