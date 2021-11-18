#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 21:19:13 2020

@author: insauer
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import time
from hhwb.agents.government import Government
from hhwb.util.constants import PI, RHO, ETA, T_RNG, DT_STEP, TEMP_RES, RECO_PERIOD
import pandas as pd

class ClimateLife():

    def __init__(self, hhs, shock, gov):
        
        self.__hhs = hhs
        self.__shock = shock
        self.__gov = gov
        self.__pt = np.array([])
        self.__dt_life = np.array([])
        
        plt.figure(figsize=(10, 5))
        plt.suptitle('Household Resilience Model', size=20)
        
        self.__info_summary = gridspec.GridSpec(3, 3)
        self.__info_gov = plt.subplot(self.__info_summary[0, :])
        self.__abs_wb = plt.subplot(self.__info_summary[1, :-1])
        self.__abs_cons = plt.subplot(self.__info_summary[1, -1])
        self.__abs_k = plt.subplot(self.__info_summary[-1, -1])
        self.__abs_inc = plt.subplot(self.__info_summary[-1, 0])
        self.__abs_inc_sp = plt.subplot(self.__info_summary[-1, -2])
        
        self.k_eff_reco = np.empty((int(RECO_PERIOD*DT_STEP/TEMP_RES + 1), len(self.__hhs)))
        self.k_eff_reco[:] = np.nan
        self.inc_reco = np.empty((int(RECO_PERIOD*DT_STEP/TEMP_RES + 1), len(self.__hhs)))
        self.inc_reco[:] = np.nan
        self.inc_sp_reco = np.empty((int(RECO_PERIOD*DT_STEP/TEMP_RES + 1), len(self.__hhs)))
        self.inc_sp_reco[:] = np.nan
        self.cons_reco = np.empty((int(RECO_PERIOD*DT_STEP/TEMP_RES + 1), len(self.__hhs)))
        self.cons_reco[:] = np.nan
        self.wb_reco = np.empty((int(RECO_PERIOD*DT_STEP/TEMP_RES + 1), len(self.__hhs)))
        self.wb_reco[:] = np.nan
        
        self.cons_reco_sm = np.empty((int(RECO_PERIOD*DT_STEP/TEMP_RES + 1), len(self.__hhs)))
        self.cons_reco_sm[:] = np.nan
        
    @property
    def dt_life(self):
        
        return self.__dt_life
    
    def __update_reco(self, t_i):

        for hh in self.__hhs:
            
            if not t_i in self.__shock.time_stemps:
                hh.update_reco(t_i, self.__gov.L_t, self.__gov.K)

            if t_i % TEMP_RES == 0:
                self.k_eff_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_k_eff_t
                self.inc_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_inc_t
                self.inc_sp_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_inc_sp_t
                self.cons_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_con_t
                self.wb_reco[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_wb_t
                self.cons_reco_sm[int(t_i/TEMP_RES), int(hh.hhid)] = hh.d_wb_t
            hh._t += hh._dt

        return

    def start(self):
        
        print('Life started')
        
        pt = np.linspace(0, RECO_PERIOD, RECO_PERIOD*DT_STEP+1)
        dt_reco = np.diff(pt)[0]
        
        self.__pt = pt[::4]
        
        self.__dt_life = np.arange(0, RECO_PERIOD*DT_STEP+1)
        
        save_spots=[260,520,1040, 2080]

        print('Tax rate: ' + str(self.__gov))
        print('Total expenditure on social programs: ' + str(self.__gov.sp_cost))
        print('Total national capital stock: ' + str(self.__gov.K))

        for hh in self.__hhs:
            print('Capital stock of HH ' + str(int(hh.hhid))+': '+str(hh.k_eff_0))


        plot_ids = self.__get_plot_hhs()

        #gov.shock()

        dt_s = 0

        plt.ion()
        n_shock = 0

        for t_i in self.__dt_life:

            print(t_i)

            if t_i in self.__shock.time_stemps:
                print('shock start')
                self.__shock.shock(n_shock, self.__gov, self.__hhs, dt_reco)
                n_shock += 1
                #dt_s = 0
            self.__gov.update_reco(t_i, self.__hhs)

            self.__update_reco(t_i)

            if t_i in save_spots:

                self.write_output_files(t_i)

            #dt_s += 1

            self.__plot_info(n_plot_hhs=5, plot_hhs=plot_ids)

            if t_i%12 ==0:
                #plt.tight_layout()
                plt.show(block=False)
                plt.pause(0.01)
        return

    def _set_agents(self):

        return
    
    def __get_plot_hhs(self, n_plot_hhs=10):
        
        sort_aff = np.argsort(self.__shock.aff_ids.sum(axis=1))
        plot_hhs = sort_aff[-n_plot_hhs:]
        #print(plot_hhs)
        return [0,1,2,3,4,5]
    
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
    
    def write_output_files(self, t_i):
        
        path = '/home/insauer/projects/WB_model/hhwb/data/ouput/'

        k_eff = pd.DataFrame(data = self.k_eff_reco)
        k_eff.to_csv(path + 'k_eff.csv')
        del k_eff
        inc = pd.DataFrame(data = self.inc_reco)
        inc.to_csv(path + 'inc.csv')
        del inc
        inc_sp = pd.DataFrame(data = self.inc_sp_reco)
        inc_sp.to_csv(path + 'inc_sp.csv')
        del inc_sp
        hh_wb = pd.DataFrame(data = self.wb_reco)
        hh_wb.to_csv(path + 'hh_wb.csv')
        del hh_wb
        cons = pd.DataFrame(data = self.cons_reco)
        cons.to_csv(path + 'cons.csv')
        del cons

        return
