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

class ClimateLife():

    def __init__(self, hhs, shock, gov):
        
        self.__hhs = hhs
        self.__shock = shock
        self.__gov = gov
        self.__pt = np.array([])
        self.__dt_life = np.array([])
        
        plt.figure(figsize=(20, 10))
        plt.suptitle('Household Resilience Model', size=20)
        
        self.__info_summary = gridspec.GridSpec(3, 3)
        self.__info_gov = plt.subplot(self.__info_summary[0, :])
        self.__abs_wb = plt.subplot(self.__info_summary[1, :-1])
        self.__abs_cons = plt.subplot(self.__info_summary[1, -1])
        self.__abs_k = plt.subplot(self.__info_summary[-1, -1])
        self.__abs_inc = plt.subplot(self.__info_summary[-1, 0])
        self.__abs_inc_sp = plt.subplot(self.__info_summary[-1, -2])

    def start(self):
        
        print('Life started')
        
        pt = np.linspace(0, RECO_PERIOD, RECO_PERIOD*DT_STEP+1)
        dt_reco = np.diff(pt)[0]
        
        self.__pt = pt[::4]
        
        self.__dt_life = np.arange(0, RECO_PERIOD*DT_STEP+1)

        print('Tax rate: ' + str(self.__gov))
        print('Total expenditure on social programs: ' + str(self.__gov.sp_cost))
        print('Total national capital stock: ' + str(self.__gov.K))

        for hh in self.__hhs:
            print('Capital stock of HH ' + str(int(hh.hhid))+': '+str(hh.k_eff_0))

        print('Lambda HH 1: ' + str(self.__hhs[1].lmbda) + ' per year')
        print('Tau HH 1: ' + str(self.__hhs[1].tau) + ' years')

        plot_ids = self.__get_plot_hhs()

        #gov.shock()

        dt_s = 0

        plt.ion()
        n_shock = 0

        for t_i in self.__dt_life:
            
            if t_i in self.__shock.time_stemps:
                self.__shock.shock(n_shock, self.__gov, self.__hhs, dt_reco)
                n_shock += 1
                dt_s = 0

            self.__gov.update_reco(t_i, self.__hhs)
            for hh in self.__hhs:
                hh.update_reco(t_i, self.__gov.L_t, self.__gov.K)

            dt_s += 1

            self.__plot_info(n_plot_hhs=5, plot_hhs=plot_ids)

            if t_i%12 ==0:
                #plt.tight_layout()
                plt.show(block=False)
                plt.pause(0.01)
        return

    def _set_agents(self):

        return
    
    def __get_plot_hhs(self, n_plot_hhs=5):
        
        sort_aff = np.argsort(self.__shock.aff_ids.sum(axis=1))
        plot_hhs = sort_aff[-n_plot_hhs:]

        return plot_hhs
    
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
        self.__abs_cons.set_ylim((0, 30000))
        self.__abs_cons.set_title('HH total consumption loss USD')
        self.__abs_cons.set_ylabel('Absolute loss in USD')
        #self.__abs_cons.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']

        for h, phh in enumerate(plot_hhs):
            self.__abs_cons.plot(self.__pt, self.__hhs[phh].cons_reco, color = colours[h])


    def __plot_inc(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_inc.clear()
        self.__abs_inc.set_xlim((-1, RECO_PERIOD))
        self.__abs_inc.set_ylim((0, 20000))
        self.__abs_inc.set_title('HH total income loss USD')
        self.__abs_inc.set_ylabel('Absolute loss in USD')
        self.__abs_inc.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']
        
        for h, phh in enumerate(plot_hhs):
            self.__abs_inc.plot(self.__pt, self.__hhs[phh].inc_reco, color = colours[h],
                               label='HH' + str(self.__hhs[phh].hhid))
        
        self.__abs_inc.legend()
        
    def __plot_inc_sp(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_inc_sp.clear()
        self.__abs_inc_sp.set_xlim((-1, RECO_PERIOD))
        self.__abs_inc_sp.set_ylim((0, 3000))
        self.__abs_inc_sp.set_title('HH income loss from social programs USD')
        self.__abs_inc_sp.set_ylabel('Absolute loss in USD')
        self.__abs_inc_sp.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']
        
        for h, phh in enumerate(plot_hhs):
            self.__abs_inc_sp.plot(self.__pt, self.__hhs[phh].inc_sp_reco, color = colours[h],
                               label='HH' + str(self.__hhs[phh].hhid))
        
        
    def __plot_k_eff(self, n_plot_hhs=4, plot_hhs=None):
        self.__abs_k.clear()
        self.__abs_k.set_xlim((-1, RECO_PERIOD))
        self.__abs_k.set_ylim((0, 80000))
        self.__abs_k.set_title('HH capital stock damage')
        self.__abs_k.set_ylabel('Absolute damage in USD')
        self.__abs_k.set_xlabel('time in years')
        colours = ['blue', 'red', 'green', 'gold', 'brown', 'black', 'purple', 'pink','seagreen', 'firebrick']
        
        for h, phh in enumerate(plot_hhs):
            self.__abs_k.plot(self.__pt, self.__hhs[phh].k_eff_reco, color = colours[h],
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
            self.__abs_wb.plot(self.__pt, self.__hhs[phh].wb_reco, color = colours[h],
                               label='HH' + str(self.__hhs[phh].hhid) + ' dec: ' + str(self.__hhs[phh].decile))
        
        self.__abs_wb.legend()
        
    def __plot_info_gov(self):
        self.__info_gov.clear()
        self.__info_gov.set_xlim((-1, RECO_PERIOD))
        self.__info_gov.set_ylim((0, 250000))
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
