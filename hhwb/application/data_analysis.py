#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 19:31:24 2022

@author: insauer
"""
import numpy as np
import pandas as pd

class DataAnalysis():
    
    def __init__(self, survey_data_path='', shock_data_path='', output_data_path='', column_id='', run_name=''):

        self.__output_data_path = output_data_path
        self.__survey_data_path = survey_data_path
        self.__hhs = pd.read_csv(self.__survey_data_path)
        self.__shocks = pd.read_csv(shock_data_path)
        del_cols = [col for col in self.__shocks.columns if 'Unnamed' in col]
        self.__shocks = self.__shocks.drop(columns=del_cols)
        self.__shocks['fhhid'] = self.__shocks.index
        self.__shocks['n_events'] = self.__shocks.iloc[:,:-2].sum(axis=1)
        self.__hhs['n_events']= self.__shocks['n_events']
        self.__hhs['hh_weight']= self.__hhs['weight']/self.__hhs['n_individuals']
        self.__column_id = column_id
        self.__run_name = run_name
        
    def analyse_time(self, step=20000):
        
        col=0
        add=step
        while add > 0:
            print(col)
            cols=np.arange(col,col+add)
            sub_list=list(self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'subsistence_line']/13)
            cons_list=list(self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'income']/13)
            df=pd.read_csv(self.__output_data_path + 'cons.csv', usecols=cols, header=None)
            
            for i,hhid in enumerate(cols):
                #print(hhid)
                self.__hhs.loc[self.__hhs['fhhid']==hhid,
                               'time_under_sub{}'.format(self.__column_id)]=np.array(cons_list[i]-df.loc[:,hhid]<sub_list[i]).sum()
                self.__hhs.loc[self.__hhs['fhhid']==hhid,
                               'reco_time{}'.format(self.__column_id)]=np.array(cons_list[i]-df.loc[:,hhid]<0.95*cons_list[i]).sum()
            col+=add
            if col+add>=self.__hhs.shape[0]:
                add=self.__hhs.shape[0]-col
            else:
                add=step
                
        col=0
        add=step
        while add > 0:
            print(col)
            cols=np.arange(col,col+add)
            sub_list=list(self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'subsistence_line']/13)
            cons_list=list(self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'income']/13)
            df=pd.read_csv(self.__output_data_path + 'cons_sm.csv', usecols=cols, header=None)
            
            for i,hhid in enumerate(cols):
                #print(hhid)
                self.__hhs.loc[self.__hhs['fhhid']==hhid,
                               'time_under_sub_sm{}'.format(self.__column_id)]=np.array(cons_list[i]-df.loc[:,hhid]<sub_list[i]).sum()
                self.__hhs.loc[self.__hhs['fhhid']==hhid,
                               'reco_time_sm{}'.format(self.__column_id)]=np.array(cons_list[i]-df.loc[:,hhid]<0.95*cons_list[i]).sum()
            col+=add
            if col+add>=self.__hhs.shape[0]:
                add=self.__hhs.shape[0]-col
            else:
                add=step
                
        self.__hhs.to_csv('survey_'+self.__run_name+'.csv')
        
    def analyse_wb(self, step=20000):
        
        col=0
        add=step
        while add > 0:
            print(col)
            cols=np.arange(col,col+add)

            df=pd.read_csv(self.__output_data_path + 'wb.csv', usecols=cols, header=None)
            
            self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'wb_loss{}'.format(self.__column_id)]=np.array(df.max())
            col+=add
            if col+add>=self.__hhs.shape[0]:
                add=self.__hhs.shape[0]-col
            else:
                add=step
        col=0
        add=step
        while add > 0:
            print(col)
            cols=np.arange(col,col+add)

            df=pd.read_csv(self.__output_data_path + 'wb_sm.csv', usecols=cols, header=None)
            
            self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'wb_loss_sm{}'.format(self.__column_id)]=np.array(df.max())
            col+=add
            if col+add>=self.__hhs.shape[0]:
                add=self.__hhs.shape[0]-col
            else:
                add=step
            
        self.__hhs.to_csv('survey_'+self.__run_name+'.csv')
            
                    
        
        