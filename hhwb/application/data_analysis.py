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
       
        self.__hhs.to_csv('survey_'+self.__run_name+'_analysed.csv')
        
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
            
        self.__hhs.to_csv('survey_'+self.__run_name+'_analysed.csv')

        
    def analyse_time_steps(self, step=10000):
        
        col=0
        add=step
        
        #month=[32, 34, 38, 43, 53, 78, 80, 101, 160, 169, 182, 240]

        month = [  7,  19,  20,  26,  28,  32,  34,  38, 43, 47,  48,  53,  59,  64,
             66,  76, 78, 80,  83,  86,  91,  92,  99, 101, 110, 117, 118, 123,
            127, 130, 132, 133, 142, 144, 160, 162, 166, 169, 180, 182, 201, 240]
        while add > 0:
            print(col)
            cols=np.arange(col,col+add)
            df_keff=pd.read_csv(self.__output_data_path + 'keff.csv', usecols=cols, header=None)
            #df_inc=pd.read_csv(self.__output_data_path + 'inc.csv', usecols=cols, header=None)
            #df_inc_sp=pd.read_csv(self.__output_data_path + 'inc_sp.csv', usecols=cols, header=None)
            for m in month:
                
                self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'keff_{}'.format(str(m))]=np.array(df_keff.iloc[m,:])
                self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'keff_diff_{}'.format(str(m))]=np.array(df_keff.diff().iloc[m,:])
                #self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'inc_{}'.format(str(m))]=np.array(df_inc.iloc[m,:])
                #self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'inc_sp_{}'.format(str(m))]=np.array(df_inc_sp.iloc[m,:])

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
            #df_cons=pd.read_csv(self.__output_data_path + 'cons.csv', usecols=cols, header=None)
            df_cons_sm=pd.read_csv(self.__output_data_path + 'cons_sm.csv', usecols=cols, header=None)
            df_wb=pd.read_csv(self.__output_data_path + 'wb.csv', usecols=cols, header=None)
            df_wb_sm=pd.read_csv(self.__output_data_path + 'wb_sm.csv', usecols=cols, header=None)
            for m in month:
                
                #self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'cons_{}'.format(str(m))]=np.array(df_cons.iloc[m,:])
                self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'cons_sm_{}'.format(str(m))]=np.array(df_cons_sm.iloc[m,:])
                self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'wb_{}'.format(str(m))]=np.array(df_wb.iloc[m,:])
                self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'wb_sm_{}'.format(str(m))]=np.array(df_wb_sm.iloc[m,:])
                self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'wb_inc_{}'.format(str(m))]=np.array(df_wb.diff().iloc[m+1,:])
                self.__hhs.loc[self.__hhs['fhhid'].isin(cols),'wb_inc_sm_{}'.format(str(m))]=np.array(df_wb_sm.diff().iloc[m+1,:])
            col+=add
            if col+add>=self.__hhs.shape[0]:
                add=self.__hhs.shape[0]-col
            else:
                add=step
            
        self.__hhs.to_csv('survey_'+self.__run_name+'_analysed.csv')
        
            
                    
        
        