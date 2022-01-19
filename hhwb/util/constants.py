#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 22:52:45 2020

@author: insauer
"""

__all__ = ['PI',
           'RHO',
           'ETA',
           'T_RNG',
           'DT_STEP',
           'TEMP_RES',
           'RECO_PERIOD']


PI = 0.33
"""Productivity of capital (value currently copied from the original model)"""
RHO = 0.06
"""Utility discount rate (value currently copied from the original model)"""
ETA = 1.5
"""Elasticity of the marginal utility of consumption
   (value currently copied from the original model)"""
T_RNG = 15
"""Time frame after disaster in years where reconstruction rate is optimized
   (value currently copied from the original model)"""
DT_STEP = 13
"""Time step used for optimization (52 : week)
   (value currently copied from the original model)"""
TEMP_RES = 1
"""Temporal resolution after recovery in weeks
   (value currently copied from the original model)"""
RECO_PERIOD = 146
"""Time frame after first disaster in years where reconstruction is modeled (in years)"""

DT = RECO_PERIOD / (RECO_PERIOD * DT_STEP)
"""Time differnece in years between two time stemps"""
SUBS_SAV_RATE = 3339/13. # 1853 PHL2  #4941/52. MWI
"""Minimum recoery and savings rate"""
OPT_DYN = False
"""Savings and recovery optimization in each timestemp"""