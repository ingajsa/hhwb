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
           'DT_STEP']


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
DT_STEP = 52
"""Time step used for optimization (52 : week)
   (value currently copied from the original model)"""
