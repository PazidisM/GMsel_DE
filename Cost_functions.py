# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 14:39:00 2018

@author: Marios D. Pazidis
"""

import numpy as np

def Obj_RMSE(Sa_ave,Sa_Tgt,sf_ave):
    
    diff=Sa_ave+sf_ave-Sa_Tgt
    z=np.sqrt(np.sum(np.power(diff,2))/np.size(Sa_Tgt,1))
    
    return z

def Obj_std(Sa_suite,sf_suite,Sa_unsc_suite,sf_suite_mean):
    
    sf_suite.shape=(len(sf_suite),1)
    diff=Sa_suite+sf_suite-(Sa_unsc_suite+sf_suite_mean)
    var=np.power(diff,2)/(len(sf_suite)-1)
    var_mean=np.mean(var)
    z=np.sqrt(var_mean)
    
    return z