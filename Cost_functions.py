# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 14:39:00 2018

@author: Marios D. Pazidis
"""

import numpy as np
import numba

# RMSE
def CF_0(Sa_ave,Sa_Tgt,sf_suite):
    
    sf_ave=np.mean(sf_suite)
    diff=Sa_ave+sf_ave-Sa_Tgt
    z=np.sqrt(np.sum(np.power(diff,2))/np.size(Sa_Tgt,1))
    
    return z

# Standard deviation
def CF_1(Sa_suite,sf_suite,Sa_ave):
    
    sf_suite_mean=np.mean(sf_suite)
    sf_suite.shape=(len(sf_suite),1)
    diff=Sa_suite+sf_suite-(Sa_ave+sf_suite_mean)
    var=np.power(diff,2)/(len(sf_suite)-1)
    var_mean=np.mean(var)
    z=np.sqrt(var_mean)
    
    return z

# RMSE
#@numba.jit(nopython=True, parallel=True)
#def CF_0(Sa_unsc_ave_suite,Sa_Tgt,u):
#    
#    z=(np.sum((Sa_unsc_ave_suite+np.sum(u[:,0])/len(u)-Sa_Tgt)**2)/Sa_Tgt.size)**0.5
#    
#    return z
#
## Standard deviation
#@numba.jit(nopython=True, parallel=True)
#def CF_1(Sa_suite,u,Sa_unsc_ave_suite):
#    
#    z=(np.sum((Sa_suite+u-(Sa_unsc_ave_suite+np.sum(u[:,0])/len(u)))**2/(len(u)-1))/len(u))**0.5
#    
#    return z
#





