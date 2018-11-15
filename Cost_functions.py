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
