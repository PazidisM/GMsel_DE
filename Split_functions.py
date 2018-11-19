# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 20:15:12 2018

@author: Marios D. Pazidis
"""

import numpy as np

def Pre_run_split(split_data,NSeed,Combs,Sa_unsc_ave,folders,formats):
    
    split_size=split_data['split_size']
    
    index=0
    while index<NSeed:
        split_start=index%split_size+index//split_size*split_size
        split_end=split_start+split_size
        
        if split_end>NSeed:
            split_end=NSeed
            
        Combs_split=Combs[split_start:split_end]
        fileName=folders['Combinations']+'\Combs_'+str(index//split_size).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Combs_split,formats['fmt_combs'])
        
        Sa_unsc_ave_split=Sa_unsc_ave[split_start:split_end]
        fileName=folders['Sa_unsc_ave']+'\Sa_unsc_ave_'+str(index//split_size).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Sa_unsc_ave_split,formats['fmt_Sa_unsc_ave'])
        
        index=split_end