# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 22:11:59 2018

@author: Marios D. Pazidis
"""
import os
import math
import numpy as np

def folder_init(mainFolder):
    Combinations=mainFolder+'\Combinations'
    Sa_unsc_ave=mainFolder+'\Sa_unsc_ave'
    Cost_Obj_01=mainFolder+'\Cost_obj_01'
    Cost_Obj_02=mainFolder+'\Cost_obj_02'
    Par_F=mainFolder+'\Par_F'
    Par_CR=mainFolder+'\Par_CR'
    Scaling_factors=mainFolder+'\Scaling_factors'
    key=['Combinations','Sa_unsc_ave', 'Cost_Obj_01', 'Cost_Obj_02', 'Par_F', 'Par_CR', 'Scaling_factors']
    value=[Combinations,Sa_unsc_ave, Cost_Obj_01, Cost_Obj_02, Par_F, Par_CR, Scaling_factors]
    temp=list(zip(key, value))
    folders=dict(temp)
    
    for folder in folders:
        os.makedirs(folders[folder], exist_ok=True)
    
    return folders


def fmt(mainFolder,split_data,NSeed,Rec_db_metadata,DE_par):
    nPop=DE_par['nPop']
    split_size=split_data['split_size']
    split_N=NSeed//split_size
    fill_fn_split=int(math.log10(split_N))+2  # at least 2 digits for file names
    fill_fn_pop=int(math.log10(nPop))+1
    fill_fn_all=int(math.log10(NSeed))+1
    fill_Par_F_CR=3
    fill_combs=int(math.log10(np.max(Rec_db_metadata['idx'])))+2  # space for printed values
    fmt_combs='%-'+str(fill_combs)+'.1d'
    fmt_Sa_unsc_ave='%-.5f'
    fmt_Par_F_CR='%-.3f'
    fmt_Cost='%-.5f'
    fmt_sf='%-10.6f'

    key=['fill_fn_split','fill_fn_all','fill_combs','fmt_combs','fmt_Sa_unsc_ave','fill_fn_pop','fmt_Par_F_CR','fill_Par_F_CR','fmt_Cost','fmt_sf']
    value=[fill_fn_split,fill_fn_all,fill_combs,fmt_combs,fmt_Sa_unsc_ave,fill_fn_pop,fmt_Par_F_CR,fill_Par_F_CR,fmt_Cost,fmt_sf]
    temp=list(zip(key, value))
    formats=dict(temp)
    
    return formats