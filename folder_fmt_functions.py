# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 22:11:59 2018

@author: Marios D. Pazidis
"""
import os
import math
import numpy as np

def folder_init(mainFolder):
    
    # general folders
    Combinations=mainFolder+'\Combinations'
    Sa_unsc_ave=mainFolder+'\Sa_unsc_ave'
    
    # population folders
    CF_0=mainFolder+'\Population\All\CF_0'
    CF_1=mainFolder+'\Population\All\CF_1'
    Par_F=mainFolder+'\Population\All\Par_F'
    Par_CR=mainFolder+'\Population\All\Par_CR'
    Scaling_factors=mainFolder+'\Population\All\Scaling_factors'
    
    # 1st Pareto front folders for every combination
    CF_0_Pareto=mainFolder+'\Population\First_Pareto_front\CF_0'
    CF_1_Pareto=mainFolder+'\Population\First_Pareto_front\CF_1'
    Par_F_Pareto=mainFolder+'\Population\First_Pareto_front\Par_F'
    Par_CR_Pareto=mainFolder+'\Population\First_Pareto_front\Par_CR'
    Scaling_factors_Pareto=mainFolder+'\Population\First_Pareto_front\Scaling_factors'
    
    # Global 1st Pareto front folders of all combinations
    CF_0_Pareto_GL=mainFolder+'\Global_First_Pareto_front\CF_0'
    CF_1_Pareto_GL=mainFolder+'\Global_First_Pareto_front\CF_1'
    Par_F_Pareto_GL=mainFolder+'\Global_First_Pareto_front\Par_F'
    Par_CR_Pareto_GL=mainFolder+'\Global_First_Pareto_front\Par_CR'
    Scaling_factors_Pareto_GL=mainFolder+'\Global_First_Pareto_front\Scaling_factors'
    
    key=['Combinations','Sa_unsc_ave', 'CF_0', 'CF_1', 'Par_F', 'Par_CR', 'Scaling_factors', 'CF_0_Pareto', 'CF_1_Pareto', 'Par_F_Pareto', 'Par_CR_Pareto', 'Scaling_factors_Pareto']
    value=[Combinations,Sa_unsc_ave, CF_0, CF_1, Par_F, Par_CR, Scaling_factors, CF_0_Pareto, CF_1_Pareto, Par_F_Pareto, Par_CR_Pareto, Scaling_factors_Pareto, CF_0_Pareto_GL, CF_1_Pareto_GL, Par_F_Pareto_GL, Par_CR_Pareto_GL, Scaling_factors_Pareto_GL]
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