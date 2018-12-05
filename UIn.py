# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:59:30 2018

@author: Marios
"""
import time
import folder_fmt_functions
Inf=1.8e308

def User_Input():
    '''
    #################################   User Input    #################################
    '''
    # selectionParams
    selectionParams={}

    
    databaseFile='ESD_Rexel_meta_data.mat'
    #databaseFile='Resorce2013_meta_data.mat'
    #databaseFile='NGA_W2_meta_data.mat'
    
    nGM=10
    comp_num=2    # search in both components
    comp_idx=1    # used only if comp_num==1
    Tmin=0.08
    Tmax=2
    minScale=0.5
    maxScale=6
    a_g=0.25
    zeta=0.05
    soiltype='A'
    w_CF=[1.0 ,0.0]    # weights of cost functions for initial suites selection
    sameEvent=1     # number of same event records
    sameEventStation=0     # 0: can only contain 1 component per suite, 1: can contain both components
    max_diff_max=6
    max_diff_min=0.25
    nSeed=1
    
    
    # allowedRecs
    Vs30=[850, Inf]
    Mag=[5.0, Inf]
    D=[55, Inf]
    
    # DE parameters
    MaxGen=100
    nPop=100
    F_l=0.2
    F_u=1
    tau_1=0.1
    tau_2=0.1
    CR_l=0.05
    CR_u=1
    F_in=0.5
    CR_in=0.8
    
    # Split data in batches for efficient RAM usage
    split_bool=0    # 1 for split / 0 for no split
    split_size=1000
    
    key=['databaseFile','nGM', 'comp_num', 'comp_idx', 'Tmin', 'Tmax', 'minScale', 'maxScale', 'a_g', 'zeta', 'soiltype', 'sameEvent','max_diff_max','max_diff_min','nSeed','w_CF','sameEventStation']
    value=[databaseFile,nGM, comp_num, comp_idx, Tmin, Tmax, minScale, maxScale, a_g, zeta, soiltype,sameEvent,max_diff_max,max_diff_min,nSeed,w_CF,sameEventStation]
    temp=list(zip(key, value))
    selectionParams=dict(temp)
    
    # delete variables
    for x in key:
        exec("del(%s)" % (x))
    
    key=['Vs30', 'Mag', 'D']
    value=[Vs30, Mag, D]
    temp=list(zip(key, value))
    allowedRecs=dict(temp)
    
    # delete variables
    for x in key:
        exec("del(%s)" % (x))
    
    key=['MaxGen', 'nPop','F_l','F_u','tau_1','tau_2','CR_l','CR_u','F_in','CR_in']
    value=[MaxGen, nPop,F_l, F_u,tau_1,tau_2,CR_l,CR_u,F_in,CR_in]
    temp=list(zip(key, value))
    DE_par=dict(temp)
    
    # delete variables
    for x in key:
        exec("del(%s)" % (x))
    
    key=['split_bool', 'split_size']
    value=[split_bool, split_size]
    temp=list(zip(key, value))
    split_data=dict(temp)
    
    # delete variables
    for x in key:
        exec("del(%s)" % (x))
    del key,value,temp,x
    return selectionParams,allowedRecs,DE_par,split_data