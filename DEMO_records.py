# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 22:08:10 2018

@author: Marios D. Pazidis
"""

#%reset -f
from time import time
## User Input
# Selection parameters
#import os
#import sys
#from math import log as log
#import numpy as np
#import scipy
#import scipy.io
#import itertools
import Database_functions
#import matplotlib
#import matplotlib.pyplot as plt
import Split_functions
import folder_fmt_functions
import DE_functions
Inf=1.8e308
'''
#################################   User Input    #################################
'''
start_time=time()
# selectionParams
selectionParams={}
SaveFolder='C:\Marios\Research\Tests\GM_selection_DEMO\Python\GMsel_DE\TEST'
folders=folder_fmt_functions.folder_init(SaveFolder)

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
sameEventStation=1     # 0: can only contain 1 component per suite, 1: can contain both components
max_diff_max=6
max_diff_min=0.25
nSeed=3


# allowedRecs
Vs30=[1000, Inf]
Mag=[5, Inf]
D=[10, Inf]

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
split_bool=1    # 1 for split / 0 for no split
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

'''
#################################   Database screening / Combinations    #################################
'''

## screen database & define target spectrum
Ndatabase, Rec_db_metadata, Sa, Periods=Database_functions.screen_database(selectionParams,allowedRecs)
Sa_Tgt=Database_functions.EC8_Elastic_Spectrum_Type_1(selectionParams['soiltype'],selectionParams['a_g'],selectionParams['zeta'],Periods['T_all'])

## Define globals for simplicity / truncate arrays
T=Periods['T_match']
Sa_Tgt=Sa_Tgt[:,Periods['idx_T_match'][0]]
Sa=Sa[:,Periods['idx_T_match'][0]]

## Calculate max/min scaling factors and further screen the database
Sa,Rec_db_metadata,sf_ind, selectionParams, Ndatabase=Database_functions.Ind_sc_factor(Sa,Rec_db_metadata,Sa_Tgt,selectionParams)

## Calculate combinations / Create GM suites
Combs, Sa_unsc_ave, NSeed, split_data=Database_functions.Combinations(selectionParams,Rec_db_metadata,Ndatabase,Sa_Tgt,Sa,split_data,sf_ind)

## Printing formats
formats=folder_fmt_functions.fmt(SaveFolder,split_data,NSeed,Rec_db_metadata,DE_par)

## Save data to files
Split_functions.Pre_run_split(split_data,NSeed,Combs,Sa_unsc_ave,folders,formats)
del Combs, Sa_unsc_ave

'''
#################################   DE    #################################
'''

## Initialize populations
DE_functions.Initialization(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)

## Differential Evolution
DE_functions.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)
end_time=time()
dur=end_time-start_time
print(dur)