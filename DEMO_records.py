# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 22:08:10 2018

@author: Marios D. Pazidis
"""

%reset -f
## User Input
# Selection parameters
import os
import sys
from math import log as log
import numpy as np
import scipy
import scipy.io
import itertools
import Database_functions
import matplotlib
import matplotlib.pyplot as plt
os.chdir("C:\Marios\Research\Tests\GM_selection_DEMO\Python\GMsel_DE")
Inf=1.8e308
'''
###########User Input##########
'''
# selectionParams
selectionParams={}

databaseFile='ESD_Rexel_meta_data.mat'
#databaseFile='Resorce2013_meta_data.mat'
#databaseFile='NGA_W2_meta_data.mat'

nGM=10
comp_num=2
comp_idx=1    # used only if comp_num==1
Tmin=0.08
Tmax=2
minScale=0.5
maxScale=2
a_g=0.25
zeta=0.05
soiltype='A'
minScaleRec=log(0.15)
maxScaleRec=log(12)
std_w=0.0
sameEvent=1
max_diff_max=2
max_diff_min=0.25
nSeed=2


key=['databaseFile','nGM', 'comp_num', 'comp_idx', 'Tmin', 'Tmax', 'minScale', 'maxScale', 'a_g', 'zeta', 'soiltype', 'minScaleRec','maxScaleRec', 'std_w', 'sameEvent','max_diff_max','max_diff_min','nSeed']
value=[databaseFile,nGM, comp_num, comp_idx, Tmin, Tmax, minScale, maxScale, a_g, zeta, soiltype, minScaleRec, maxScaleRec,std_w,sameEvent,max_diff_max,max_diff_min,nSeed]
temp=list(zip(key, value))
selectionParams=dict(temp)

# delete variables
for x in key:
    exec("del(%s)" % (x))

# allowedRecs
Vs30=[1000,  Inf]
Mag=[5, Inf]
D=[10, Inf]

key=['Vs30', 'Mag', 'D']
value=[Vs30, Mag, D]
temp=list(zip(key, value))
allowedRecs=dict(temp)

# delete variables
for x in key:
    exec("del(%s)" % (x))

# DE parameters
MaxIt=100
nPop=50
beta_min=0.6
beta_max=0.6
pCR=0.8
VarSize=[1, selectionParams['nGM']]
VarMin=log(selectionParams['minScale'])
VarMax=log(selectionParams['maxScale'])

key=['MaxIt', 'nPop', 'beta_min', 'beta_max', 'pCR', 'VarSize', 'VarMin', 'VarMax']
value=[MaxIt, nPop, beta_min, beta_max, pCR, VarSize, VarMin, VarMax]
temp=list(zip(key, value))
DE_par=dict(temp)

# delete variables
for x in key:
    exec("del(%s)" % (x))
del key
del value
del temp
del x

## screen database & define target spectrum
Ndatabase, Rec_db_metadata, Sa, Periods=Database_functions.screen_database(selectionParams,allowedRecs)
Sa_Tgt=Database_functions.EC8_Elastic_Spectrum_Type_1(selectionParams['soiltype'],selectionParams['a_g'],selectionParams['zeta'],Periods['T_all'])

## Define globals for simplicity / truncate arrays
T=Periods['T_match']
Sa_Tgt=Sa_Tgt[:,Periods['idx_T_match'][0]]
Sa=Sa[:,Periods['idx_T_match'][0]]

## Calculate max/min scaling factors and further screen the database
Sa,Rec_db_metadata,sf_ind, Max_sf_ind, Min_sf_ind, Ndatabase=Database_functions.Ind_sc_factor(Sa,Rec_db_metadata,Sa_Tgt,selectionParams)

## Calculate combinations / Create GM suites
Combs, NSeed=Database_functions.Combinations(selectionParams,Rec_db_metadata)




















