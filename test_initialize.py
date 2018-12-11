# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 11:06:13 2018

@author: Marios
"""


import UIn
import Database_functions
import Split_functions
import folder_fmt_functions
import DE_functions
import DE_functions_cython_test
#start_time=time()

SaveFolder='C:\Marios\Research\Tests\GM_selection_DEMO\Python\GMsel_DE\TEST'
folders=folder_fmt_functions.folder_init(SaveFolder)
selectionParams,allowedRecs,DE_par,split_data=UIn.User_Input()
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
#'''
#



## Initialize populations
DE_functions.Initialization(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)
#end_time=time()
#dur=end_time-start_time
#print(dur)













def two_obj_dominance(p_obj_01,p_obj_02,q_obj_01,q_obj_02,cond):
    case_01=((p_obj_01<q_obj_01)&(p_obj_02<=q_obj_02))|((p_obj_01<=q_obj_01)&(p_obj_02<q_obj_02))
    case_02=((p_obj_01>q_obj_01)&(p_obj_02>=q_obj_02))|((p_obj_01>=q_obj_01)&(p_obj_02>q_obj_02))
    if cond=='min':
        if case_01:
            case=0    # p dominates q
        elif case_02:
            case=1    # q dominates p
        else:
            case=2    # nondominated
    elif cond=='max':
        if case_01:
            case=1    # q dominates p
        elif case_02:
            case=0    # p dominates q
        else:
            case=2    # nondominated
    return case

import math
import numba
import numpy as np
import Cost_functions 
from tqdm import tqdm
import random
#import math
#import time
#from NSGAII_functions import fast_non_dominated_sorting as fast_non_dominated_sorting
#from NSGAII_functions import crowding_distance_assignment as crowding_distance_assignment

split_size=split_data['split_size']
Max_sf_ind=selectionParams['Max_sf_ind']
Min_sf_ind=selectionParams['Min_sf_ind']
nGM=selectionParams['nGM']
nPop=DE_par['nPop']
F_l=DE_par['F_l']
F_u=DE_par['F_u']
tau_1=DE_par['tau_1']
tau_2=DE_par['tau_2']
cond='min'
MaxGen=DE_par['MaxGen']

index=0
index_spl=0

split_start=index%split_size+index//split_size*split_size
split_end=split_start+split_size
if split_end>NSeed:
    split_end=NSeed


P={}
fileName=folders['Combinations']+'\Combs_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
Combs_split=np.loadtxt(fileName).astype(int)
fileName=folders['Sa_unsc_ave']+'\Sa_unsc_ave_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
Sa_unsc_ave_split=np.loadtxt(fileName)
fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['Par_F']=np.loadtxt(fileName)
fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['Par_CR']=np.loadtxt(fileName)
fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['CF_0']=np.loadtxt(fileName)
fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['CF_1']=np.loadtxt(fileName)    

cmb=0

comb=Combs_split[cmb]
fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
w=np.loadtxt(fileName)
Sa_unsc_ave_suite=Sa_unsc_ave_split[cmb,:]
Sa_suite=Sa[comb,:]

gen=0

C={}
C['w']=np.zeros((nGM, nPop))
C['CF_0']=np.full((1, nPop),np.inf)
C['CF_1']=np.full((1, nPop),np.inf)
C['Par_F']=np.full((1, nPop),np.inf)
C['Par_CR']=np.full((1, nPop),np.inf)

print('\n')
print('Optimizing Batch '+str(index_spl)+' of '+str(split_data['split_num']))

split_start=index%split_size+index//split_size*split_size
split_end=split_start+split_size
if split_end>NSeed:
    split_end=NSeed

P={}
fileName=folders['Combinations']+'\Combs_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
Combs_split=np.loadtxt(fileName).astype(int)
fileName=folders['Sa_unsc_ave']+'\Sa_unsc_ave_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
Sa_unsc_ave_split=np.loadtxt(fileName)
fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['Par_F']=np.loadtxt(fileName)
fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['Par_CR']=np.loadtxt(fileName)
fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['CF_0']=np.loadtxt(fileName)
fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
P['CF_1']=np.loadtxt(fileName)

P_CF_0=P['CF_0']
P_CF_1=P['CF_1']
P_Par_F=P['Par_F']
P_Par_CR=P['Par_CR']
formats_fmt_sf=formats['fmt_sf']
formats_fill_fn_all=formats['fill_fn_all']
folders_Scaling_factors=folders['Scaling_factors']









for x in range(nPop):
    rand_1=random.uniform(0,1)
    rand_2=random.uniform(0,1)
    rand_3=random.uniform(0,1)
    rand_4=random.uniform(0,1)
    
    ## Mutation (DE/rand/1 for now) ##
    foo=list(range(0,nPop))
    foo.remove(x)
    Xr=random.sample(foo,3)
    #Xr=np.random.choice(foo,3)
    if rand_2<tau_1:
        u_F=F_l+rand_1*F_u
    else:
        u_F=P['Par_F'][cmb,x]
    
    u=w[:,Xr[0]]+u_F*(w[:,Xr[1]]-w[:,Xr[2]])   #DE/rand/1
    
    # enforce constraints in individual values
    case=u>Max_sf_ind[comb]
    u[case]=Max_sf_ind[comb[np.where(case)]]
    case=u<Min_sf_ind[comb]
    u[case]=Min_sf_ind[comb[np.where(case)]]
    
    ## Crossover ##
    if rand_4<tau_2:
        u_CR=rand_3
    else:
        u_CR=P['Par_CR'][cmb,x]
    
    #case_01=random.randint(0, len(u)-1)
    case_01=np.random.randint(0, len(u))
    for k in range(len(u)):
        case_02=random.uniform(0,1)<=u_CR
        #case_02=np.random.uniform(0,1)<=u_CR
        if not (case_01==k or case_02):
            u[k]=w[k,x]
    u.shape=(nGM,1)
    
    # Domination check
    u_c_01=(np.sum((Sa_unsc_ave_suite+np.sum(u)/len(u)-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5
    u_c_02=(np.sum((Sa_suite+u-(Sa_unsc_ave_suite+np.sum(u)/len(u)))**2/(len(u)-1))/len(u))**0.5
    
    case=two_obj_dominance(P['CF_0'][cmb,x],P['CF_1'][cmb,x],u_c_01,u_c_02,cond)
    
    if case==1:
        w[:,x]=u[:,0]     # DEMO
        P['CF_0'][cmb,x]=u_c_01
        P['CF_1'][cmb,x]=u_c_02
        P['Par_F'][cmb,x]=u_F
        P['Par_CR'][cmb,x]=u_CR
    elif case==2:
        C['w'][:,x]=u[:,0]
        C['CF_0'][0,x]=u_c_01
        C['CF_1'][0,x]=u_c_02
        C['Par_F'][0,x]=u_F
        C['Par_CR'][0,x]=u_CR
    
# combine parent and offspring populations
R={}
case=C['CF_0']!=np.inf
for key in [ v for v in C if v != 'w' ]:
    R[key]=np.concatenate((P[key][cmb,:],C[key][case]),axis=0)
R['w']=np.concatenate((w,C['w'][:,np.where(case)[1]]),axis=1)

# NSGAII - sort and truncate
F=fast_non_dominated_sorting(R['CF_0'],R['CF_1'],nPop,cond)
F2={}
for i in range(len(F)):
    F2[i+1]=list(F[i].astype(int))
F=F2
P_n=crowding_distance_assignment(F,R['CF_0'],R['CF_1'],nPop)

# form next generation population
for key in [ v for v in C if v != 'w' ]:
    P[key][cmb,:]=R[key][P_n]
w=R['w'][:,P_n]






