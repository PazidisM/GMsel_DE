# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 18:53:38 2018

@author: Marios D. Pazidis
"""

import math
import numpy as np
import Cost_functions 
from tqdm import tqdm
import random
#import math
#import time
from NSGAII_functions import fast_non_dominated_sorting as fast_non_dominated_sorting
from NSGAII_functions import crowding_distance_assignment as crowding_distance_assignment


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


def Initialization(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa):
    
    split_size=split_data['split_size']
    
    minSF=math.log(selectionParams['minScale'])
    maxSF=math.log(selectionParams['maxScale'])
    Max_sf_ind=selectionParams['Max_sf_ind']
    Min_sf_ind=selectionParams['Min_sf_ind']
    nGM=selectionParams['nGM']
    nPop=DE_par['nPop']
    F_in=DE_par['F_in']
    CR_in=DE_par['CR_in']
    
    index=0
    index_spl=0
    while index<NSeed:
        split_start=index%split_size+index//split_size*split_size
        split_end=split_start+split_size
        if split_end>NSeed:
            split_end=NSeed
        
        print('\n')
        print('Initializing Batch '+str(index_spl+1)+' of '+str(split_data['split_num']+1))
        
        fileName=folders['Combinations']+'\Combs_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Combs_split=np.loadtxt(fileName).astype(int)
        fileName=folders['Sa_unsc_ave']+'\Sa_unsc_ave_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Sa_unsc_ave_split=np.loadtxt(fileName)
        Par_F=np.matlib.repmat(F_in,len(Combs_split),nPop)
        Par_CR=np.matlib.repmat(CR_in,len(Combs_split),nPop)
        CF_0=np.empty([len(Combs_split),nPop])
        CF_1=np.empty([len(Combs_split),nPop])
        
        #for cmb in range(len(Combs_split)):
        for cmb in tqdm(range(len(Combs_split)),miniters =round(len(Combs_split)*0.10),desc='% of batch'):
            sf=np.random.uniform(minSF,maxSF,(nGM, nPop))
            Sa_suite=Sa[Combs_split[cmb],:]
            Sa_unsc_ave_suite=Sa_unsc_ave_split[cmb,:]
            
            # Enforce sf limits
            Max_sf=Max_sf_ind[Combs_split[cmb]]
            Min_sf=Min_sf_ind[Combs_split[cmb]]
            Max_sf.shape=(nGM,1)
            Min_sf.shape=(nGM,1)
            Max_sf=np.matlib.repmat(Max_sf,1,nPop)
            Min_sf=np.matlib.repmat(Min_sf,1,nPop)
            case=sf>Max_sf
            sf[case]=Max_sf[case]
            case=sf<Min_sf
            sf[case]=Min_sf[case]
            fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
            np.savetxt(fileName, sf,formats['fmt_sf'])
            
            for x in range(nPop):
                CF_0[cmb,x]=Cost_functions.CF_0(Sa_unsc_ave_split[cmb],Sa_Tgt,sf[:,x])
                CF_1[cmb,x]=Cost_functions.CF_1(Sa_suite,sf[:,x],Sa_unsc_ave_suite)
            
        fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Par_F,formats['fmt_Par_F_CR'])
        fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Par_CR,formats['fmt_Par_F_CR'])
        fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, CF_0,formats['fmt_Cost'])
        fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, CF_1,formats['fmt_Cost'])        
        
        index_spl=index_spl+1
        index=split_end

def jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa):

    split_size=split_data['split_size']
#    minSF=math.log(selectionParams['minScale'])
#    maxSF=math.log(selectionParams['maxScale'])
    Max_sf_ind=selectionParams['Max_sf_ind']
    Min_sf_ind=selectionParams['Min_sf_ind']
    nGM=selectionParams['nGM']
    nPop=DE_par['nPop']
    F_l=DE_par['F_l']
    F_u=DE_par['F_u']
#    CR_l=DE_par['CR_l']
#    CR_u=DE_par['CR_u']
    tau_1=DE_par['tau_1']
    tau_2=DE_par['tau_2']
    cond='min'
    MaxGen=DE_par['MaxGen']
    
    index=0
    index_spl=0
    while index<NSeed:
        
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
        
        for cmb in tqdm(range(len(Combs_split)),miniters =int(len(Combs_split)*0.05),desc='% of batch'):
            fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
            P['sf_cmb']=np.loadtxt(fileName)
            Sa_unsc_ave_suite=Sa_unsc_ave_split[cmb,:]
            Sa_suite=Sa[Combs_split[cmb],:]
            
            for gen in range(MaxGen):
            
                C={}
                C['sf_cmb']=np.zeros((nGM, nPop))
                C['CF_0']=np.full((1, nPop),np.inf)
                C['CF_1']=np.full((1, nPop),np.inf)
                C['Par_F']=np.full((1, nPop),np.inf)
                C['Par_CR']=np.full((1, nPop),np.inf)
                
                for x in range(nPop):
                    rand_1=random.uniform(0,1)
                    rand_2=random.uniform(0,1)
                    rand_3=random.uniform(0,1)
                    rand_4=random.uniform(0,1)
                    
                    ## Mutation (DE/rand/1 for now) ##
                    Xr=random.sample([i for i in range(0,nPop) if i != x],3)
                    if rand_2<tau_1:
                        u_F=F_l+rand_1*F_u
                    else:
                        u_F=P['Par_F'][cmb,x]
                    
                    u=P['sf_cmb'][:,Xr[0]]+u_F*(P['sf_cmb'][:,Xr[1]]-P['sf_cmb'][:,Xr[2]])   #DE/rand/1
                    
                    # enforce constraints in individual values
                    case=u>Max_sf_ind[Combs_split[cmb]]
                    u[case]=Max_sf_ind[Combs_split[cmb,np.where(case)][0]]
                    case=u<Min_sf_ind[Combs_split[cmb]]
                    u[case]=Min_sf_ind[Combs_split[cmb,np.where(case)][0]]
                    
                    ## Crossover ##
                    if rand_4<tau_2:
                        u_CR=rand_3
                    else:
                        u_CR=P['Par_CR'][cmb,x]
                    
                    case_01=random.choice(u)==u
                    case_02=np.array([random.random() for i in range(nGM)])<=u_CR
                    case=~np.array(case_01|case_02)
                    u[case]=P['sf_cmb'][case,x]
                    u.shape=(nGM,1)
                    
                    # Domination check
                    u_c_01=Cost_functions.CF_0(Sa_unsc_ave_suite,Sa_Tgt,u)
                    u_c_02=Cost_functions.CF_1(Sa_suite,u,Sa_unsc_ave_suite)
                    
                    case=two_obj_dominance(P['CF_0'][cmb,x],P['CF_1'][cmb,x],u_c_01,u_c_02,cond)
                    
                    if case==1:
                        P['sf_cmb'][:,x]=u[:,0]     # DEMO
                        P['CF_0'][cmb,x]=u_c_01
                        P['CF_1'][cmb,x]=u_c_02
                        P['Par_F'][cmb,x]=u_F
                        P['Par_CR'][cmb,x]=u_CR
                    elif case==2:
                        C['sf_cmb'][:,x]=u[:,0]
                        C['CF_0'][0,x]=u_c_01
                        C['CF_1'][0,x]=u_c_02
                        C['Par_F'][0,x]=u_F
                        C['Par_CR'][0,x]=u_CR
                    
                # combine parent and offspring populations
                R={}
                case=C['CF_0']!=np.inf
                for key in [ v for v in C if v != 'sf_cmb' ]:
                    R[key]=np.concatenate((P[key][cmb,:],C[key][case]),axis=0)
                R['sf_cmb']=np.concatenate((P['sf_cmb'],C['sf_cmb'][:,np.where(case)[1]]),axis=1)
                
                # NSGAII - sort & truncate
                F=fast_non_dominated_sorting(R['CF_0'],R['CF_1'],cond)
                P_n=crowding_distance_assignment(F,R['CF_0'],R['CF_1'],nPop)
                
                # form next generation population
                for key in [ v for v in C if v != 'sf_cmb' ]:
                    P[key][cmb,:]=R[key][P_n]
                P['sf_cmb']=R['sf_cmb'][:,P_n]
            
            #print(len(F[1]))
            # 1st Pareto front
#            if len(F[1])<nPop:
#                fileName=folders['Par_F_Pareto']+'\Par_F_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, R['Par_F'][F[1]],formats['fmt_Par_F_CR'])
#                fileName=folders['Par_CR_Pareto']+'\Par_CR_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, R['Par_CR'][F[1]],formats['fmt_Par_F_CR'])
#                fileName=folders['CF_0_Pareto']+'\CF_0_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, R['CF_0'][F[1]],formats['fmt_Cost'])
#                fileName=folders['CF_1_Pareto']+'\CF_1_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, R['CF_1'][F[1]],formats['fmt_Cost'])
#                fileName=folders['Scaling_factors_Pareto']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, R['sf_cmb'][:,F[1]],formats['fmt_sf'])
#            else:
#                fileName=folders['Par_F_Pareto']+'\Par_F_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, P['Par_F'][cmb,:],formats['fmt_Par_F_CR'])
#                fileName=folders['Par_CR_Pareto']+'\Par_CR_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, P['Par_CR'][cmb,:],formats['fmt_Par_F_CR'])
#                fileName=folders['CF_0_Pareto']+'\CF_0_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, P['CF_0'][cmb,:],formats['fmt_Cost'])
#                fileName=folders['CF_1_Pareto']+'\CF_1_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, P['CF_1'][cmb,:],formats['fmt_Cost'])
#                fileName=folders['Scaling_factors_Pareto']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
#                np.savetxt(fileName, P['sf_cmb'],formats['fmt_sf'])
#            
            fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
            np.savetxt(fileName, P['sf_cmb'],formats['fmt_sf'])
        
        fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['Par_F'],formats['fmt_Par_F_CR'])
        fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['Par_CR'],formats['fmt_Par_F_CR'])
        fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['CF_0'],formats['fmt_Cost'])
        fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['CF_1'],formats['fmt_Cost'])
