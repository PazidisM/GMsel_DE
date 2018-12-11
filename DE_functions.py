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
from CDA import crowding_distance_assignment as crowding_distance_assignment
from TOD import two_obj_dominance as two_obj_dominance
from Cost_functions import CF_0
from Cost_functions import CF_1

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
        
        for cmb in tqdm(range(split_end-split_start),miniters =round((split_end-split_start)*0.10),desc='% of batch'):
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
        
        for cmb in tqdm(range(split_end-split_start),miniters =int((split_end-split_start)*0.05),desc='% of batch',maxinterval=3600):
        #for cmb in range(split_end-split_start):
            comb=Combs_split[cmb]
            fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
            w=np.loadtxt(fileName)
            Sa_unsc_ave_suite=Sa_unsc_ave_split[cmb,:]
            Sa_suite=Sa[comb,:]
            
            for gen in range(MaxGen):
            
                C={}
                C['w']=np.zeros((nGM, nPop))
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
#                    Xr=random.sample([i for i in range(0,nPop) if i != x],3)
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
                    u=np.repeat(u, 56,axis=1)
                    # Domination check
#                    u_c_01=(np.sum((Sa_unsc_ave_suite+np.sum(u)/len(u)-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5
#                    u_c_02=(np.sum((Sa_suite+u-(Sa_unsc_ave_suite+np.sum(u)/len(u)))**2/(len(u)-1))/len(u))**0.5
#                    
                    u_c_01=CF_0(Sa_unsc_ave_suite,Sa_Tgt,u)
                    u_c_02=CF_1(Sa_suite,u,Sa_unsc_ave_suite)
                    
                    u=u[:,0]
                    u.shape=(nGM,1)
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
#                np.savetxt(fileName, R['w'][:,F[1]],formats['fmt_sf'])
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
#                np.savetxt(fileName, w,formats['fmt_sf'])
            
            fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
            np.savetxt(fileName, w,formats['fmt_sf'])
        
        fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['Par_F'],formats['fmt_Par_F_CR'])
        fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['Par_CR'],formats['fmt_Par_F_CR'])
        fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['CF_0'],formats['fmt_Cost'])
        fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, P['CF_1'],formats['fmt_Cost'])
        
        index_spl=index_spl+1
        index=split_end
        
        
        