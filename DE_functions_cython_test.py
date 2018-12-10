# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 14:40:28 2018

@author: Marios
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 18:53:38 2018

@author: Marios D. Pazidis
"""
import numpy as np
import random
from CDA import crowding_distance_assignment as crowding_distance_assignment
from TOD import two_obj_dominance as two_obj_dominance
from NSGAII_functions import fast_non_dominated_sorting as fast_non_dominated_sorting
from jDE_main import jDE_main as jDE_main


import math
import numba
import Cost_functions 
from tqdm import tqdm
#import math
#import time




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
            
            P=jDE_main(cmb,Sa_Tgt,formats,folders,Max_sf_ind,Min_sf_ind,nGM,nPop,F_l,F_u,tau_1,tau_2,cond,MaxGen,index,P,comb,w,Sa_unsc_ave_suite,Sa_suite)
        
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
        
        
        