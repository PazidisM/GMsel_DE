# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 18:53:38 2018

@author: Marios D. Pazidis
"""

import math
import numpy as np
from Cost_functions import Obj_RMSE as Obj_RMSE
from Cost_functions import Obj_std as Obj_std
from tqdm import tqdm

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
        print('Initializing Batch '+str(index_spl)+' of '+str(split_data['split_num']))
        
        fileName=folders['Combinations']+'\Combs_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Combs_split=np.loadtxt(fileName).astype(int)
        fileName=folders['Sa_unsc_ave']+'\Sa_unsc_ave_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Sa_unsc_ave_split=np.loadtxt(fileName)
        Par_F=np.matlib.repmat(F_in,len(Combs_split),nPop)
        Par_CR=np.matlib.repmat(CR_in,len(Combs_split),nPop)
        Cost_Obj_01=np.empty([len(Combs_split),nPop])
        Cost_Obj_02=np.empty([len(Combs_split),nPop])
        
        #for cmb in range(len(Combs_split)):
        for cmb in tqdm(range(len(Combs_split)),miniters =round(len(Combs_split)*0.05),desc='% of batch'):
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
            sf_mean_pop=np.mean(sf,axis=0)
            
            for x in range(nPop):
                Cost_Obj_01[cmb,x]=Obj_RMSE(Sa_unsc_ave_split[cmb],Sa_Tgt,sf_mean_pop[x])
                Cost_Obj_02[cmb,x]=Obj_std(Sa_suite,sf[:,x],Sa_unsc_ave_suite,sf_mean_pop[x])
            
        fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Par_F,formats['fmt_Par_F_CR'])
        fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Par_CR,formats['fmt_Par_F_CR'])
        fileName=folders['Cost_Obj_01']+'\Cost_Obj_01_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Cost_Obj_01,formats['fmt_Cost'])
        fileName=folders['Cost_Obj_02']+'\Cost_Obj_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        np.savetxt(fileName, Cost_Obj_02,formats['fmt_Cost'])        
            
        index_spl=index_spl+1
        index=split_end

