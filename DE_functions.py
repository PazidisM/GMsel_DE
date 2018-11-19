# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 18:53:38 2018

@author: Marios D. Pazidis
"""
import math
import numpy as np
from Cost_functions import Obj_RMSE as Obj_RMSE
from tqdm import tqdm

def Initialization(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt):
    
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
        
        #for ii in range(len(Combs_split)):
        for ii in tqdm(range(len(Combs_split)),miniters =int(len(Combs_split)*0.05),desc='% of batch'):
            sf=np.random.uniform(minSF,maxSF,(nGM, nPop))
            Sample_sf=np.mean(sf,axis=0)
            
            # Enforce sf limits
            Max_sf=Max_sf_ind[Combs_split[ii]]
            Min_sf=Min_sf_ind[Combs_split[ii]]
            Max_sf.shape=(nGM,1)
            Min_sf.shape=(nGM,1)
            Max_sf=np.matlib.repmat(Max_sf,1,nPop)
            Min_sf=np.matlib.repmat(Min_sf,1,nPop)
            Max_sf_idx=np.where(sf>Max_sf)
            Min_sf_idx=np.where(sf<Min_sf)
            sf[Max_sf_idx]=Max_sf[Max_sf_idx]
            sf[Min_sf_idx]=Min_sf[Min_sf_idx]
            fileName=folders['Scaling_factors']+'\SF_'+str(index+ii).zfill(formats['fill_fn_all'])+'.out'
            np.savetxt(fileName, sf,formats['fmt_sf'])
            
            for jj in range(nPop):
                Cost_Obj_01[ii,jj]=Obj_RMSE(Sa_unsc_ave_split[ii],Sa_Tgt,Sample_sf[jj])
                Cost_Obj_02[ii,jj]=Obj_RMSE(Sa_unsc_ave_split[ii],Sa_Tgt,Sample_sf[jj])
            
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

