# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 13:05:29 2018

@author: Marios D. Pazidis
"""

import numpy as np
import scipy
import scipy.io
import itertools
import math
import sys
from tqdm import tqdm


"""
DEFINITIONS
"""

def dictionary_fn(key):
    value=[]
    for x in key:
        print(eval(x))
        #value.append(eval(x))
    temp=list(zip(key, value))
    z=dict(temp)
    return z


def screen_database(selectionParams,allowedRecs):

    # load database
    mat = scipy.io.loadmat(selectionParams['databaseFile'])

    Rec_db_metadata={}

    if selectionParams['comp_num']==2:

        ## Metadata
        Event_name=np.tile(mat['dirLocation'],[2,1])
        EQID=np.matlib.repmat(mat['EQID'],2,1)
        EQidx=np.transpose(np.matlib.repmat(np.array(range(len(mat['EQID']))),1,2))
        Distance=np.matlib.repmat(mat['closest_D'],2,1)
        Vs30=np.matlib.repmat(mat['soil_Vs30'],2,1)
        Magnitude=np.matlib.repmat(mat['magnitude'],2,1)
        Component=np.array(list(itertools.repeat(1,len(mat['EQID'])))+list(itertools.repeat(2,len(mat['EQID']))))

        ## Spectral Accelerations
        Sa_1=mat['Sa_1']
        Sa_2=mat['Sa_2']
        Sa=np.concatenate((Sa_1,Sa_2), axis=0)

    else:

        ## Metadata
        Event_name=mat['dirLocation']
        EQID=mat['EQID']
        EQidx=np.array(range(len(mat['EQID'])))
        Distance=mat['closest_D']
        Vs30=mat['soil_Vs30']
        Magnitude=mat['magnitude']
        Component=np.array(list(itertools.repeat(1,len(mat['EQID']))))

        if selectionParams['comp_idx']==1:
            ## Spectral Accelerations
            Sa=mat['Sa_1']
        else:
            ## Spectral Accelerations
            Sa=mat['Sa_2']

    Sa=np.log(Sa)

    key=['EQID','EQidx','Event_name', 'Distance', 'Vs30', 'Magnitude', 'Component']
    value=[EQID,EQidx, Event_name, Distance, Vs30, Magnitude, Component]
    temp=list(zip(key, value))
    Rec_db_metadata=dict(temp)

    # Periods
    idx_T_all=np.where(mat['Periods']<=10)[1]    # throw out periods > 10s
    idx_T_all=idx_T_all.reshape((1,len(idx_T_all)))
    T_all=np.array(mat['Periods'])[0][idx_T_all]

    idx_T_match=np.intersect1d(np.where(T_all>=selectionParams['Tmin'])[1],np.where(T_all<=selectionParams['Tmax'])[1])
    idx_T_match=idx_T_match.reshape((1,len(idx_T_match)))
    T_match=T_all[0][idx_T_match]

    key=['idx_T_all','T_all','idx_T_match','T_match']
    value=[idx_T_all,T_all,idx_T_match,T_match]
    temp=list(zip(key, value))
    Periods=dict(temp)

    # Screening
    Vs_screen=np.intersect1d(np.where(Rec_db_metadata['Vs30']>=allowedRecs['Vs30'][0])[0],np.where(Rec_db_metadata['Vs30']<=allowedRecs['Vs30'][1])[0])
    Magnitude_screen=np.intersect1d(np.where(Rec_db_metadata['Magnitude']>=allowedRecs['Mag'][0])[0],np.where(Rec_db_metadata['Magnitude']<=allowedRecs['Mag'][1])[0])
    Distance_screen=np.intersect1d(np.where(Rec_db_metadata['Distance']>=allowedRecs['D'][0])[0],np.where(Rec_db_metadata['Distance']<=allowedRecs['D'][1])[0])
    
    Fin_1_screen=np.intersect1d(Vs_screen,Magnitude_screen)
    Fin_screen=np.intersect1d(Fin_1_screen,Distance_screen)
    
    Sa=Sa[Fin_screen]
    
    for x in Rec_db_metadata:
        Rec_db_metadata[x]=Rec_db_metadata[x][Fin_screen]
    
    Ndatabase=len(Fin_screen)
    Rec_db_metadata['idx']=np.array(range(0,Ndatabase))
    
    print('\n')
    print('1st screening: Number of allowed GMs: ' + str(Ndatabase))
    print('1st screening: Number of unique events = '+ str(len(np.unique(Rec_db_metadata['EQID']))))
    print('\n')
    if (len(np.unique(Rec_db_metadata['EQID']))<selectionParams['nGM']) & (selectionParams['sameEvent']==1):
        sys.exit('Not enough Ground Motions - Required : '+str(selectionParams['nGM'])+' - Available : '+str(len(np.unique(Rec_db_metadata['EQID']))))

    return (Ndatabase, Rec_db_metadata, Sa, Periods)


def EC8_Elastic_Spectrum_Type_1(soiltype_EC8,a_g_EC8,zeta_EC8,periods_EC8):

    eta_EC8 =(10/(5+zeta_EC8*100))**0.5
    
    if soiltype_EC8=='A':
        S=1
        TB=0.15
        TC=0.4
        TD=2.0
    elif soiltype_EC8=='B':
        S=1.2
        TB=0.15
        TC=0.5
        TD=2.0
    elif soiltype_EC8=='C':
        S=1.15
        TB=0.2
        TC=0.6
        TD=2.0
    elif soiltype_EC8=='D':
        S=1.35
        TB=0.2
        TC=0.8
        TD=2.0
    elif soiltype_EC8=='E':
        S=1.4
        TB=0.15
        TC=0.5
        TD=2.0
    
    Sa_Tgt=np.empty([1, periods_EC8.shape[1]])
    
    per_cond=np.array(np.where(periods_EC8<=TB))
    Sa_Tgt[0][per_cond[1]]=a_g_EC8*S*(1+periods_EC8[0][per_cond[1]]/TB*(eta_EC8*2.5-1))
    
    per_cond=np.array(np.where((periods_EC8<=TC)&(periods_EC8>TB)))
    Sa_Tgt[0][per_cond[1]]=a_g_EC8*S*eta_EC8*2.5
    
    per_cond=np.array(np.where((periods_EC8<=TD)&(periods_EC8>TC)))
    Sa_Tgt[0][per_cond[1]]=a_g_EC8*S*eta_EC8*2.5*TC/periods_EC8[0][per_cond[1]]
    
    per_cond=np.array(np.where(periods_EC8>TD))
    Sa_Tgt[0][per_cond[1]]=a_g_EC8*S*eta_EC8*2.5*TC*TD/periods_EC8[0][per_cond[1]]**2
    
    Sa_Tgt=np.log(Sa_Tgt)
    
    return Sa_Tgt


def Ind_sc_factor(Sa,Rec_db_metadata,Sa_Tgt,selectionParams):
    ## Individual scaling factor to match target spectrum
    sf_ind=np.sum(Sa_Tgt-Sa,1)/len(Sa_Tgt[0])
    
    ## Calclulate maximum individual scaling factors
    Max_sf_ind=selectionParams['max_diff_max']/np.exp(np.max(Sa-Sa_Tgt,1))
    Max_sf_ind[np.where(Max_sf_ind>selectionParams['maxScale'])]=selectionParams['maxScale']
    
    Min_sf_ind=selectionParams['max_diff_min']/np.exp(np.min(Sa-Sa_Tgt,1))
    Min_sf_ind[np.where(Min_sf_ind<selectionParams['minScale'])]=selectionParams['minScale']
    
    del_idx=np.concatenate((np.array(np.where(Max_sf_ind<selectionParams['minScale'])),np.array(np.where(Min_sf_ind>selectionParams['maxScale']))),axis=1)
    
    Max_sf_ind=np.delete(Max_sf_ind,del_idx,axis=0)
    Min_sf_ind=np.delete(Min_sf_ind,del_idx,axis=0)
    selectionParams['Max_sf_ind']=np.log(Max_sf_ind)
    selectionParams['Min_sf_ind']=np.log(Min_sf_ind)
    
    sf_ind=np.delete(sf_ind,del_idx,axis=0)
    Sa=np.delete(Sa,del_idx,axis=0)
    
    for x in Rec_db_metadata:
        Rec_db_metadata[x]=np.delete(Rec_db_metadata[x],del_idx,axis=0)
    
    Ndatabase=len(Sa)
    Rec_db_metadata['idx']=np.array(range(0,Ndatabase))
    
    print('2nd screening: Number of allowed GMs: ' + str(Ndatabase))
    print('2nd screening: Number of unique events = '+ str(len(np.unique(Rec_db_metadata['EQID']))))
    print('\n')
    if (len(np.unique(Rec_db_metadata['EQID']))<selectionParams['nGM']) & (selectionParams['sameEvent']==1):
        sys.exit('Not enough Ground Motions - Required : '+str(selectionParams['nGM'])+' - Available : '+str(len(np.unique(Rec_db_metadata['EQID']))))
    
    return  Sa,Rec_db_metadata,sf_ind, selectionParams, Ndatabase

def Combinations(selectionParams,Rec_db_metadata,Ndatabase,Sa_Tgt,Sa,split_data,sf_ind):

    case_st=bool(selectionParams['sameEventStation'])
    # Calculate combinations
    Combs=np.array(list(itertools.combinations(Rec_db_metadata['idx'], selectionParams['nSeed'])))
    
    # Remove same event combinations accordingly
    Unique_check=np.array([len(np.unique(Rec_db_metadata['EQID'][comb])) for comb in Combs])
    Combs=np.delete(Combs,np.where(Unique_check<=selectionParams['nSeed']-selectionParams['sameEvent']),axis=0)
    NSeed=len(Combs)
    
    # Incrementaly add records - Create trial sets
    rem=selectionParams['nGM']-selectionParams['nSeed']
    Sample_sel_idx=np.full((NSeed,rem),-1)
    #for ii in range(NSeed):
    for ii in tqdm(range(NSeed),miniters =round(NSeed*0.10),desc='Combinations'):
        for jj in range(rem):
            Sample_Cost=np.full((Ndatabase),np.inf)
            for kk in range(Ndatabase):
                Sample_add=np.append(kk,Sample_sel_idx[ii][np.where(Sample_sel_idx[ii]!=-1)])   # ignore -1 values from initialization
                Sample_suite=np.append(Combs[ii],Sample_add).astype(int)
                Sa_suite=Sa[Sample_suite,:]
                sf_ind_suite=sf_ind[Sample_suite]
                sf_ind_suite.shape=(len(sf_ind_suite),1)
                case_01=len(np.unique(Rec_db_metadata['EQID'][Sample_suite]))>len(Sample_suite)-selectionParams['sameEvent']   # same event records constraint
                case_02=len(np.unique(Rec_db_metadata['EQidx'][Sample_suite]))==len(Sample_suite)   # components allowed constraint
                if not case_02:
                    case_02=case_st
                
                if case_01 and case_02:
                    Sample_Sa=Sa[Sample_suite]
                    Sample_Sa_ave=np.mat(np.transpose(np.mean(Sample_Sa,axis=0)))
#                    Sample_sf=np.sum(Sa_Tgt-Sample_Sa_ave)/np.size(Sample_Sa_ave,1)
                    if selectionParams['w_CF'][1]==0:
#                        Sample_Cost[kk]=(np.sum((Sample_Sa_ave+np.sum(sf_ind_suite)/len(sf_ind_suite)-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5
                        Sample_Cost[kk]=(np.sum(np.power(Sample_Sa_ave+np.sum(sf_ind_suite)/len(sf_ind_suite)-Sa_Tgt,2))/np.size(Sa_Tgt,1))**0.5
#                        u_c_01=(np.sum((Sa_unsc_ave_suite+np.sum(u)/len(u)-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5
                    elif selectionParams['w_CF'][0]==0:
#                        Sample_Cost[kk]=(np.sum((Sa_suite+sf_ind_suite-(Sample_Sa_ave+np.sum(sf_ind_suite)/len(sf_ind_suite)))**2/(len(sf_ind_suite)-1))/len(sf_ind_suite))**0.5
                        Sample_Cost[kk]=(np.sum(np.power(Sa_suite+sf_ind_suite-(Sample_Sa_ave+np.sum(sf_ind_suite)/len(sf_ind_suite)),2)/(len(sf_ind_suite)-1))/len(sf_ind_suite))**0.5
                    else:
                        Sample_Cost[kk]=((np.sum(np.power(Sample_Sa_ave+np.sum(sf_ind_suite)/len(sf_ind_suite)-Sa_Tgt,2))/np.size(Sa_Tgt,1))+(np.sum(np.power(Sa_suite+sf_ind_suite-(Sample_Sa_ave+np.sum(sf_ind_suite)/len(sf_ind_suite)),2)/(len(sf_ind_suite)-1))/len(sf_ind_suite)))**0.5
            Sample_sel_idx[ii,jj]=Sample_Cost.argmin(0)
            
    Combs=np.concatenate((Combs,Sample_sel_idx),axis=1)
    Sa_combs=Sa[Combs.astype(int)]
    Sa_unsc_ave=np.mean(Sa_combs,axis=1)
    if not split_data['split_bool']:
        split_data['split_size']=NSeed
    split_data['split_num']=NSeed//(split_data['split_size']+1)
    print('\n')
    
    return Combs, Sa_unsc_ave, NSeed,split_data

