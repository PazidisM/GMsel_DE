# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 11:54:11 2018

@author: Marios
"""
import numpy as np
import numba
@numba.jit(nopython=True, parallel=True)
def crowding_distance_assignment(F,Cost_0,Cost_1,nPop):
    nPop_new=0
    f_idx=0
    foo=F[f_idx]
    P_new=[]
#    P_new.append(np.empty(0))
    while nPop_new!=nPop:
        foo=F[f_idx]
        nPop_trial=nPop_new+len(foo)
        if nPop_trial<=nPop:
            for item in foo:
                P_new.append(item)
            nPop_new=nPop_trial
            f_idx=f_idx+1
        else:
            C_0=np.empty(len(foo))
            C_1=np.empty(len(foo))
            for j in range(len(foo)):
                C_0[j]=Cost_0[int([i for i in foo][j])]
                C_1[j]=Cost_1[int([i for i in foo][j])]
#            I_sort_CF_0=sorted(C_0)
#            I_sort_CF_1=sorted(C_1)
            I_sort_CF_0=np.sort(C_0)
            I_sort_CF_1=np.sort(C_1)
#            I_sort_CF_0_idx=sorted(range(len(C_0)), key=lambda k: C_0[k])
#            I_sort_CF_1_idx=sorted(range(len(C_1)), key=lambda k: C_1[k])
            I_sort_CF_0_idx=np.argsort(C_0)
            I_sort_CF_1_idx=np.argsort(C_1)
            #I_cd_0=np.zeros((len(F[f_idx])))
            #I['cd_1']=np.zeros((len(F[f_idx])))
            I_cd_0=np.zeros(len(foo))
            I_cd_1=np.zeros(len(foo))
            I_cd_0[0]=10e16
            I_cd_0[-1]=10e16
            I_cd_1[0]=10e16
            I_cd_1[-1]=10e16
            I_cd=np.zeros(len(foo))
            rang_0=(np.max(I_sort_CF_0)-np.min(I_sort_CF_0))
            rang_1=(np.max(I_sort_CF_1)-np.min(I_sort_CF_1))
            if rang_0!=0 and rang_1!=0:
                for i in range(1,len(foo)-1):
                    I_cd_0[i]=(I_sort_CF_0[i+1]-I_sort_CF_0[i-1])/rang_0
                    I_cd_1[i]=(I_sort_CF_1[i+1]-I_sort_CF_1[i-1])/rang_1
            elif rang_0!=0:
                for i in range(1,len(foo)-1):
                    I_cd_0[i]=(I_sort_CF_0[i+1]-I_sort_CF_0[i-1])/rang_0
                    I_cd_1[i]=0
            elif rang_1!=0:
                for i in range(1,len(foo)-1):
                    I_cd_0[i]=0
                    I_cd_1[i]=(I_sort_CF_0[i+1]-I_sort_CF_0[i-1])/rang_1
            else:
                for i in range(1,len(foo)-1):
                    I_cd_0[i]=0
                    I_cd_1[i]=0
            for i in range(len(foo)):
                f0=np.where(I_sort_CF_0_idx==i)[0][0]
                f1=np.where(I_sort_CF_1_idx==i)[0][0]
                I_cd[i]=I_cd_0[f0]+I_cd_1[f1]
            
            I_sorted_final=np.argsort(I_cd)
            #P_new=P_new+[F[f_idx][i] for i in I_sorted_final[len(P_new)-nPop:]]
#            P_new.append([F[f_idx][i] for i in I_sorted_final[len(P_new)-nPop:]])
            for i in I_sorted_final[len(P_new)-nPop:]:
                #for j in P_new
                P_new.append(foo[i])
                #print(i)
            #P_new.append([F[f_idx][i] for i in I_sorted_final[len(P_new)-nPop:]])
            
            nPop_new=len(P_new)
    return P_new

ans=[F[f_idx][i] for i in I_sorted_final[len(P_new)-nPop:]]
crowding_distance_assignment(F,Cost_0,Cost_1,nPop)

