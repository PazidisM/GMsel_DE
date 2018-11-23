# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 15:56:04 2018

@author: Marios D. Pazidis
"""

import numpy as np

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

def fast_non_dominated_sorting(Cost_1,Cost_2,cond):
    nPop_n=len(Cost_1)
    F={}
    F[1]=[]
    S={}
    n={}
    for p in range(nPop_n):
        S[p]=[]
        n[p]=0
        for q in range(nPop_n):
            if p==q:
                continue
            case=two_obj_dominance(Cost_1[p],Cost_2[p],Cost_1[q],Cost_2[q],cond)
            if case==0:
                S[p].append(q)
            elif case==1:
                n[p]=n[p]+1
        if n[p]==0:
            F[1].append(p)
    
    i=1
    while F[i]:
        Q=[]
        for p in F[i]:
            for q in S[p]:
                n[q]=n[q]-1
                if n[q]==0:
                    Q.append(q)
        i=i+1
        F[i]=Q
    del F[i]
    return F


def crowding_distance_assignment(F,Cost_0,Cost_1,nPop):
    nPop_new=0
    f_idx=1
    P_new=[]
    while nPop_new!=nPop:
        nPop_trial=nPop_new+len(F[f_idx])
        if nPop_trial<=nPop:
            P_new=P_new+F[f_idx]
            nPop_new=nPop_trial
            f_idx=f_idx+1
        else:
            I={}
            I['sort_CF_0']=sorted(Cost_0[F[f_idx]])
            I['sort_CF_1']=sorted(Cost_1[F[f_idx]])
            I['sort_CF_0_idx']=sorted(range(len(Cost_0[F[f_idx]])), key=lambda k: Cost_0[F[f_idx]][k])
            I['sort_CF_1_idx']=sorted(range(len(Cost_1[F[f_idx]])), key=lambda k: Cost_1[F[f_idx]][k])
            #I['cd_0']=np.zeros((len(F[f_idx])))
            #I['cd_1']=np.zeros((len(F[f_idx])))
            I['cd_0']=[0]*len(F[f_idx])
            I['cd_1']=[0]*len(F[f_idx])
            I['cd_0'][0]=np.inf
            I['cd_0'][-1]=np.inf
            I['cd_1'][0]=np.inf
            I['cd_1'][-1]=np.inf
            I['cd']=[0]*len(F[f_idx])
            for i in range(1,len(F[f_idx])-1):
                I['cd_0'][i]=(I['sort_CF_0'][i+1]-I['sort_CF_0'][i-1])/(max(I['sort_CF_0'])-min(I['sort_CF_0']))
                I['cd_1'][i]=(I['sort_CF_1'][i+1]-I['sort_CF_1'][i-1])/(max(I['sort_CF_1'])-min(I['sort_CF_1']))
            for i in range(len(F[f_idx])):
                I['cd'][i]=I['cd_0'][I['sort_CF_0_idx'].index(i)]+I['cd_1'][I['sort_CF_1_idx'].index(i)]
            I['sorted_final']=sorted(range(len(I['cd'])), key=lambda k: I['cd'][k])
            P_new=P_new+[F[f_idx][i] for i in I['sorted_final'][len(P_new)-nPop:]]
            nPop_new=len(P_new)
    return P_new