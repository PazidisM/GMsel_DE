# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 15:56:04 2018

@author: Marios D. Pazidis
"""

import numpy as np
import numba

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

#def fast_non_dominated_sorting(Cost_1,Cost_2,cond):
#    nPop_n=len(Cost_1)
#    F={}
#    F[1]=[]
#    S={}
#    n={}
#    for p in range(nPop_n):
#        S[p]=[]
#        n[p]=0
#        for q in range(nPop_n):
#            if p==q:
#                continue
#            case=two_obj_dominance(Cost_1[p],Cost_2[p],Cost_1[q],Cost_2[q],cond)
#            if case==0:
#                S[p].append(q)
#            elif case==1:
#                n[p]=n[p]+1
#        if n[p]==0:
#            F[1].append(p)
#    
#    i=1
#    while F[i]:
#        Q=[]
#        for p in F[i]:
#            for q in S[p]:
#                n[q]=n[q]-1
#                if n[q]==0:
#                    Q.append(q)
#        i=i+1
#        F[i]=Q
#    del F[i]
#    return F



@numba.jit(nopython=True)
def fast_non_dominated_sorting(Cost_1,Cost_2,nPop):
    nPop_n=len(Cost_1)
    F=[]
    F.append(np.empty(0))
    S=[]
    n=[]
    for p in range(nPop_n):
        S.append(np.empty(0))
        n.append(np.zeros(1))
        for q in range(nPop_n):
            if p==q:
                continue
            case_01=False
            case_02=False
            if (Cost_1[int(q)]>Cost_1[int(p)]):
                if (Cost_2[int(p)]<=Cost_2[int(q)]):
                    case_01=True
            elif (Cost_1[int(p)]<=Cost_1[int(q)]):
                if (Cost_2[int(p)]<Cost_2[int(q)]):
                    case_01=True
            elif (Cost_1[int(p)]>Cost_1[int(q)]):
                if (Cost_2[int(p)]>=Cost_2[int(q)]):
                    case_02=True
            elif (Cost_1[int(p)]>=Cost_1[int(q)]):
                if (Cost_2[int(p)]>Cost_2[int(q)]):
                    case_02=True
                    
            if case_01:
                case=0    # p dominates q
            elif case_02:
                case=1    # q dominates p
            else:
                case=2    # nondominated
            
            if case==0:
                S[int(p)]=np.concatenate((S[int(p)],np.atleast_1d(np.array(int(q)))),0)
            elif case==1:
                n[int(p)]=n[int(p)]+1
        
        foo= ~np.any(n[int(p)])
        if foo:
            F[0]=np.concatenate((F[0],np.atleast_1d(np.array(int(p)))),0)
    if len(F[0])<nPop:
        i=0
        while F[i].size!=0:
            Q=np.empty(0)
            for p in F[i]:
                for q in S[int(p)]:
                    n[int(q)]=n[int(q)]-1
                    foo= ~np.any(n[int(q)])
                    if foo:
                        Q=np.concatenate((Q,np.atleast_1d(np.array(int(q)))),0)
            i=i+1
            F.append(Q)
        del F[i]
    return F


@numba.jit(nopython=True)
def fast_non_dominated_sorting_first_front(Cost_1,Cost_2,leng,nPop):
    
    Dom=np.zeros(leng)
    p=0
    while p < leng-1:
        if Dom[p]==1:
            p=p+1
            continue
        q=p+1
        while q<leng:
            if p==q or Dom[q]==1:
                q=q+1
                continue
            case_01=False
            case_02=False
            if (Cost_1[int(q)]>Cost_1[int(p)]):
                if (Cost_2[int(p)]<=Cost_2[int(q)]):
                    case_01=True
            elif (Cost_1[int(p)]<=Cost_1[int(q)]):
                if (Cost_2[int(p)]<Cost_2[int(q)]):
                    case_01=True
            elif (Cost_1[int(p)]>Cost_1[int(q)]):
                if (Cost_2[int(p)]>=Cost_2[int(q)]):
                    case_02=True
            elif (Cost_1[int(p)]>=Cost_1[int(q)]):
                if (Cost_2[int(p)]>Cost_2[int(q)]):
                    case_02=True
            
            if case_01:
                case=0    # p dominates q
            elif case_02:
                case=1    # q dominates p
            else:
                case=2    # nondominated
            
            if case==0:
                Dom[int(q)]=1
            elif case==1:
                Dom[int(p)]=1
            q=q+1
        p=p+1
        if p%100==0:
            print(p)
    foo= ~np.all(Dom[:])
    if foo:
        F_t=np.where(Dom==0)
        print(len(F_t[0]))
        l=len(F_t[0])
        Fr=np.zeros((l,2))
        foo=np.array(F_t)//nPop
        foo=foo.reshape((len(F_t[0]),1))
        Fr[:,0]=foo
#        Fr[:,1]=np.array(F_t)%nPop
#    else:
#        Fr=np.empty(0)
    F=Fr
    return F

foo=fast_non_dominated_sorting_first_front(P_CF_0,P_CF_1,len(P_CF_0),nPop)




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
            try:
                for i in range(1,len(F[f_idx])-1):
                    I['cd_0'][i]=(I['sort_CF_0'][i+1]-I['sort_CF_0'][i-1])/(max(I['sort_CF_0'])-min(I['sort_CF_0']))
                    I['cd_1'][i]=(I['sort_CF_1'][i+1]-I['sort_CF_1'][i-1])/(max(I['sort_CF_1'])-min(I['sort_CF_1']))
            except:
                print(F,f_idx,I['sort_CF_0'],I['sort_CF_1'])
            for i in range(len(F[f_idx])):
                I['cd'][i]=I['cd_0'][I['sort_CF_0_idx'].index(i)]+I['cd_1'][I['sort_CF_1_idx'].index(i)]
            I['sorted_final']=sorted(range(len(I['cd'])), key=lambda k: I['cd'][k])
            P_new=P_new+[F[f_idx][i] for i in I['sorted_final'][len(P_new)-nPop:]]
            nPop_new=len(P_new)
    return P_new