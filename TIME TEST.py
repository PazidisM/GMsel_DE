# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 10:08:23 2018

@author: Marios
"""

Cost_1=R['CF_0']
Cost_2=R['CF_1']

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



@numba.jit(nopython=True, parallel=True)
def fast_non_dominated_sorting(Cost_1,Cost_2,nPop,cond):
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


@numba.jit(nopython=True, parallel=True)
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


def two_obj_dominance(p_obj_01,p_obj_02,q_obj_01,q_obj_02,cond):
    case_01=((p_obj_01<q_obj_01)and(p_obj_02<=q_obj_02))or((p_obj_01<=q_obj_01)and(p_obj_02<q_obj_02))
    case_02=((p_obj_01>q_obj_01)and(p_obj_02>=q_obj_02))or((p_obj_01>=q_obj_01)and(p_obj_02>q_obj_02))
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


p_obj_01=R['CF_0'][0]
p_obj_02=R['CF_1'][0]

q_obj_01=R['CF_0'][1]
q_obj_02=R['CF_1'][1]





from NSGAII_functions import crowding_distance_assignment as crowding_distance_assignment

from CDA import crowding_distance_assignment as crowding_distance_assignment


for i in range(10):
    rep=1000
    tot_dur=0
    st=time.clock()
    for i in range(rep):
        #F= FNS.fast_non_dominated_sorting(Cost_1,Cost_2,cond)
        F= fast_non_dominated_sorting(Cost_1,Cost_2,nPop,cond)
#        P_n=crowding_distance_assignment(F,R['CF_0'],R['CF_1'],nPop)
    end=time.clock()
    dur=end-st
    print(dur)


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
case

p_obj_01=P['CF_0'][0,0]
p_obj_02=P['CF_1'][0,0]
q_obj_01=P['CF_0'][0,1]
q_obj_02=P['CF_1'][0,1]

case=two_obj_dominance(p_obj_01,p_obj_02,q_obj_01,q_obj_02,cond)
case




from DE_functions import jDE
%lprun -f jDE jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)




from NSGAII_functions import fast_non_dominated_sorting as fast_non_dominated_sorting


DE_functions.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)

%load_ext line_profiler

%lprun -f jDE_main jDE_main(cmb,Sa_Tgt,formats,folders,Max_sf_ind,Min_sf_ind,nGM,nPop,F_l,F_u,tau_1,tau_2,cond,MaxGen,index,P,comb,w,Sa_unsc_ave_suite,Sa_suite)


import line_profiler

jDE_main(cmb,Sa_Tgt,formats,folders,Max_sf_ind,Min_sf_ind,nGM,nPop,F_l,F_u,tau_1,tau_2,cond,MaxGen,index,P,comb,w,Sa_unsc_ave_suite,Sa_suite)

import  pyximport; pyximport.install()
from distutils.core import setup
from Cython.Build import cythonize

import FNS
from CDA import crowding_distance_assignment as crowding_distance_assignment
from TOD import two_obj_dominance as two_obj_dominance
import time
import DE_functions

import numpy as np
import numba

P_n=crowding_distance_assignment(F,R['CF_0'],R['CF_1'],nPop)

# importing the required module 
import timeit 
  
# code snippet to be executed only once 
mysetup = "import numpy as np,numba"

# code snippet whose execution time is to be measured 
mycode = '''
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
'''
# timeit statement 
timeit.timeit(setup = mysetup,stmt = mycode,number = 1000000)




