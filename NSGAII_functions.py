# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 15:56:04 2018

@author: Marios D. Pazidis
"""
from Pareto_functions import two_obj_dominance as two_obj_dominance

def fast_non_dominated_sorting(Cost_1):
    nPop_n=len(pop['Cost_Obj_01'])
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
            case=two_obj_dominance(pop['Cost_Obj_01'][p],pop['Cost_Obj_02'][p],pop['Cost_Obj_01'][q],pop['Cost_Obj_02'][q],cond)
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
    return F
