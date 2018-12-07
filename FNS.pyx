# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 15:31:06 2018

@author: Marios
"""
import numpy as np

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



