# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 15:15:21 2018

@author: Marios
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 14:28:46 2018

@author: Marios
"""
import numpy as np
import random
import numba
from CDA import crowding_distance_assignment as crowding_distance_assignment
#from TOD import two_obj_dominance as two_obj_dominance
from NSGAII_functions import fast_non_dominated_sorting as fast_non_dominated_sorting

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


@numba.jit(nopython=True, parallel=True)
def jDE_main_numba(cmb,Sa_Tgt,formats_fmt_sf,formats_fill_fn_all,folders_Scaling_factors,Max_sf_ind,Min_sf_ind,nGM,nPop,F_l,F_u,tau_1,tau_2,cond,MaxGen,index,P_CF_0,P_CF_1,P_Par_F,P_Par_CR,comb,w,Sa_unsc_ave_suite,Sa_suite):

    for gen in range(MaxGen):
    
        C_w=np.zeros((nGM, nPop))
        C_CF_0=np.full((1, nPop),np.inf)
        C_CF_1=np.full((1, nPop),np.inf)
        C_Par_F=np.full((1, nPop),np.inf)
        C_Par_CR=np.full((1, nPop),np.inf)
        
        for x in range(nPop):
            rand_1=random.uniform(0,1)
            rand_2=random.uniform(0,1)
            rand_3=random.uniform(0,1)
            rand_4=random.uniform(0,1)
            
            ## Mutation (DE/rand/1 for now) ##
            #Xr=random.sample([i for i in range(0,nPop) if i != x],3)
#            foo=list(range(0,nPop))
#            foo.remove(x)
            foo=np.arange(nPop)
            #foo.shape = (100,1)
            #foo=np.delete(foo,x,0)
            
            foo2=np.empty(nPop-1)
            f_id=0
            for kk in range(len(foo)):
                if foo[kk]!=x:
                    foo2[f_id]=foo[kk]
                    f_id=f_id+1
            #Xr=random.sample(foo,3)
            #Xr=np.random.choice(foo,3)
            Xr=foo2[np.random.choice(len(foo2),3)]
            
            if rand_2<tau_1:
                u_F=F_l+rand_1*F_u
            else:
                u_F=P_Par_F[cmb,x]
            
            u=w[:,int(Xr[0])]+u_F*(w[:,int(Xr[1])]-w[:,int(Xr[2])])   #DE/rand/1
            
            
            # enforce constraints in individual values
            case=u>Max_sf_ind[comb]
            u[case]=Max_sf_ind[comb[np.where(case)]]
            case=u<Min_sf_ind[comb]
            u[case]=Min_sf_ind[comb[np.where(case)]]
            
            ## Crossover ##
            if rand_4<tau_2:
                u_CR=rand_3
            else:
                u_CR=P_Par_CR[cmb,x]
            
            #case_01=random.randint(0, len(u)-1)
            case_01=np.random.randint(0, len(u))
            for k in range(len(u)):
                case_02=random.uniform(0,1)<=u_CR
                #case_02=np.random.uniform(0,1)<=u_CR
                if not (case_01==k or case_02):
                    u[k]=w[k,x]
#            u.shape=(nGM,1)
            u_rs=u.reshape(nGM, 1) 
            # Domination check
            u_c_01=(np.sum((Sa_unsc_ave_suite+np.sum(u_rs)/len(u_rs)-Sa_Tgt)**2)/len(Sa_Tgt[0]))**0.5
            u_c_02=(np.sum((Sa_suite+u_rs-(Sa_unsc_ave_suite+np.sum(u_rs)/len(u_rs)))**2/(len(u_rs)-1))/len(u_rs))**0.5


#            case=two_obj_dominance(P_CF_0[cmb,x],P_CF_1[cmb,x],u_c_01,u_c_02,cond)
            p_obj_01=P_CF_0[cmb,x]
            p_obj_02=P_CF_1[cmb,x]
            q_obj_01=u_c_01
            q_obj_02=u_c_02
            case_01=((p_obj_01<q_obj_01)&(p_obj_02<=q_obj_02))|((p_obj_01<=q_obj_01)&(p_obj_02<q_obj_02))
            case_02=((p_obj_01>q_obj_01)&(p_obj_02>=q_obj_02))|((p_obj_01>=q_obj_01)&(p_obj_02>q_obj_02))
            if cond=='min':
                if case_01:
                    case2=0    # p dominates q
                elif case_02:
                    case2=1    # q dominates p
                else:
                    case2=2    # nondominated
            elif cond=='max':
                if case_01:
                    case2=1    # q dominates p
                elif case_02:
                    case2=0    # p dominates q
                else:
                    case2=2    # nondominated

            if case2==1:
                #pass
                w[:,x]=u_rs[:,0]     # DEMO
                P_CF_0[cmb,x]=u_c_01
                P_CF_1[cmb,x]=u_c_02
                P_Par_F[cmb,x]=u_F
                P_Par_CR[cmb,x]=u_CR
            elif case2==2:
#                pass
                C_w[:,x]=u_rs[:,0]
                C_CF_0[0,x]=u_c_01
                C_CF_1[0,x]=u_c_02
                C_Par_F[0,x]=u_F
                C_Par_CR[0,x]=u_CR
            
        # combine parent and offspring populations
        case3=np.zeros(len(C_CF_0))
        for ll in range(len(C_CF_0)):
            if C_CF_0[0,ll]!=np.inf:
                case3[ll]=1
        C_CF_0[int(case3[:])]
#        case3=C_CF_0!=np.inf
#        R_CF_0=np.concatenate((P_CF_0[cmb,:],C_CF_0[case3]),axis=0)
#        R_CF_1=np.concatenate((P_CF_1[cmb,:],C_CF_1[case3]),axis=0)
#        R_Par_F=np.concatenate((P_Par_F[cmb,:],C_Par_F[case3]),axis=0)
#        R_Par_CR=np.concatenate((P_Par_CR[cmb,:],C_Par_CR[case3]),axis=0)
#        R_w=np.concatenate((w,C_w[:,np.where(case3)[1]]),axis=1)
#        np.where(case3)
#        C_w[:,np.where(case3)[1]]
#        # NSGAII - sort and truncate
#        F=fast_non_dominated_sorting(R_CF_0,R_CF_1,nPop,cond)
#        F2={}
#        for i in range(len(F)):
#            F2[i+1]=list(F[i].astype(int))
#        F=F2
#        P_n=crowding_distance_assignment(F,R_CF_0,R_CF_1,nPop)
#        
#        # form next generation population
#        
#        P_CF_0[cmb,:]=R_CF_0[P_n]
#        P_CF_1[cmb,:]=R_CF_1[P_n]
#        P_Par_F[cmb,:]=R_Par_F[P_n]
#        P_Par_CR[cmb,:]=R_Par_CR[P_n]
#        P_Par_CR[cmb,:]=R_w[P_n]
#
#        
#        w=R_w[:,P_n]
#    
#    fileName=folders_Scaling_factors+'\SF_'+str(index+cmb).zfill(formats_fill_fn_all)+'.out'
#    np.savetxt(fileName, w,formats_fmt_sf)
    return rand_4


jDE_main_numba(cmb,Sa_Tgt,formats_fmt_sf,formats_fill_fn_all,folders_Scaling_factors,Max_sf_ind,Min_sf_ind,nGM,nPop,F_l,F_u,tau_1,tau_2,cond,MaxGen,index,P_CF_0,P_CF_1,P_Par_F,P_Par_CR,comb,w,Sa_unsc_ave_suite,Sa_suite)

