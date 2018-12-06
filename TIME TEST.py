# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 10:08:23 2018

@author: Marios
"""


st=time()
for i in range(rep):
    u_c_01=Cost_functions.CF_0(Sa_unsc_ave_suite,Sa_Tgt,u)

end=time()
dur=end-st
print(dur)



def CF_0(Sa_ave,Sa_Tgt,sf_suite):
    
    sf_ave=np.mean(sf_suite)
    diff=Sa_ave+sf_ave-Sa_Tgt
    z=math.sqrt(np.sum(np.power(diff,2))/np.size(Sa_Tgt,1))
    
    return z


tot_dur=0
for w in range(20):
    st=time()
    for i in range(rep):
        
#        sf_ave=np.mean(u)
#        sf_ave=np.sum(u)/np.size(u,0)
#        sf_ave=np.sum(u)/len(u)
#        diff=Sa_unsc_ave_suite+sf_ave-Sa_Tgt
#        z=math.sqrt(np.sum(np.power(diff,2))/np.size(Sa_Tgt,1))
#        z=(np.sum(np.power(diff,2))/np.size(Sa_Tgt,1))**0.5
        
#        z=(np.sum(diff**2)/len(Sa_Tgt[0]))**0.5
#        z=(np.sum(diff**2)/np.size(Sa_Tgt,1))**0.5
#        z=(np.sum((Sa_unsc_ave_suite+sf_ave-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5
        z=(np.sum((Sa_unsc_ave_suite+np.sum(u)/len(u)-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5

        
#        u_c_01=Cost_functions.CF_0(Sa_unsc_ave_suite,Sa_Tgt,u)
    
    end=time()
    dur=end-st
    tot_dur=tot_dur+dur
    print(dur)
    
print('\n')
print(tot_dur/20)


np.power(diff,2)


np.size(Sa_Tgt,1)
len(Sa_Tgt[0])



tot_dur=0
for w in range(5):
    st=time()
    for i in range(rep):
        
#        u_c_02=Cost_functions.CF_1(Sa_suite,u,Sa_unsc_ave_suite)
        
        
        sf_suite_mean=np.sum(sf_suite)/len(sf_suite)
        sf_suite.shape=(len(sf_suite),1)
        
#        diff=Sa_suite+sf_suite-(Sa_unsc_ave_suite+sf_suite_mean)
        
        
        
        
#        var=diff**2/(len(sf_suite)-1)
        
#        
#       var_mean=np.sum(var,0)/len(var)
        
#       z=var_mean**0.5
        
#        z=(np.sum(var)/len(var))**0.5
        
#        z=(np.sum(diff**2/(len(sf_suite)-1))/len(var))**0.5
        
        z=(np.sum((Sa_suite+sf_suite-(Sa_unsc_ave_suite+sf_suite_mean))**2/(len(sf_suite)-1))/len(sf_suite))**0.5
#       

#         u_c_02=Cost_functions.CF_1(Sa_suite,u,Sa_unsc_ave_suite)
        
#        z=(np.sum((Sa_suite+sf_suite-(Sa_unsc_ave_suite+np.sum(sf_suite)/len(sf_suite)))**2/(len(sf_suite)-1))/len(sf_suite))**0.5
        
    end=time()
    dur=end-st
    tot_dur=tot_dur+dur
    print(dur)
    
print('\n')
print(tot_dur/5)


        z=(np.sum((Sa_suite+sf_suite-(Sa_ave+sf_suite_mean))**2/(len(sf_suite)-1))/len(sf_suite))**0.5


ans=np.sum(var,1)/len(var)
ans=np.sum(np.sum(var,1))/len(var)
ans=np.sum(var)/len(var)

np.size(sf_suite,0)



        sf_suite_mean=np.mean(sf_suite)
        sf_suite.shape=(len(sf_suite),1)
        diff=Sa_suite+sf_suite-(Sa_unsc_ave_suite+sf_suite_mean)
        var=np.power(diff,2)/(len(sf_suite)-1)
        var_mean=np.mean(var)
        z=np.sqrt(var_mean)
    

ans=(P['sf_cmb'][:,Xr[0]]+u_F*(P['sf_cmb'][:,Xr[1]]-P['sf_cmb'][:,Xr[2]]))

u_c_02=Cost_functions.CF_1(Sa_suite,u,Sa_unsc_ave_suite)

sf_suite=ans

        case_01=((P['CF_0'][cmb,x]<u_c_01)and(P['CF_1'][cmb,x]<=u_c_02))or((P['CF_0'][cmb,x]<=u_c_01)and(P['CF_1'][cmb,x]<u_c_02))
        case_02=((P['CF_0'][cmb,x]>u_c_01)and(P['CF_1'][cmb,x]>=u_c_02))or((P['CF_0'][cmb,x]>=u_c_01)and(P['CF_1'][cmb,x]>u_c_02))
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



import Cost_functions 
from time import time as time
import numpy as np
import math
import statistics

import scipy as sp

ans=[random.random() for i in range(nGM)]
ans=[random.random() for i in range(nGM)]

w=P['sf_cmb']











tot_dur=0
tot_dur2=0
#st=time.clock()

for gen in range(MaxGen):

    C={}
    C['w']=np.zeros((nGM, nPop))
    C['CF_0']=np.full((1, nPop),np.inf)
    C['CF_1']=np.full((1, nPop),np.inf)
    C['Par_F']=np.full((1, nPop),np.inf)
    C['Par_CR']=np.full((1, nPop),np.inf)
    st=time.clock()

    for x in range(nPop):
        rand_1=random.uniform(0,1)
        rand_2=random.uniform(0,1)
        rand_3=random.uniform(0,1)
        rand_4=random.uniform(0,1)
        
        ## Mutation (DE/rand/1 for now) ##
        #Xr=random.sample([i for i in range(0,nPop) if i != x],3)
        ans=list(range(0,nPop))
        ans.remove(x)
        Xr=random.sample(ans,3)
        if rand_2<tau_1:
            u_F=F_l+rand_1*F_u
        else:
            u_F=P['Par_F'][cmb,x]
        
        u=w[:,Xr[0]]+u_F*(w[:,Xr[1]]-w[:,Xr[2]])   #DE/rand/1
        
        # enforce constraints in individual values
        case=u>Max_sf_ind[comb]
        u[case]=Max_sf_ind[comb[np.where(case)]]
        case=u<Min_sf_ind[comb]
        u[case]=Min_sf_ind[comb[np.where(case)]]
        
        ## Crossover ##
        if rand_4<tau_2:
            u_CR=rand_3
        else:
            u_CR=P['Par_CR'][cmb,x]
        
        case_01=random.randint(0, len(u)-1)
        for k in range(len(u)):
            case_02=random.uniform(0,1)<=u_CR
            if not (case_01==k or case_02):
                u[k]=w[k,x]
        u.shape=(nGM,1)
        
        # Domination check
        u_c_01=(np.sum((Sa_unsc_ave_suite+np.sum(u)/len(u)-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5
        u_c_02=(np.sum((Sa_suite+u-(Sa_unsc_ave_suite+np.sum(u)/len(u)))**2/(len(u)-1))/len(u))**0.5
        
        case=two_obj_dominance(P['CF_0'][cmb,x],P['CF_1'][cmb,x],u_c_01,u_c_02,cond)
        
        if case==1:
            w[:,x]=u[:,0]     # DEMO
            P['CF_0'][cmb,x]=u_c_01
            P['CF_1'][cmb,x]=u_c_02
            P['Par_F'][cmb,x]=u_F
            P['Par_CR'][cmb,x]=u_CR
        elif case==2:
            C['w'][:,x]=u[:,0]
            C['CF_0'][0,x]=u_c_01
            C['CF_1'][0,x]=u_c_02
            C['Par_F'][0,x]=u_F
            C['Par_CR'][0,x]=u_CR
    end=time.clock()
    dur=end-st
    tot_dur=tot_dur+dur
    
    st2=time.clock()
    # combine parent and offspring populations
    R={}
    case=C['CF_0']!=np.inf
    for key in [ v for v in C if v != 'w' ]:
        R[key]=np.concatenate((P[key][cmb,:],C[key][case]),axis=0)
    R['w']=np.concatenate((w,C['w'][:,np.where(case)[1]]),axis=1)
    end2=time.clock()

    # NSGAII - sort and truncate
    F=fast_non_dominated_sorting(R['CF_0'],R['CF_1'],cond)
    P_n=crowding_distance_assignment(F,R['CF_0'],R['CF_1'],nPop)
    
    # form next generation population
    for key in [ v for v in C if v != 'w' ]:
        P[key][cmb,:]=R[key][P_n]
    w=R['w'][:,P_n]
    
    
    
    dur2=end2-st2
    tot_dur2=tot_dur2+dur2
    
    fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
np.savetxt(fileName, w,formats['fmt_sf'])

#end=time.clock()
#dur=end-st
#tot_dur=tot_dur+dur


print('\n')
print(tot_dur)

print('\n')
print(tot_dur2)



tot_dur=0
tot_dur2=0
#st=time.clock()

for gen in range(MaxGen):

    C={}
    C['w']=np.zeros((nGM, nPop))
    C['CF_0']=np.full((1, nPop),np.inf)
    C['CF_1']=np.full((1, nPop),np.inf)
    C['Par_F']=np.full((1, nPop),np.inf)
    C['Par_CR']=np.full((1, nPop),np.inf)
    st=time.clock()

    for x in range(nPop):
        rand_1=random.uniform(0,1)
        rand_2=random.uniform(0,1)
        rand_3=random.uniform(0,1)
        rand_4=random.uniform(0,1)
        
        ## Mutation (DE/rand/1 for now) ##
        #Xr=random.sample([i for i in range(0,nPop) if i != x],3)
        ans=list(range(0,nPop))
        ans.remove(x)
        Xr=random.sample(ans,3)
        if rand_2<tau_1:
            u_F=F_l+rand_1*F_u
        else:
            u_F=P['Par_F'][cmb,x]
        
        u=w[:,Xr[0]]+u_F*(w[:,Xr[1]]-w[:,Xr[2]])   #DE/rand/1
        
        # enforce constraints in individual values
        case=u>Max_sf_ind[comb]
        u[case]=Max_sf_ind[comb[np.where(case)]]
        case=u<Min_sf_ind[comb]
        u[case]=Min_sf_ind[comb[np.where(case)]]
        
        ## Crossover ##
        if rand_4<tau_2:
            u_CR=rand_3
        else:
            u_CR=P['Par_CR'][cmb,x]
        
        case_01=random.randint(0, len(u)-1)
        for k in range(len(u)):
            case_02=random.uniform(0,1)<=u_CR
            if not (case_01==k or case_02):
                u[k]=w[k,x]
        u.shape=(nGM,1)
        
        # Domination check
        u_c_01=(np.sum((Sa_unsc_ave_suite+np.sum(u)/len(u)-Sa_Tgt)**2)/np.size(Sa_Tgt,1))**0.5
        u_c_02=(np.sum((Sa_suite+u-(Sa_unsc_ave_suite+np.sum(u)/len(u)))**2/(len(u)-1))/len(u))**0.5
        
        case=two_obj_dominance(P['CF_0'][cmb,x],P['CF_1'][cmb,x],u_c_01,u_c_02,cond)
        
        if case==1:
            w[:,x]=u[:,0]     # DEMO
            P['CF_0'][cmb,x]=u_c_01
            P['CF_1'][cmb,x]=u_c_02
            P['Par_F'][cmb,x]=u_F
            P['Par_CR'][cmb,x]=u_CR
        elif case==2:
            C['w'][:,x]=u[:,0]
            C['CF_0'][0,x]=u_c_01
            C['CF_1'][0,x]=u_c_02
            C['Par_F'][0,x]=u_F
            C['Par_CR'][0,x]=u_CR
    end=time.clock()
    dur=end-st
    tot_dur=tot_dur+dur
    
    # combine parent and offspring populations
    R={}
    case=C['CF_0']!=np.inf
    for key in [ v for v in C if v != 'w' ]:
        R[key]=np.concatenate((P[key][cmb,:],C[key][case]),axis=0)
    R['w']=np.concatenate((w,C['w'][:,np.where(case)[1]]),axis=1)
    

    # NSGAII - sort and truncate
    st2=time.clock()

    F=fast_non_dominated_sorting(R['CF_0'],R['CF_1'],cond)
#
    end2=time.clock()
    dur2=end2-st2
    print(dur2)
    tot_dur2=tot_dur2+dur2
    
    
    

    
    P_n=crowding_distance_assignment(F,R['CF_0'],R['CF_1'],nPop)
    
    # form next generation population
    for key in [ v for v in C if v != 'w' ]:
        P[key][cmb,:]=R[key][P_n]
    w=R['w'][:,P_n]
    
    
    

    
    fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
np.savetxt(fileName, w,formats['fmt_sf'])

#end=time.clock()
#dur=end-st
#tot_dur=tot_dur+dur


print('\n')
print(tot_dur)

print('\n')
print(tot_dur2)

st2=time.clock()
end2=time.clock()
dur2=end2-st2
print(dur2)
tot_dur2=tot_dur2+dur2



F=fast_non_dominated_sorting(R['CF_0'],R['CF_1'],cond)


def fast_non_dominated_sorting():
    global R
    global cond
    Cost_1=R['CF_0']
    Cost_2=R['CF_1']
    cond2=cond
    
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
            case=two_obj_dominance(Cost_1[p],Cost_2[p],Cost_1[q],Cost_2[q],cond2)
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

def fast_non_dominated_sorting():
    global R
    global cond
    Cost_1=R['CF_0']
    Cost_2=R['CF_1']
    cond2=cond
    
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
            p_obj_01=Cost_1[p]
            p_obj_02=Cost_2[p]
            q_obj_01=Cost_1[q]
            q_obj_02=Cost_2[q]
            case_011=((p_obj_01<q_obj_01)and(p_obj_02<=q_obj_02))
            case_012=((p_obj_01<=q_obj_01)and(p_obj_02<q_obj_02))
            case_021=((p_obj_01>q_obj_01)and(p_obj_02>=q_obj_02))
            case_022=((p_obj_01>=q_obj_01)and(p_obj_02>q_obj_02))
            
            case_01=case_011 or case_012
            case_02=case_021 or case_022
                        
            if case_01:
                case=0    # p dominates q
            elif case_02:
                case=1    # q dominates p
            else:
                case=2    # nondominated

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

%lprun -f fast_non_dominated_sorting fast_non_dominated_sorting(C1,C2,cond)
%lprun -f jDE jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)
%lprun -f two_obj_dominance two_obj_dominance(p1,p2,q1,q2,cond)

fast_non_dominated_sorting(R['CF_0'],R['CF_1'],cond)
C1=R['CF_0']
C2=R['CF_1']

%load_ext line_profiler
load_ext line_profiler
p1=R['CF_0'][0]
p2=R['CF_1'][0]
q1=R['CF_0'][1]
q2=R['CF_1'][1]



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



def two_obj_dominance(p_obj_01,p_obj_02,q_obj_01,q_obj_02,cond):
    
    case_01=False
    case_02=False
    if (p_obj_01<q_obj_01):
        if (p_obj_02<=q_obj_02):
            case_01=True
        else:
            case_01=False
            case_02=False
    elif (p_obj_01<=q_obj_01):
        if (p_obj_02<q_obj_02):
            case_01=True
        else:
            case_01=False
            case_02=False
    elif (p_obj_01>q_obj_01):
        if (p_obj_02>=q_obj_02):
            case_02=True
        else:
            case_01=False
            case_02=False
    elif (p_obj_01>=q_obj_01):
        if (p_obj_02>q_obj_02):
            case_02=True
        else:
            case_01=False
            case_02=False
    
    if case_01:
        case=0    # p dominates q
    elif case_02:
        case=1    # q dominates p
    else:
        case=2    # nondominated
    
    return case


def two_obj_dominance(p_obj_01,p_obj_02,q_obj_01,q_obj_02,cond):
    
    case_01=False
    case_02=False
    if (q_obj_01>p_obj_01):
        if (p_obj_02<=q_obj_02):
            case_01=True
        else:
            case_01=False
            case_02=False
    elif (p_obj_01<=q_obj_01):
        if (p_obj_02<q_obj_02):
            case_01=True
        else:
            case_01=False
            case_02=False
    elif (p_obj_01>q_obj_01):
        if (p_obj_02>=q_obj_02):
            case_02=True
        else:
            case_01=False
            case_02=False
    elif (p_obj_01>=q_obj_01):
        if (p_obj_02>q_obj_02):
            case_02=True
        else:
            case_01=False
            case_02=False
    
    if case_01:
        case=0    # p dominates q
    elif case_02:
        case=1    # q dominates p
    else:
        case=2    # nondominated
    
    return case









fast_non_dominated_sorting(C1,C2,cond)

%lprun -f fast_non_dominated_sorting fast_non_dominated_sorting(C1,C2,cond)
%lprun -f fast_non_dominated_sorting_numbla fast_non_dominated_sorting_numbla(C1,C2,cond)

@jit
fast_non_dominated_sorting_numbla=jit()(fast_non_dominated_sorting)
fast_non_dominated_sorting_numbla(C1,C2,cond)


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






















p_obj_01=Cost_1[p]
p_obj_02=Cost_2[p]
q_obj_01=Cost_1[q]
q_obj_02=Cost_2[q]
case_01=((p_obj_01<q_obj_01)and(p_obj_02<=q_obj_02))or((p_obj_01<=q_obj_01)and(p_obj_02<q_obj_02))
case_02=((p_obj_01>q_obj_01)and(p_obj_02>=q_obj_02))or((p_obj_01>=q_obj_01)and(p_obj_02>q_obj_02))
if cond2=='min':
    if case_01:
        case=0    # p dominates q
    elif case_02:
        case=1    # q dominates p
    else:
        case=2    # nondominated
elif cond2=='max':
    if case_01:
        case=1    # q dominates p
    elif case_02:
        case=0    # p dominates q
    else:
        case=2    # nondominated





import line_profiler

rep=100
tot_dur2=0
tot_dur=0
for rrrr in range(5):
    st=time.clock()
    for i in range(rep):
        
#        F=fast_non_dominated_sorting(R['CF_0'],R['CF_1'],cond)
        
        st2=time.clock()
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
                
        end2=time.clock()
        dur2=end2-st2
        tot_dur2=tot_dur2+dur2
        
        
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

        
    end=time.clock()
    dur=end-st
    tot_dur=tot_dur+dur
    print(dur)
    
print('\n')
print(tot_dur/5)

print('\n')
print(tot_dur2)




rep=100
tot_dur=0
st=time.clock()
for i in range(rep):
    #fast_non_dominated_sorting(C1,C2,cond)
    fast_non_dominated_sorting_numbla(C1,C2,cond)
end=time.clock()
dur=end-st
print(dur)


%lprun -f two_obj_dominance two_obj_dominance(p1,p2,q1,q2,cond)

import time



two_obj_dominance_numbla=jit(two_obj_dominance)

DE_functions.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)

st=time.clock()
for i in range(rep):

    #two_obj_dominance(p1,p2,q1,q2,cond)
    #two_obj_dominance_numbla(p1,p2,q1,q2,cond)
    #fast_non_dominated_sorting(C1,C2,cond)
    fast_non_dominated_sorting_numbla(C1,C2,cond)
#DE_functions.jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)
#jDE_numbla(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa)

end=time.clock()
dur=end-st
print(dur)

fast_non_dominated_sorting(C1,C2,cond)
fast_non_dominated_sorting_numbla(C1,C2,cond)

%lprun -f fast_non_dominated_sorting fast_non_dominated_sorting(C1,C2,cond)
%lprun -f fast_non_dominated_sorting_numbla fast_non_dominated_sorting_numbla(C1,C2,cond)

@jit
fast_non_dominated_sorting_numbla=jit(fast_non_dominated_sorting)
fast_non_dominated_sorting_numbla(C1,C2,cond)






