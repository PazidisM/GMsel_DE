# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 15:42:32 2018

@author: Marios
"""
import random
import math
import numpy as np
import tqdm
import time
from Cost_functions import Obj_RMSE as Obj_RMSE
from Cost_functions import Obj_std as Obj_std
from Pareto_functions import two_obj_dominance as two_obj_dominance


def jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa):
    
    split_size=split_data['split_size']
    minSF=math.log(selectionParams['minScale'])
    maxSF=math.log(selectionParams['maxScale'])
    Max_sf_ind=selectionParams['Max_sf_ind']
    Min_sf_ind=selectionParams['Min_sf_ind']
    nGM=selectionParams['nGM']
    nPop=DE_par['nPop']
    F_l=DE_par['F_l']
    F_u=DE_par['F_u']
    CR_l=DE_par['CR_l']
    CR_u=DE_par['CR_u']
    tau_1=DE_par['tau_1']
    tau_2=DE_par['tau_2']
    cond='min'
    
    
    index=0
    index_spl=0
    while index<NSeed:
        
        print('\n')
        print('Optimizing Batch '+str(index_spl)+' of '+str(split_data['split_num']))
        
        split_start=index%split_size+index//split_size*split_size
        split_end=split_start+split_size
        if split_end>NSeed:
            split_end=NSeed
        
        P={}
        fileName=folders['Combinations']+'\Combs_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Combs_split=np.loadtxt(fileName).astype(int)
        fileName=folders['Sa_unsc_ave']+'\Sa_unsc_ave_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Sa_unsc_ave_split=np.loadtxt(fileName)
        fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        P['Par_F']=np.loadtxt(fileName)
        fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        P['Par_CR']=np.loadtxt(fileName)
        fileName=folders['Cost_Obj_01']+'\Cost_Obj_01_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        P['Cost_Obj_01']=np.loadtxt(fileName)
        fileName=folders['Cost_Obj_02']+'\Cost_Obj_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        P['Cost_Obj_02']=np.loadtxt(fileName)
        
        for cmb in tqdm(range(len(Combs_split)),miniters =int(len(Combs_split)*0.05),desc='% of batch'):
            fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
            P['sf_cmb']=np.loadtxt(fileName)
            Sa_unsc_ave_suite=Sa_unsc_ave_split[cmb,:]
            Sa_suite=Sa[Combs_split[cmb],:]
            
            for gen in reang(MaxGen):
                
                Q={}
                Q['sf_cmb']=np.zeros((nGM, nPop))
                Q['Cost_Obj_01']=np.full((1, nPop),np.inf)
                Q['Cost_Obj_02']=np.full((1, nPop),np.inf)
                Q['Par_F']=np.full((1, nPop),np.inf)
                Q['Par_CR']=np.full((1, nPop),np.inf)
                
                for x in range(nPop):
                    rand_1=random.uniform(0,1)
                    rand_2=random.uniform(0,1)
                    rand_3=random.uniform(0,1)
                    rand_4=random.uniform(0,1)
                    
                    ## Mutation (DE/rand/1 for now) ##
                    Xr=random.sample([i for i in range(0,nPop) if i != x],3)
                    if rand_2<tau_1:
                        u_F=F_l+rand_1*F_u
                    else:
                        u_F=P['Par_F'][cmb,x]
                    
                    u=P['sf_cmb'][:,Xr[0]]+u_F*(P['sf_cmb'][:,Xr[1]]-P['sf_cmb'][:,Xr[2]])   #DE/rand/1
                    
                    # enforce constraints in individual values
                    case=u>Max_sf_ind[Combs_split[cmb]]
                    u[case]=Max_sf_ind[Combs_split[cmb,np.where(case)][0]]
                    case=u<Min_sf_ind[Combs_split[cmb]]
                    u[case]=Min_sf_ind[Combs_split[cmb,np.where(case)][0]]
                    
                    ## Crossover ##
                    if rand_4<tau_2:
                        u_CR=rand_3
                    else:
                        u_CR=P['Par_CR'][cmb,x]
                    
                    case_01=random.choice(u)==u
                    case_02=np.array([random.random() for i in range(nGM)])<=u_CR
                    case=~np.array(case_01|case_02)
                    u[case]=P['sf'][case,x]
                    u.shape=(nGM,1)
                    u_mean=np.mean(u)
                    
                    # Domination check
                    u_c_01=Obj_RMSE(Sa_unsc_ave_suite,Sa_Tgt,u_mean)
                    u_c_02=Obj_std(Sa_suite,u,Sa_unsc_ave_suite,u_mean)
                    
                    case=two_obj_dominance(P['Cost_Obj_01'],P['Cost_Obj_02'],u_c_01,u_c_02,cond)
                    
                    if case==1:
                        P['sf_cmb'][:,x]=u[:,0]     # DEMO
                        P['Cost_Obj_01'][cmb,x]=u_c_01
                        P['Cost_Obj_02'][cmb,x]=u_c_02
                        P['Par_F'][cmb,x]=u_F
                        P['Par_CR'][cmb,x]=u_CR
                    elif case==2:
                        C['sf_cmb'][:,x]=u[:,0]
                        C['Cost_Obj_01'][0,x]=u_c_01
                        C['Cost_Obj_02'][0,x]=u_c_02
                        C['Par_F'][0,x]=u_F
                        C['Par_CR'][0,x]=u_CR
                    
                # combine parent and offspring populations
                R={}
                case=C['Cost_Obj_01']!=np.inf
                for key in [ v for v in Q if v != 'sf_cmb' ]:
                    R[key]=np.concatenate((P[key][cmb,:],C[key][case]),axis=0)
                R['sf_cmb']=np.concatenate((P['sf_cmb'],C['sf_cmb'][:,np.where(case)[1]]),axis=1)
                
                
                # Non-dominated sorting (NSGA-II)
                def fast_non_dominated_sorting(pop):
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
                        
                                
                                
                                
                                
                                
                                S[F[f_i]]
                            S[9,1]
                            F[f_i]
                            a=[i for i in F[f_i]]
                        S[a]
                        F[[i for i in F[f_i]]]
                        
                        a=[S[i] for i in F[f_i]]
                        a2=[[k for k in S[i][:]] for i in F[f_i]]
                        f_i=f_i+1
                    a3=S[9][:]
                    
                            
                            
                            
                            
                            
                            

                            

                            
                            
                            
                        
                        
                        
                        
                        
                        
                    
                ind=0
                while ind
                
                
                len(R['Cost_Obj_01'])
                
                
                
                
                % Non-dominated sorting (NSGA-II)
            Sp=zeros(size(R,2),size(R,2)-1);
            np=zeros(size(R,2),1);
            F=zeros(size(R,2),size(R,2));    % rows represent fronts, columns are members
            
            for indR_1=1:size(R,2)
                dom_ind=0;
                p_f1=R_Cost_mean(indR_1);
                p_f2=R_Cost_std(indR_1);
                for indR_2=1:size(R,2)
                    q_f1=R_Cost_mean(indR_2);
                    q_f2=R_Cost_std(indR_2);
                    if (p_f1<q_f1 && p_f2<=q_f2) || (p_f2<q_f2 && p_f1<=q_f1)
                        dom_ind=dom_ind+1;
                        Sp(indR_1,dom_ind)=indR_2;
                    elseif (q_f1<p_f1 && q_f2<=p_f2) || (q_f2<p_f2 && q_f1<=p_f1)
                        np(indR_1)=np(indR_1)+1;    % increment domination counter of p
                    end
                end
            end
            
            Sp=sort(Sp,2,'descend');    %non-zero elements first
            F(1,1:length(find(np==0)))=transpose(find(np==0));  % determine first non-dominated front
            f_ind=1;    % initialize front counter
            
            while any(F(f_ind,:))
                nfr_ind=0;
                Q_f=zeros(1,size(R,2));
                for p_ele_ind=1:nnz(F(f_ind,:))
                    for q_ele_ind=1:nnz(Sp(F(f_ind,p_ele_ind),:))
                        np(Sp(F(f_ind,p_ele_ind),q_ele_ind))= np(Sp(F(f_ind,p_ele_ind),q_ele_ind))-1;
                        if  np(Sp(F(f_ind,p_ele_ind),q_ele_ind))==0
                            nfr_ind=nfr_ind+1;
                            Q_f(1,nfr_ind)=Sp(F(f_ind,p_ele_ind),q_ele_ind);
                        end
                    end
                end
                f_ind=f_ind+1;
                F(f_ind,1:length(Q_f))=Q_f;
            end
            
            F(~any(F,2),:) = [];  %rows
            F(:,~any(F,1)) = [];  %columns
                
                
                
                
                
                
                
start=time.time()
summm=0
for kkk in range(50):
    for klk in range(10000):
        Xr=random.sample([i for i in range(0,nPop) if i != pm],3)
        if rand_2<tau_1:
            Par_F[cmb,pm]=F_l+rand_1*F_u
            
        u=sf[:,Xr[0]]+Par_F[cmb,pm]*(sf[:,Xr[1]]-sf[:,Xr[2]])   #DE/rand/1
        # enforce constraints in individual values
        case=u>Max_sf_ind[Combs_split[cmb]]
        u[case]=Max_sf_ind[Combs_split[cmb,[i for i in np.where(case)]]][0]
end=time.time()
dur_1=end-start
print(dur_1)


start=time.time()
summm=0
for kkk in range(50):
    for klk in range(10000):
        Xr=random.sample([i for i in range(0,nPop) if i != pm],3)
        if rand_2<tau_1:
            Par_F[cmb,pm]=F_l+rand_1*F_u
            
        u=sf[:,Xr[0]]+Par_F[cmb,pm]*(sf[:,Xr[1]]-sf[:,Xr[2]])   #DE/rand/1
        # enforce constraints in individual values
        case=u>Max_sf_ind[Combs_split[cmb]]
        u[case]=Max_sf_ind[Combs_split[cmb,np.where(case)][0]]
end=time.time()
dur_2=end-start
print(dur_2)
print(dur_1/dur_2)








        for i=1:nPop

                rand_1=unifrnd(0,1);
                rand_2=unifrnd(0,1);
                rand_3=unifrnd(0,1);
                rand_4=unifrnd(0,1);
                
                p=pop_S_sf(:,i);

                A=randperm(nPop);
                A(A==i)=[];
                a=A(1);
                b=A(2);
                c=A(3);

                % Mutation
                if rand_2<tau_1
                    F_mut(jj,i)=F_l+rand_1*F_u;
                else
                    F_mut(jj,i)=F_mut(jj,i);
                end
                
                y=pop_S_sf(:,a)+F_mut(jj,i)*(pop_S_sf(:,b)-pop_S_sf(:,c));
                %y = max(y, VarMin);
                %y = min(y, VarMax);
                
                y = max(y, VarMin);
                y = min(y, Max_ind_Sf(GM_Trials_GMs_SLICED,1));

                % Crossover
                if rand_4<tau_2
                    CR_par(jj,i)=rand_3;
                else
                    CR_par(jj,i)=CR_par(jj,i);
                end
                z=zeros(size(p));
                j0=randi([1 numel(p)]);
                
                for j=1:numel(p)
                    if j==j0 || rand<=CR_par(jj,i)
                        z(j)=y(j);
                    else
                        z(j)=p(j);
                    end
                end            

                q=z;
                q_Cost_mean=Cost_mean_log(Sa_Tgt_log,GM_Trials_Sa_ave_unsc_log_SLICED,mean(q),T_match);
                q_Cost_std=Cost_std_log(GM_Trials_Sa_ave_unsc_log_SLICED,Sa_db_log_SLICED,q,T_match);

    
    
    
    
    
     
    
    
    
        % select file folder for GM arrays
    foldername=sprintf('%1$s\\%2$s\\Split_arrays\\spl_%3$04.f',currentFolder,spl_dir,ii);
    
    % load files
    filename=sprintf('%1$s\\spl_%2$04.f_GMs.out',foldername,ii);
    GM_Trials_GMs=load(filename);
    filename=sprintf('%1$s\\spl_%2$04.f_Sa_ave_unsc_log.out',foldername,ii);
    GM_Trials_Sa_ave_unsc_log=load(filename);
    
    % Compute size
    NSeed_split=size(GM_Trials_GMs,1);
    
    % select file folder for population arrays
    foldername=sprintf('%1$s\\%2$s\\Split_pop\\splpop_%3$04.f',currentFolder,spl_dir,ii);
    
    % load pop cost files
    filename=sprintf('%1$s\\splpop_%2$04.f_Cost_mean.out',foldername,ii);
    pop_Cost_mean=load(filename);
    
    % load pop cost files
    filename=sprintf('%1$s\\splpop_%2$04.f_F_mut.out',foldername,ii);
    F_mut=load(filename);
    
    % load pop cost files
    filename=sprintf('%1$s\\splpop_%2$04.f_CR_par.out',foldername,ii);
    CR_par=load(filename);
            
    % load pop cost files
    filename=sprintf('%1$s\\splpop_%2$04.f_Cost_std.out',foldername,ii);
    pop_Cost_std=load(filename);
    
    % select file folder for scaling factors arrays
    foldername=sprintf('%1$s\\%2$s\\Split_pop\\splpop_%3$04.f\\popSF',currentFolder,spl_dir,ii);

    for jj=1:NSeed_split
        
        filename=sprintf('%1$s\\splpop_%2$04.f_suite_%3$06.f.out',foldername,ii,jj);
        pop_S_sf=load(filename);
        
        GM_Trials_Sa_ave_unsc_log_SLICED=GM_Trials_Sa_ave_unsc_log(jj,:);
        GM_Trials_GMs_SLICED=GM_Trials_GMs(jj,:);
        Sa_db_log_SLICED=Sa_db_log(GM_Trials_GMs_SLICED,:);
        

        if  mod(jj,floor(split.N/20))==0
            waitbar(jj/NSeed_split,hw,sprintf('Optimizing sets from split %1$g of %2$g  ...%3$g %%',ii,split.div,jj/NSeed_split*100)); % update waitbar
        end
    
        for it=1:MaxIt
            
        

            Q=zeros(selectionParams.nGM,nPop);
            Q_Cost_mean=inf(1,nPop);
            Q_Cost_std=inf(1,nPop);
            Q_F_mut=inf(1,nPop);
            Q_CR_par=inf(1,nPop);
            
            for i=1:nPop
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
