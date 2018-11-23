# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 15:42:32 2018

@author: Marios
"""
import random
import math
import numpy as np
from tqdm import tqdm as tqdm
import time
import Cost_functions
from NSGAII_functions import fast_non_dominated_sorting as fast_non_dominated_sorting
from NSGAII_functions import crowding_distance_assignment as crowding_distance_assignment


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

#def jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt,Sa):
    
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
MaxGen=DE_par['MaxGen']

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
    fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
    P['CF_0']=np.loadtxt(fileName)
    fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
    P['CF_1']=np.loadtxt(fileName)
    
    for cmb in tqdm(range(len(Combs_split)),miniters =int(len(Combs_split)*0.05),desc='% of batch'):
        fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
        P['sf_cmb']=np.loadtxt(fileName)
        Sa_unsc_ave_suite=Sa_unsc_ave_split[cmb,:]
        Sa_suite=Sa[Combs_split[cmb],:]
        
        for gen in range(MaxGen):
        
            C={}
            C['sf_cmb']=np.zeros((nGM, nPop))
            C['CF_0']=np.full((1, nPop),np.inf)
            C['CF_1']=np.full((1, nPop),np.inf)
            C['Par_F']=np.full((1, nPop),np.inf)
            C['Par_CR']=np.full((1, nPop),np.inf)
            
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
                u[case]=P['sf_cmb'][case,x]
                u.shape=(nGM,1)
                u_mean=np.mean(u)
                
                # Domination check
                u_c_01=Cost_functions.CF_0(Sa_unsc_ave_suite,Sa_Tgt,u)
                u_c_02=Cost_functions.CF_1(Sa_suite,u,Sa_unsc_ave_suite)
                
                case=two_obj_dominance(P['CF_0'][cmb,x],P['CF_1'][cmb,x],u_c_01,u_c_02,cond)
                
                if case==1:
                    P['sf_cmb'][:,x]=u[:,0]     # DEMO
                    P['CF_0'][cmb,x]=u_c_01
                    P['CF_1'][cmb,x]=u_c_02
                    P['Par_F'][cmb,x]=u_F
                    P['Par_CR'][cmb,x]=u_CR
                elif case==2:
                    C['sf_cmb'][:,x]=u[:,0]
                    C['CF_0'][0,x]=u_c_01
                    C['CF_1'][0,x]=u_c_02
                    C['Par_F'][0,x]=u_F
                    C['Par_CR'][0,x]=u_CR
                
            # combine parent and offspring populations
            R={}
            case=C['CF_0']!=np.inf
            for key in [ v for v in C if v != 'sf_cmb' ]:
                R[key]=np.concatenate((P[key][cmb,:],C[key][case]),axis=0)
            R['sf_cmb']=np.concatenate((P['sf_cmb'],C['sf_cmb'][:,np.where(case)[1]]),axis=1)
            
            F=fast_non_dominated_sorting(R['CF_0'],R['CF_1'],cond)
            P_n=crowding_distance_assignment(F,R['CF_0'],R['CF_1'],nPop)
            
            for key in [ v for v in C if v != 'sf_cmb' ]:
                P[key][cmb,:]=R[key][P_n]
            P['sf_cmb']=R['sf_cmb'][:,P_n]
                            
                            
            
                            
                        
                        
                        
                        
                
                
                
                
                
                
                
                
                
                
                    
                    
                        ans=I['sorted_final'][nPop-nPop_trial:]
                        F[f_idx][I['sorted_final'][nPop-nPop_trial:]]
                        ans2=F[f_idx][k for k in ans]
                        
                        searchList=F[f_idx]
                        ans3=lambda searchList, ans: [searchList[i] for i in ans]
                        ans3(F[f_idx],ans)
                        
                        ans2= [ F[f_idx][i] for i in I['sorted_final'][nPop-nPop_trial:]]
                        
                        ans3=slice(I['sorted_final'][nPop-nPop_trial:])
                            I['sort_cd_1_idx'].index(i)
                            np.where(I['sort_cd_1_idx']==i)
                            
                            
                            
                            sorted(range(len(R['CF_0'][F[f_idx]])), key=lambda k: R['CF_0'][F[f_idx]][k])
                            I[i,1]=R['CF_0'][F[f_idx]]
                            R['CF_0'][F[f_idx]]
                            
                            
                            
                            
                        ans='Cost_Obj_'+str(1).zfill(2)
                        
                        I=np.stack((R['CF_0'][F[f_idx]], R['CF_1'][F[f_idx]]),axis=1)
                        I
                        
                        
                        
                        
                        
                        
                    
                    
                    
                    
                        
                    
                    
                
                
             % Crowding distance sorting
            nNewPop=0;
            NewPop=zeros(1,nPop);
            for kk=1:size(F,1)
                if nNewPop+nnz(F(kk,:))<=nPop
                    NewPop(nNewPop+1:nNewPop+nnz(F(kk,:)))=F(kk,1:nnz(F(kk,:)));
                    nNewPop=nNewPop+nnz(F(kk,:));
                else
                    n_acc=nPop-nNewPop;
                    front_mean=R_Cost_mean(F(kk,1:nnz(F(kk,:))));
                    [front_mean,I_mean]=sort(front_mean,'ascend');
                    front_indices_mean=F(kk,I_mean);
                    
                    front_std=R_Cost_std(F(kk,1:nnz(F(kk,:))));
                    [front_std,I_std]=sort(front_std,'ascend');
                    front_indices_std=F(kk,I_std);
                    
                    Icd_mean=zeros(length(front_mean),2);
                    Icd_mean(1,1)=Inf;
                    Icd_mean(end,1)=Inf;
                    Icd_mean(:,2)=front_indices_mean;
     
                    Icd_std=zeros(length(front_std),2);
                    Icd_std(1,1)=Inf;
                    Icd_std(end,1)=Inf;                    
                    Icd_std(:,2)=front_indices_std;
                    
                    Icd=zeros(length(front_std),2);

                    if length(front_mean)>2
                        for mm=2:length(front_mean)-1
                            Icd_mean(mm,1)=(front_mean(mm+1)-front_mean(mm-1))/(max(front_mean)-min(front_mean));
                            Icd_std(mm,1)=(front_std(mm+1)-front_std(mm-1))/(max(front_std)-min(front_std));
                        end
                    else
                        NewPop(nNewPop+1:end)=front_indices_mean(1);    %if only 2 elements choose 1 of them
                    end
                    Icd_mean=sortrows(Icd_mean,2);
                    Icd_std=sortrows(Icd_std,2);
                    
                    Icd(:,1)=Icd_mean(:,1)+Icd_std(:,1);
                    Icd(:,2)=Icd_mean(:,2);
                    Icd=sortrows(Icd,1);
                    NewPop(nNewPop+1:end)=Icd(end+1-n_acc:end,2);
                end
            end
            
            pop_S_sf=R(:,NewPop);
            F_mut(jj,:)=R_F_mut(:,NewPop);
            CR_par(jj,:)=R_CR_par(:,NewPop);
            F_mut_iter(ii,jj,:,it)=F_mut(jj,:);
            CR_par_iter(ii,jj,:,it)=CR_par(jj,:);

                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                            

                            

                            
                            
                            
                        
                        
                        
                        
                        
                        
                    
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
for kkk in range(100):
    for klk in range(10000):
        ans3=lambda searchList, ans: [searchList[i] for i in ans]
        ans3(F[f_idx],ans)

end=time.time()
dur_1=end-start
print(dur_1)


start=time.time()
for kkk in range(10):
    for klk in range(10000):
        ans2= [ F[f_idx][i] for i in ans]
end=time.time()
dur_2=end-start
print(dur_2)


print(dur_1/dur_2)

                        ans3(F[f_idx],ans)
                        
                        ans2= [ F[f_idx][i] for i in ans]


from sklearn.metrics import mean_squared_error
from math import sqrt

Sa_unsc_ave_suite2=Sa_unsc_ave_suite+np.mean(sf_ind_suite)
rms = sqrt(mean_squared_error(Sa_unsc_ave_suite2, Sa_Tgt))


Cost_functions.CF_0(Sa_unsc_ave_suite,Sa_Tgt,sf_ind_suite)







Sa_unsc_ave_suite.shape=(1,56)

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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
