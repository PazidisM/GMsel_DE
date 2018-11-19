# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 15:42:32 2018

@author: Marios
"""
import random

def jDE(selectionParams,DE_par,NSeed,folders,formats,split_data,Sa_Tgt):
    
    split_size=split_data['split_size']
    
    minSF=math.log(selectionParams['minScale'])
    maxSF=math.log(selectionParams['maxScale'])
    Max_sf_ind=selectionParams['Max_sf_ind']
    Min_sf_ind=selectionParams['Min_sf_ind']
    nGM=selectionParams['nGM']
    nPop=DE_par['nPop']
    F_in=DE_par['F_in']
    CR_in=DE_par['CR_in']
    
    index=0
    index_spl=0
    while index<NSeed:
        split_start=index%split_size+index//split_size*split_size
        split_end=split_start+split_size
        if split_end>NSeed:
            split_end=NSeed
        
        print('\n')
        print('Optimizing Batch '+str(index_spl)+' of '+str(split_data['split_num']))
        
        fileName=folders['Combinations']+'\Combs_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Combs_split=np.loadtxt(fileName).astype(int)
        fileName=folders['Sa_unsc_ave']+'\Sa_unsc_ave_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Sa_unsc_ave_split=np.loadtxt(fileName)
        
        fileName=folders['Par_F']+'\Par_F_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Par_F=np.loadtxt(fileName)
        fileName=folders['Par_CR']+'\Par_CR_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Par_CR=np.loadtxt(fileName)
        fileName=folders['Cost_Obj_01']+'\Cost_Obj_01_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Cost_Obj_01=np.loadtxt(fileName)
        fileName=folders['Cost_Obj_02']+'\Cost_Obj_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
        Cost_Obj_02=np.loadtxt(fileName)
        
        for ii in tqdm(range(len(Combs_split)),miniters =int(len(Combs_split)*0.05),desc='% of batch'):
            fileName=folders['Scaling_factors']+'\SF_'+str(index+ii).zfill(formats['fill_fn_all'])+'.out'
            sf=np.loadtxt(fileName)
            Sample_sf=np.mean(sf,axis=0)
            
            for gen in reang(MaxGen):
                
                Q={}
                Q['sf']=np.zeros((nGM, nPop))
                Q['Cost_Obj_01']=np.full((1, nPop),np.inf)
                Q['Cost_Obj_02']=np.full((1, nPop),np.inf)
                Q['Par_F']=np.full((1, nPop),np.inf)
                Q['Par_CR']=np.full((1, nPop),np.inf)
                
                for pm in range(nPop):
                    rand_1=random.uniform(0,1)
                    rand_2=random.uniform(0,1)
                    rand_3=random.uniform(0,1)
                    rand_4=random.uniform(0,1)
                    
                    [random.sample([i for i in range(0,nPop) if i != pm],3)]
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    








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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
