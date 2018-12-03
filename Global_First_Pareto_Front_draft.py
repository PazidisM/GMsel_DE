# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:27:25 2018

@author: Marios
"""
import numpy as np

split_size=split_data['split_size']
nPop=DE_par['nPop']
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
    fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
    P['CF_0']=np.loadtxt(fileName)
    fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
    P['CF_1']=np.loadtxt(fileName)
    
    Dom=np.zeros((split_end-split_start,nPop),dtype=int)
    for cmb in range(split_end-split_start):
        
        p=0
        while p<nPop-1:
            if Dom[cmb,p]==1:
                p=p+1
                continue
            q=p+1
            while q<nPop:
                if Dom[cmb,q]==1:
                    q=q+1
                    continue
                case=two_obj_dominance(P['CF_0'][cmb,p],P['CF_1'][cmb,p],P['CF_0'][cmb,q],P['CF_1'][cmb,q],cond)
                if case==0:
                    Dom[cmb,q]=1
                elif case==1 and Dom[cmb,p]==0:
                    Dom[cmb,p]=1
                q=q+1
            p=p+1
    index_spl=index_spl+1
    index=split_end
    
    
    



    Dom=np.zeros((split_end-split_start,nPop),dtype=int)
    for cmb in range(split_end-split_start):
        
        p=0
        while p<nPop-1:
            if Dom[cmb,p]==1:
                p=p+1
                continue
            q=p+1
            while q<nPop:
                if Dom[cmb,q]==1:
                    q=q+1
                    continue
                case=two_obj_dominance(P['CF_0'][cmb,p],P['CF_1'][cmb,p],P['CF_0'][cmb,q],P['CF_1'][cmb,q],cond)
                if case==0:
                    Dom[cmb,q]=1
                elif case==1 and Dom[cmb,p]==0:
                    Dom[cmb,p]=1
                q=q+1
            p=p+1




P['CF_1']
Dom==0


[i for i, e in enumerate(a) if e != 0]

ans=[i for i, e in enumerate(Dom[0])]

 e in enumerate(Dom[0])

import matplotlib.pyplot as plt

# unit area ellipse
rx, ry = 3., 1.
area = rx * ry * np.pi
theta = np.arange(0, 2 * np.pi + 0.01, 0.1)
verts = np.column_stack([rx / area * np.cos(theta), ry / area * np.sin(theta)])

x, y, s, c = np.random.rand(4, 70)
s=np.ones(70)
s *= 2**2.

x=Summary[:,0]
y=Summary[:,1]
=

fig, ax = plt.subplots()
ax.scatter(x, y,s)
ax.scatter(x2, y2,s)
ax.scatter(x3, y3,s)

plt.show()







