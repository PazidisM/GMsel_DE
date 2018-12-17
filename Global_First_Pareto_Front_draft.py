# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:27:25 2018

@author: Marios
"""

import numpy as np
from TOD import two_obj_dominance as two_obj_dominance

split_size=split_data['split_size']
nPop=DE_par['nPop']
cond='min'

index=0
index_spl=0
while index<NSeed:
    

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
        fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
        np.savetxt(fileName, w,formats['fmt_sf'])
#        p=0
#        while p<nPop-1:
#            if Dom[cmb,p]==1:
#                p=p+1
#                continue
#            q=p+1
#            while q<nPop:
#                if Dom[cmb,q]==1:
#                    q=q+1
#                    continue
#                case=two_obj_dominance(P['CF_0'][cmb,p],P['CF_1'][cmb,p],P['CF_0'][cmb,q],P['CF_1'][cmb,q],cond)
#                if case==0:
#                    Dom[cmb,q]=1
#                elif case==1 and Dom[cmb,p]==0:
#                    Dom[cmb,p]=1
#                q=q+1
#            p=p+1



Dom=np.zeros((split_end-split_start,nPop),dtype=int)
p_idx=0
while p_idx<Dom.size-1:
    p=p_idx%nPop
    cmb_p=p_idx//nPop
    if Dom[cmb_p,p]==1:
        p_idx=p_idx+1
        continue
    q_idx=p_idx+1
    q=q_idx%nPop
    cmb_q=q_idx//nPop
    while q_idx<Dom.size:
        if Dom[cmb_q,q]==1:
            q_idx=q_idx+1
            continue
        case=two_obj_dominance(P['CF_0'][cmb_p,p],P['CF_1'][cmb_p,p],P['CF_0'][cmb_q,q],P['CF_1'][cmb_q,q],cond)
        if case==0:
            Dom[cmb_q,q]=1
        elif case==1 and Dom[cmb_p,p]==0:
            Dom[cmb_p,p]=1
        q_idx=q_idx+1
        q=q_idx%nPop
        cmb_q=q_idx//nPop
    p_idx=p_idx+1
    if p_idx%1000==0:
        print(p_idx)


Front_CF_0=P['CF_0'][np.where(Dom==0)]
Front_CF_1=P['CF_1'][np.where(Dom==0)]
Front=np.column_stack((Front_CF_0,Front_CF_1))


arr=np.row_stack((foo_r,foo_c))



foo2=tuple(map(tuple, arr))

foo=Fr[0]
foo_r=Fr[0]//nPop
foo_c=Fr[0]%nPop
foo2=tuple([foo_r.astype(int),foo_c.astype(int)])

Front_CF_0=P['CF_0'][foo2]
Front_CF_1=P['CF_1'][foo2]
Front=np.column_stack((Front_CF_0,Front_CF_1))




P_CF_0=np.zeros((62,nPop),dtype=int)
P_CF_0=P['CF_0'].flatten()
P_CF_1=np.zeros((62,nPop),dtype=int)
P_CF_1=P['CF_1'].flatten()
Fr=fast_non_dominated_sorting(P_CF_0,P_CF_1,Dom.size)



from matplotlib import pyplot as plt

# unit area ellipse
rx, ry = 3., 1.
area = rx * ry * np.pi
theta = np.arange(0, 2 * np.pi + 0.01, 0.1)
verts = np.column_stack([rx / area * np.cos(theta), ry / area * np.sin(theta)])

x, y, s, c = np.random.rand(4, 70)
s=np.ones(70)
s *= 2**2.

x=Front[:,0]
y=Front[:,1]


fig, ax = plt.subplots()
ax.scatter(x, y,s)
#ax.scatter(x2, y2,s)
#ax.scatter(x3, y3,s)

plt.show()



fig = plt.figure()
ax = plt.axes()

x = T[0]
y=Sa_Tgt[0]

plt.plot(x, y);






