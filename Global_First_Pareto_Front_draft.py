# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 14:27:25 2018

@author: Marios
"""

import numpy as np
from TOD import two_obj_dominance as two_obj_dominance
from NSGAII_functions import fast_non_dominated_sorting as fast_non_dominated_sorting
from NSGAII_functions import fast_non_dominated_sorting_first_front as fast_non_dominated_sorting_first_front

split_size=split_data['split_size']
nPop=DE_par['nPop']
cond='min'

index=0
index_spl=0
while index<NSeed:
    
    Fr=np.empty(0)
    Fr=np.reshape(Fr,[0,1])
    split_start=index%split_size+index//split_size*split_size
    split_end=split_start+split_size
    if split_end>NSeed:
        split_end=NSeed
    
    P={}
    fileName=folders['CF_0']+'\CF_0_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
    P['CF_0']=np.loadtxt(fileName)
    fileName=folders['CF_1']+'\CF_1_'+str(index_spl).zfill(formats['fill_fn_split'])+'.out'
    P['CF_1']=np.loadtxt(fileName)
    P_CF_0=P['CF_0'].flatten()
    P_CF_1=P['CF_1'].flatten()
    foo=fast_non_dominated_sorting_first_front(P_CF_0,P_CF_1,nPop)
    foo[:,0]=foo[:,0]+split_start
    Fr=np.vstack((Fr,np.reshape(foo[:,0],[len(foo),1])))
    index_spl=index_spl+1
    index=split_end


leng=len(P_CF_0)


foo=Fr_buffer[:,0]

foo2=np.reshape(foo,[40,1])
fooo=np.reshape(Fr,[0,1])
len(foo)
Fr=np.vstack((fooo,foo2))

foo=F

Cost_1=P_CF_0
Cost_2=P_CF_1
nPop=len(P_CF_0)

Fr_buffer=F

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






