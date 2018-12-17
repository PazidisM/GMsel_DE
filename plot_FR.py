# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 10:46:14 2018

@author: Marios
"""



FR_cmb=[int(i//nPop) for i in Fr[0]]
FR_p=[int(i%nPop) for i in Fr[0]]

SF=np.empty((nGM,len(Fr[0])))
Sa_un=np.empty((len(Fr[0]),56))


index=0
while index<NSeed:
    
    split_start=index%split_size+index//split_size*split_size
    split_end=split_start+split_size
    if split_end>NSeed:
        split_end=NSeed
    
    for cmb in range(split_end-split_start):
        for ii in range(len(Fr[0])):
            if index==FR_cmb[ii]:
                fileName=folders['Scaling_factors']+'\SF_'+str(index+cmb).zfill(formats['fill_fn_all'])+'.out'
                w=np.loadtxt(fileName)
                SF[:,ii]=w[:,FR_p[ii]]
                Sa_un[ii,:]=Sa_unsc_ave_split[index,:]
#                print(index)
                
    index_spl=index_spl+1
    index=index+1



Sa_sc=Sa_un+np.expand_dims(np.mean(SF,axis=0), axis=1)



fig = plt.figure()
ax = plt.axes()

x = T[0]
y=Sa_Tgt[0]

plt.plot(x, y)

for i in range(len(Fr[0])):
    x = T[0]
    y=Sa_sc[i,:]
    plt.plot(x, y)













