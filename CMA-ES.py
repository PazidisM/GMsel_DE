# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 11:11:46 2018

@author: Marios
"""
import numpy as np
from Cost_functions import CF_0
from math import log
#from Cost_functions import CF_1

def cigar(x):
    f = x[1]**2+ 1e6*sum(x[1::]**2)
    return f





#--------------------------------------------------
# User defined input parameters 
#--------------------------------------------------

N=10    # number of objective variables/problem dimension
xmeanw = np.random.rand(N,1)     # objective variables initial point
sigma=1    # coordinate wise standard deviation (step size)
minsigma= 1e-15     # minimal step size
stopfitness = 1e-10    # stop if fitness < stopfitness (minimization)
maxeval = 300*(N+2)**2    # stop after stopeval number of function 


#--------------------------------------------------
# Strategy parameter setting: Selection
#--------------------------------------------------

lambd = 4+int(3*log(N))    # population size, offspring number
mu = int(lambd/2)    # number of parents/points for recombination
arweights =np.repeat(log((lambd+1)/2),mu)-np.array([log(i) for i in range(1,mu+1)])   # muXone array for weighted recombination

#--------------------------------------------------
# Strategy parameter setting: Adaptation
#--------------------------------------------------

cc = 4/(N+4)    # time constant for cumulation for C
cs = 4/(N+4)   # t-const for cumulation for sigma control
ccov= 2/(N+2**0.5)**2
damp= 1/cs+ 1

#--------------------------------------------------
# Initialize dynamic strategy parameters and constants
#--------------------------------------------------

pc = np.zeros(N)    # evolution paths for C
ps = np.zeros(N)    # evolution paths for sigma
B = np.eye(N)    # B defines the coordinate system
D =np.eye(N)    # diagonal D defines the scaling
BD=B*D
C = BD*np.transpose(BD)    # covariance matrix C
cw=np.sum(arweights)/np.linalg.norm(arweights)
chiN= N**0.5*(1-1/(4*N)+1/(21*N**2))

#--------------------------------------------------
# Generation loop
#--------------------------------------------------

arfitness=np.empty(lambd)
arfitness[0]=2*abs(stopfitness)+1
arz=np.empty((N,lambd))
arx=np.empty((N,lambd))


# Generate and evaluate lambd offspring
counteval= 0
while arfitness[0]>stopfitness and counteval< maxeval:
    for k in range(lambd):
        arz[:,[k]]=np.random.randn(N,1)
        arx[:,[k]]=xmeanw.reshape(N,1)+sigma*(BD @ arz[:,[k]])
#        arfitness[k]=CF_0(Sa_ave,Sa_Tgt,arz[:,k])
        arfitness[k]=cigar(arz[:,k])
        counteval= counteval+1
    
    # Sort by fitness and compute weighted mean
    arindex=np.argsort(arfitness)
    xmeanw= arx[:,arindex[range(mu)]] @ arweights/sum(arweights)
    zmeanw= arz[:,arindex[range(mu)]] @ arweights/sum(arweights)
    
    # Adapt covariance matrix
    pc = (1-cc)*pc+ ((cc*(2-cc))**0.5*cw)*(BD @ zmeanw)
    C = (1-ccov)*C+ ccov*pc @ np.transpose(pc)
    
    # adapt sigma
    ps = (1-cs)*ps+ ((cs*(2-cs))**0.5*cw)* (B @ zmeanw)
    sigma= sigma* np.exp((np.linalg.norm(ps)-chiN)/chiN/damp)
    
    # Update B and D from C
    if (counteval/lambd)%(N/10):
        C=np.triu(C)+np.transpose(np.triu(C,1)) # enforce symmetry
        B, D = np.linalg.eig(C)
        
        # limit condition of C to 1e14+1
        if max(np.diag(D))>1e14*min(np.diag(D)):
            tmp=max(np.diag(D))/1e14-min(np.diag(D))
            C=C+tmp * np.eye(N)
            D=D+tmp * np.eye(N)
        B=np.diag(np.diag(B)**0.5)      # D contains standard deviations now
        BD=D*B
    # Adjust minimal step size
    if sigma*np.min(D)<minsigma or arfitness[0]==arfitness[min(mu+1,lambd)-1] or xmeanw==xmeanw+0.2*sigma*BD[:,1+int((counteval/lambd)%N)]:
        sigma=1.4*sigma


print(counteval)


foo=np.transpose(xmeanw)+sigma*(BD @ arz[:,[k]])
foo=xmeanw+sigma*(BD @ arz[:,[k]])
xmeanw.reshape(N,1)

foo2=BD @ arz[:,[k]]
