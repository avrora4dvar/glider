# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 13:51:11 2018

Calculate BV-frequency

@author: ipasmans
"""


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import numpy as np
import scipy.io as sio
import mod_plotter as plotter
import netCDF4
import scipy.io as sio
from datetime import datetime, timedelta
import imp
from matplotlib import animation
from matplotlib import ticker
import seawater as sw
from matplotlib import patches
import gsw
imp.reload(plotter)

plt.close('all')
fig=plt.figure(figsize=(7.4,8.2))
#fig=plt.figure(figsize=(9,6))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':10})

#%% Load

mat=sio.loadmat('V:/ipasmans/along_glider.mat',squeeze_me=True,struct_as_record=False)

#Time to Python
mat['obs'].t=mat['obs'].t-366
    
def dateFormatter(x,pos):
    date=plotter.num2time(x).strftime('%m/%d')
    return date

tLim=plotter.time2num([datetime(2011,7,21),datetime(2011,8,11)])

#%% 

def BV2(S,T,depth,lon,lat):
    p1=gsw.p_from_z(depth,np.array(lat))
    sa=gsw.SA_from_SP(S,p1,lon,lat)
    ct=gsw.CT_from_pt(sa,T)
    N2,pOut=gsw.Nsquared(sa,ct,p1,lat)
    return N2,pOut
    
mat1=mat['model'][3]
N2=None; p=None
for (salt1,temp1,lon1,lat1) in zip(mat1.salt,mat1.temp,mat['obs'].lon,mat['obs'].lat):
    depth1=-np.array(mat['obs'].depth,dtype=float)
    N1,p1=BV2(salt1,temp1,depth1,lon1,lat1)
    if N2 is None:
        N2=N1; p=p1
    else:
        N2=np.vstack((N2,N1))
        p=np.vstack((p,p1))
        
    
z=np.reshape(mat['obs'].depth[:-1],(1,-1))+np.zeros(np.shape(N2))
t=np.reshape(mat['obs'].t,(-1,1))+np.zeros(np.shape(N2))


fig=plt.figure()
ax=fig.add_subplot(1,1,1)
n_instable=0
for N in N2:
    if np.any(N<0):
        ax.plot(.5*mat['obs'].depth[1:]+.5*mat['obs'].depth[:-1],N)
        n_instable=n_instable+1
print(n_instable)