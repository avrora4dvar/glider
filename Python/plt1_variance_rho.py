# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 09:49:47 2018

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
imp.reload(plotter)

#%% Figure

plt.rcParams.update({'font.size':8})
matplotlib.rcParams['font.family']='Arial'
plt.close('all')
fig=plt.figure(figsize=(5.61,8.92*.35))

dateRef=plotter.time2num(datetime(2005,1,1))
tLim=np.array([2392,2413])+dateRef

ax=[]
for iax in np.arange(0,3):
    ax1=fig.add_subplot(1,3,iax+1)
    
    #xaxis
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
    ax1.set_xlabel('Longitude')
    ax1.set_xlim(-129,-123.5)
    
    #yaxis
    ax1.set_ylim(0,60)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(5))
    if iax==0:
        ax1.set_ylabel('Std. isopycnal depth [m]')
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
    else:
        ax1.yaxis.set_major_formatter(ticker.NullFormatter())
        None
    
    ax1.grid(color=(.5,.5,.5),linestyle='--',linewidth=.5)
    
    ax.append(ax1)
    

#%%

#mat=sio.loadmat('V:\ipasmans\cross_NH10_for.mat',squeeze_me=True,struct_as_record=False)
mat=sio.loadmat('V:\ipasmans\cross_NH10_his_zrho_ana.mat',squeeze_me=True,struct_as_record=False)
loc=mat['loc']; mat=mat['model']
mat=mat[[0,1,2,4]]

for mat1 in mat:
    mat1.t=mat1.t-366.

#Read observations
obs=plotter.read_obslist()
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6]
lon1,lat1=zip(*obsll)




#%% Extract contours

markers=['*','o','x',None]
colors=['m','r','g','b']
labels=['Surface','Glider','Combined','No DA']
linewidths=[1,1,1,1]

incon=mat[0].rho==1024.5
ax[0].set_title(r'a) $\sigma$=24.5 $\mathrm{kg m^{-3}}$')
plots=[]
for (mat1,marker1,color1,label1) in zip(mat,markers,colors,labels):
    zstd=np.nanstd(mat1.z,axis=2)
    plt1=ax[0].plot(mat1.lon[:,0],zstd[:,incon],marker=marker1,color=color1,markevery=10,label=label1)
    plots.append(plt1)
gplot1=ax[0].plot([np.min(lon1),np.min(lon1)],[0,100],'k--')

incon=mat[0].rho==1025.5
ax[1].set_title(r'b) $\sigma$=25.5 $\mathrm{kg m^{-3}}$')
plots=[]
for (mat1,marker1,color1,label1) in zip(mat,markers,colors,labels):
    zstd=np.nanstd(mat1.z,axis=2)
    plt1=ax[1].plot(mat1.lon[:,0],zstd[:,incon],marker=marker1,color=color1,markevery=10,label=label1)
    plots.append(plt1)
gplot1=ax[1].plot([np.min(lon1),np.min(lon1)],[0,100],'k--')

incon=mat[0].rho==1026.5
ax[2].set_title(r'c) $\sigma$=26.5 $\mathrm{kg m^{-3}}$')
plots=[]
for (mat1,marker1,color1,label1) in zip(mat,markers,colors,labels):
    zstd=np.nanstd(mat1.z,axis=2)
    plt1=ax[2].plot(mat1.lon[:,0],zstd[:,incon],marker=marker1,color=color1,markevery=10,label=label1)
    plots.append(plt1)
gplot1=ax[2].plot([np.min(lon1),np.min(lon1)],[0,100],'k--')

leg=ax[0].legend(loc='upper left')
leg.get_frame().set_alpha(1)

#%% Save

fig.subplots_adjust(top=.92,bottom=.14,left=.08,right=.97)
fig.savefig('cross_var_zrho.pdf')



