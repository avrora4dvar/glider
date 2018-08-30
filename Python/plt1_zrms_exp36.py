# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 11:43:24 2017

Create map of study area

@author: ipasmans

"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import numpy as np
import scipy.io as sio
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import netCDF4
from os import listdir
import re as re
import mod_plotter as plotter
import datetime
import matplotlib.ticker as ticker


mat=sio.loadmat('V:\ipasmans\zprofile_plots_ensfor.mat',squeeze_me=True,struct_as_record=False)
mat=mat['plots']

#%% Load data

plt.close('all')
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})
fig=plt.figure(figsize=(5.61,8.92*.6))


ax=[]
for i in np.arange(0,6):
    ax1=fig.add_subplot(2,3,i+1)
    ax1.set_ylim(200,0)
    ax1.grid(color=(.7,.7,.7),linestyle=':',linewidth=.5)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(.5))
    if np.mod(i,3)==0:
        ax1.set_ylabel(r'Depth [m]')
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
    else:
        ax1.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.append(ax1)

#%% Plot

markers=['*','o','x',None]
colors=['m','r','g','b']
labels=['Surface','Glider','Combined','No DA']
linewidths=[1,1,1,1]

#Salinity
pplot=[]
val0=mat[2].val[:,-1]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[2].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths): 
    pplot1=ax[0].plot(val1,mat[2].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[0].text(.93,.05,'a)',transform=ax[0].transAxes)
ax[0].set_xlabel(r'Analysis $RMSE_S\,\mathrm{[ppt]}$')
ax[0].set_xlim(0,2.0)

#Temperature
pplot=[]
val0=mat[1].val[:,-1]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[1].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):    
    pplot1=ax[1].plot(val1,mat[1].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[1].text(.93,.05,'b)',transform=ax[1].transAxes)
ax[1].set_xlabel(u'Analysis $RMSE_T$ [\u00B0C]')
ax[1].set_xlim(0,2.5)

#Density
pplot=[]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[0].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[2].plot(val1,mat[0].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
    
ax[2].legend(loc='lower left',bbox_to_anchor=(.29,.3),framealpha=1)
ax[2].text(.93,.05,'c)',transform=ax[2].transAxes)
ax[2].set_xlabel(r'Analysis $RMSE_\rho\,\mathrm{[kgm^{-3}]}$')
ax[2].set_xlim(0,1.7)

#%% Forecasts

mat=sio.loadmat('V:\ipasmans\zprofile_plots_for.mat',squeeze_me=True,struct_as_record=False)
mat=mat['plots']

#%% Plot

markers=['*','o','x',None]
colors=['m','r','g','b']
labels=['Surface','Glider','Combined','No DA']
linewidths=[1,1,1,1]

#Salinity
pplot=[]
val0=mat[2].val[:,-1]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[2].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths): 
    pplot1=ax[3].plot(val1,mat[2].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[3].text(.93,.05,'d)',transform=ax[3].transAxes)
ax[3].set_xlabel(r'Forecast $RMSE_S\,\mathrm{[ppt]}$')
ax[3].set_xlim(0,2.0)

#Temperature
pplot=[]
val0=mat[1].val[:,-1]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[1].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[4].plot(val1,mat[1].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[4].text(.93,.05,'e)',transform=ax[4].transAxes)
ax[4].set_xlabel(u'Forecast $RMSE_T$ [\u00B0C]')
ax[4].set_xlim(0,2.5)

#Density
pplot=[]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[0].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[5].plot(val1,mat[0].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
    
#ax[5].legend(loc='lower left',bbox_to_anchor=(.29,.01),framealpha=1)
ax[5].text(.93,.05,'f)',transform=ax[5].transAxes)
ax[5].set_xlabel(r'Forecast $RMSE_\rho\,\mathrm{[kgm^{-3}]}$')
ax[5].set_xlim(0,1.7)


#%%

plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=.05)
#fig.subplots_adjust(top=.98,bottom=.09)

#fig.savefig('z_rms.pdf')

