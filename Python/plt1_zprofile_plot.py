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


mat=sio.loadmat('V:\ipasmans\zprofile_plots_ana.mat',squeeze_me=True,struct_as_record=False)
mat=mat['plots']

#%% Load data

plt.close('all')
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})
fig=plt.figure(figsize=(5.61,8.92*.58))


ax=[]
for i in np.arange(0,3):
    ax1=fig.add_subplot(2,3,i+1)
    ax1.set_ylim(200,0)
    ax1.grid(color=(.7,.7,.7),linestyle=':',linewidth=.5)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(.5))
    if i==0:
        ax1.set_ylabel(r'Depth [m]')
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
    else:
        ax1.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.append(ax1)

for i in np.arange(0,3):
    ax1=fig.add_subplot(2,3,i+4)
    ax1.set_ylim(200,0)
    ax1.grid(color=(.7,.7,.7),linestyle=':',linewidth=.5)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(.2))
    ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.1f}'))
    if i==0:
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
for (val1,label1,marker1,color1,linewidth1) in zip(mat[2].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[0].plot(val1,mat[2].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[0].text(.1,195,'a)')
ax[0].set_xlabel(r'$RMSE_S\,\mathrm{[ppt]}$')
ax[0].set_xlim(0,2.0)

#Temperature
pplot=[]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[1].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[1].plot(val1,mat[1].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[1].text(.1,195,'b)')
ax[1].set_xlabel(u'$RMSE_T$ [\u00B0C]')
ax[1].set_xlim(0,2.5)

#Density
pplot=[]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[0].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[2].plot(val1,mat[0].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
    
ax[2].legend(loc='lower left',bbox_to_anchor=(.29,.01),framealpha=1)
ax[2].text(.1,195,'c)')
ax[2].set_xlabel(r'$RMSE_\rho\,\mathrm{[kgm^{-3}]}$')
ax[2].set_xlim(0,1.7)

#%% Plot

markers=['*','o','x',None]
colors=['m','r','g','b']
labels=['Surface','Glider','Combined','No DA']
linewidths=[1,1,1,1]

#Salinity
pplot=[]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[5].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[5].plot(val1,mat[5].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[3].text(.1,195,'d)')
ax[3].set_xlabel(r'Skill$_S$')
ax[3].set_xlim(0,1)

#Temperature
for (val1,label1,marker1,color1,linewidth1) in zip(mat[4].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[4].plot(val1,mat[4].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
ax[4].text(.1,195,'e)')
ax[4].set_xlabel(r'Skill$_T$')
ax[4].set_xlim(0,1)

#Density
for (val1,label1,marker1,color1,linewidth1) in zip(mat[3].val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    pplot1=ax[3].plot(val1,mat[3].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    
ax[5].text(.1,195,'f)')
ax[5].set_xlabel(r'Skill$_\rho$')
ax[5].set_xlim(0,1)

#%%

fig.subplots_adjust(top=.98,bottom=.09)
fig.savefig('z_rms_skill_ana.pdf')

