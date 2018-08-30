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


mat=sio.loadmat('V:\ipasmans\zprofile_plots_exp384142.mat',squeeze_me=True,struct_as_record=False)
mat=mat['plots']


#%% Load data

plt.close('all')
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})
fig=plt.figure(figsize=(5.61,8.92*.3))


ax=[]
for i in np.arange(0,3):
    ax1=fig.add_subplot(1,3,i+1)
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
colors=['r','b','g','k']
labels=['Glider Only','Glider T','Glider+SST','Glider+HFR']
linewidths=[1,1,1,1]

#Salinity
pplot=[]
val0=mat[2].val[:,2]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[2].val[:,[1,3,5,6]].T,labels,markers,colors,linewidths): 
    pplot1=ax[0].plot(val1/val0,mat[2].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[0].text(.03,.95,'a)',transform=ax[0].transAxes)
ax[0].set_xlabel(r'Glider salinity')
ax[0].set_xlim(0,3)

#Temperature
pplot=[]
val0=mat[1].val[:,2]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[1].val[:,[1,3,5,6]].T,labels,markers,colors,linewidths):    
    pplot1=ax[1].plot(val1/val0,mat[1].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
ax[1].text(.03,.95,'b)',transform=ax[1].transAxes)
ax[1].set_xlabel(u'Glider temperature')
ax[1].set_xlim(0,2)

#Density
pplot=[]
mat0=mat[0].val[:,2]
for (val1,label1,marker1,color1,linewidth1) in zip(mat[0].val[:,[1,3,5,6]].T,labels,markers,colors,linewidths):
    pplot1=ax[2].plot(val1/val0,mat[0].z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)
    
ax[2].legend(loc='lower left',bbox_to_anchor=(.29,.3),framealpha=1)
ax[2].text(.03,.95,'c)',transform=ax[2].transAxes)
ax[2].set_xlabel(r'Glider density')
ax[2].set_xlim(0,2.1)


#%%

plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=.05)
#fig.subplots_adjust(top=.98,bottom=.09)

fig.savefig('z_rms_ratio.pdf')

