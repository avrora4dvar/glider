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

plt.rcParams.update({'font.size':12})

#%% Load data

mat=sio.loadmat('V:\ipasmans\zprofile_plots_for.mat',squeeze_me=True,struct_as_record=False)
mat=mat['plots'][0]

plt.close('all')
fig=plt.figure(figsize=(3,4.5),dpi=300)
plt.rcParams.update({'font.size': 12})

ax=fig.add_subplot(1,1,1)
ax.set_ylim(200,0)
ax.grid(color=(.7,.7,.7),linestyle=':',linewidth=.5)
ax.set_ylabel(r'Depth [m]')

fig.subplots_adjust(left=.22,bottom=.16,top=.98)

#%% Plot

markers=['*','o','x',None]
colors=['m','r','g','b']
labels=['Surface','Glider','Combined','No DA']
linewidths=[1.5,3.5,3.5,1.5]

pplot=[]
for (val1,label1,marker1,color1,linewidth1) in zip(mat.val[:,[1,2,3,-1]].T,labels,markers,colors,linewidths):
    print(label1)
    pplot1=ax.plot(val1,mat.z,color=color1,marker=marker1,label=label1,markevery=6,linewidth=linewidth1)
    pplot.append(pplot1)

#plt.legend(loc='upper left',bbox_to_anchor=(.01,.99),framealpha=1)
plt.legend(loc='lower left',bbox_to_anchor=(.3,.01),framealpha=1)
ax.set_xlim(0,1.5)
#ax.set_xlabel(r'$RMSE_S\;\mathrm{[ppt]}$')
#ax.set_xlabel(r'$RMSE_T\;\mathrm{[^\circ C]}$')
ax.set_xlabel(r'$RMSE_\rho\;\mathrm{[kgm^{-3}]}$')
#ax.set_xlabel(r'$skill_\rho$')

#fig.savefig('rms_z_for_rho.png',dpi=600)

