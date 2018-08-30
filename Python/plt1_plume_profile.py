# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 14:06:57 2018

Plot 

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


mat=sio.loadmat('V:\ipasmans\TSz_plume.mat',squeeze_me=True,struct_as_record=False)

#%% Load data

plt.close('all')
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})
fig=plt.figure(figsize=(5.61*.9,8.92*.4))

#%% 

ax=[]
labels=['a)','b)']
for iax in np.arange(0,2):
    ax1=fig.add_subplot(1,2,iax+1)
    ax1.set_ylim(200,0)
    if iax==0:
        ax1.set_ylabel('Depth [m]')
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
    else:
        ax1.yaxis.set_major_formatter(ticker.NullFormatter())
    ax1.text(.02,.94,s=labels[iax],transform=ax1.transAxes)
    ax1.grid()
    ax.append(ax1)
        
#%% Plot salinity
    
ax1=ax[0]
ax1.set_xlabel('Salinity [ppt]')
ax1.plot(mat['cr'].salt[:,0],mat['cr'].depth,'k-',label='Obs (plume)')
ax1.plot(mat['cr'].salt[:,1],mat['cr'].depth,'b-',label='No DA (plume)')
ax1.plot(mat['nocr'].salt[:,0],mat['nocr'].depth,'k--',label='Obs (no plume)')
ax1.plot(mat['nocr'].salt[:,1],mat['nocr'].depth,'b--',label='No DA (no plume)')
leg=ax1.legend(loc='lower left')
leg.get_frame().set_alpha(1)

#%% Plot temperature
    
ax1=ax[1]
ax1.set_xlabel(u'Temperature [\u00B0C]')
ax1.plot(mat['cr'].temp[:,0],mat['cr'].depth,'k-',label='Obs (plume)')
ax1.plot(mat['cr'].temp[:,1],mat['cr'].depth,'b-',label='No DA (plume)')
ax1.plot(mat['nocr'].temp[:,0],mat['nocr'].depth,'k--',label='Obs (no plume)')
ax1.plot(mat['nocr'].temp[:,1],mat['nocr'].depth,'b--',label='No DA (no plume)')
leg=ax1.legend(loc='lower right')
leg.get_frame().set_alpha(1)

#%% Save

fig.savefig('profile_plume.pdf')