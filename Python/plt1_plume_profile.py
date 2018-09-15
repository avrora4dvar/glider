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


mat=sio.loadmat('E:\ipasmans\TS_plume_exp26.mat',squeeze_me=True,struct_as_record=False)

#%% Load data

plt.close('all')
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})
fig=plt.figure(figsize=(5.61,8.92*.39))

#%% 

ax=[]
labels=['a)','b)']
for iax in np.arange(0,2):
    ax1=fig.add_subplot(1,2,iax+1)
    ax1.set_ylim(150,0)
    if iax==0:
        ax1.set_ylabel('Depth [m]')
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(25))
    else:
        ax1.yaxis.set_major_formatter(ticker.NullFormatter())
    ax1.text(.02,.94,s=labels[iax],transform=ax1.transAxes)
    ax1.grid()
    ax.append(ax1)
        
#%% Plot salinity
    
ax1=ax[0]
ax1.set_xlabel('Salinity')
ax1.plot(mat['cr'].salt[:,0],mat['cr'].depth,'k-',label='Obs (plume)')
ax1.plot(mat['cr'].salt[:,1],mat['cr'].depth,'-',color='dodgerblue',label='Long Free (plume)')
ax1.plot(mat['nocr'].salt[:,0],mat['nocr'].depth,'k--',label='Obs (no plume)')
ax1.plot(mat['nocr'].salt[:,1],mat['nocr'].depth,'--',color='dodgerblue',label='Long Free (no plume)')
leg=ax1.legend(loc='lower left',bbox_to_anchor=(-.05,.02))
leg.get_frame().set_alpha(1)

#%% Plot temperature
    
ax1=ax[1]
ax1.set_xlabel(u'Temperature [\u00B0C]')
ax1.plot(mat['cr'].temp[:,0],mat['cr'].depth,'k-',label='Obs (plume)')
ax1.plot(mat['cr'].temp[:,1],mat['cr'].depth,'-',color='dodgerblue',label='Long Free (plume)')
ax1.plot(mat['nocr'].temp[:,0],mat['nocr'].depth,'k--',label='Obs (no plume)')
ax1.plot(mat['nocr'].temp[:,1],mat['nocr'].depth,'--',color='dodgerblue',label='Long Free (no plume)')
leg=ax1.legend(loc='lower left',bbox_to_anchor=(.3,.02))
leg.get_frame().set_alpha(1)

#%% Save

fig.subplots_adjust(top=.98,left=.1,right=.98)
fig.savefig('profile_plume.pdf')