# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:34:36 2018

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

plt.close('all')
fig=plt.figure(figsize=(3.66,8.92*.35))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})


#%% Read

mat=sio.loadmat('V:\ipasmans\cov_z.mat',squeeze_me=True,struct_as_record=False)

#%% Plot

ax1=fig.add_subplot(111)
ax2=ax1.twiny()

ax1.plot(mat['profile'].sig_temp,mat['profile'].z,'k:',linewidth=1)
ax1.plot(mat['c1'],-mat['z1'],'k-',linewidth=1)
ax1.plot(mat['c2'],-mat['z1'],'k--',linewidth=1)

ax2.plot(mat['profile'].temp,mat['profile'].z,linestyle='-',color='gray',linewidth=1)

#Layout
ax1.set_ylim(250,0); ax1.set_xlim(-.05,.9)
ax2.set_ylim(250,0); ax2.set_xlim(7-.38888889,14)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(.1))
ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.1f}'))
ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
ax1.yaxis.grid(color='gray',linestyle='--',linewidth=.5)
ax1.xaxis.grid(color='gray',linestyle='--',linewidth=.5)
ax2.xaxis.set_major_locator(ticker.MultipleLocator(.777777))
ax2.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.1f}'))
ax2.xaxis.label.set_color('gray')
ax2.xaxis.set_tick_params(labelcolor='gray')
ax1.set_ylabel(r'Depth [m]')
ax1.set_xlabel(u'Temperature,temperature-covariance [\u00B0C]')
ax2.set_xlabel(u'Average glider temperature [\u00B0C]')

#%% Save

fig.subplots_adjust(bottom=.15,top=.88)
fig.savefig('cov_z.pdf')