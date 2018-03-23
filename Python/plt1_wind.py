# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 17:36:37 2018

@author: ipasmans
"""

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

plt.rcParams.update({'font.size':8})
matplotlib.rcParams['font.family']='Arial'

plt.close('all')
#fig=plt.figure(figsize=(5,3))
fig=plt.figure(figsize=(3.66,8.92*.35))

#%%

mat=sio.loadmat('V:\ipasmans\wind_stress.mat',squeeze_me=True,struct_as_record=False)
wind=mat['wind']

wind.t=wind.t-366.
dateRef=datetime(2005,1,1).toordinal()
tLim=np.array([2380,2392,2410,2434])+dateRef

#%% Plot



ax=fig.add_subplot(1,1,1)

ax.plot(wind.t,wind.v,linewidth=1,color='b')
ax.plot([tLim[0],tLim[-1]],[0,0],color='k',linewidth=1)
ax.plot([tLim[1],tLim[1]],[-10,10],color='k',linewidth=1,linestyle='--')
ax.plot([tLim[2],tLim[2]],[-10,10],color='k',linewidth=1,linestyle='--')

ax.text(tLim[0]+.5,8,'{0}\n{1}'.format('down-','welling'))
ax.text(tLim[0]+.5,-9,'upwelling')

ax.grid(color=(.7,.7,.7),linestyle='--',linewidth=.5)
#ax.spines['bottom'].set_linewidth(.5)
#ax.spines['top'].set_linewidth(.5)
#ax.spines['left'].set_linewidth(.5)
#ax.spines['right'].set_linewidth(.5)

ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(wind.t[0],wind.t[-1]+.01,12)))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:plotter.num2time(x).strftime('%m/%d')))
ax.yaxis.set_major_locator(ticker.MultipleLocator(.05))
ax.set_ylim(-.15,.15)
ax.set_xlim(tLim[0],tLim[-1])

ax.set_xlabel('2011')
ax.set_ylabel('{0} {1}'.format(r'Meridional wind stress',r'velocity [$\mathrm{Nm^{-2}}$]'))
fig.subplots_adjust(left=.2,bottom=.14,top=.97,right=.96)

#%%

#fig.savefig('wind.png',dpi=600)
fig.savefig('wind.pdf')