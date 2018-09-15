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
from matplotlib import dates as dates
import matplotlib
import seawater as sw
imp.reload(plotter)

plt.rcParams.update({'font.size':10})
matplotlib.rcParams['font.family']='Arial'

plt.close('all')
#fig=plt.figure(figsize=(5,3))
fig=plt.figure(figsize=(5.61,8.92*.5))

#%%

mat=sio.loadmat('E:\ipasmans\wind_stress.mat',squeeze_me=True,struct_as_record=False)
wind=mat['wind']


wind.t=wind.t-366.
dateRef=datetime(2005,1,1).toordinal()
tLim=np.array([2380,2392,2413,2434])+dateRef

#%% Plot

ax=[]
labels=['a)','b)','c)']
for iax in np.arange(0,2):
    ax1=fig.add_subplot(2,1,iax+1)
    ax1.grid(color=(.7,.7,.7),linestyle='--',linewidth=.5)
    ax1.xaxis.set_major_locator(dates.MonthLocator())
    if iax==1:
        ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:plotter.num2time(x).strftime('%m/%d')))
        ax1.set_xlabel('2011')
    else:
        ax1.xaxis.set_major_formatter(ticker.NullFormatter())
    ax1.text(.02,.92,s=labels[iax],transform=ax1.transAxes)
    ax1.set_ylim(24.,34.)
    ax1.set_xlim(datetime(2011,5,1),datetime(2011,10,1))
    ax.append(ax1)

#%% Plot wind
iax=0
#wind.vf=plotter.running_mean(wind.t,wind.v,1.)
#wind.vf=plotter.running_mean(wind.t,wind.vf,1.)
ax[iax].plot(wind.t,wind.v,linewidth=1,color='b')
ax[iax].plot([datetime(2011,1,1),datetime(2011,12,1)],[0,0],color='k',linewidth=1)
ax[iax].plot([tLim[1],tLim[1]],[-10,10],color='k',linewidth=1,linestyle='--')
ax[iax].plot([tLim[2],tLim[2]],[-10,10],color='k',linewidth=1,linestyle='--')

ax[iax].text(.7,.9,s='{0}'.format('downwelling'),transform=ax[iax].transAxes)
ax[iax].text(.7,.04,'upwelling',transform=ax[iax].transAxes)

ax[iax].yaxis.set_major_locator(ticker.MultipleLocator(.05))
ax[iax].set_ylim(-.15,.15)
ax[iax].set_ylabel('{0} {1}'.format(r'Meridional wind stress',r'[$\mathrm{Nm^{-2}}$]'))

#%% Plot No DA, observations salinity NH10

mat=sio.loadmat('E:\ipasmans\\NH10_exp26.mat',squeeze_me=True,struct_as_record=False)
mat['obs'].t=mat['obs'].t-366.

#Apply running mean
mat['obs'].saltf=plotter.running_mean(mat['obs'].t,mat['obs'].salt[0,:],1.)
mat['model'].saltf=plotter.running_mean(mat['model'].t,mat['model'].salt[0,:],1.)
mat['obs'].saltf=plotter.running_mean(mat['obs'].t,mat['obs'].saltf,1.)
mat['model'].saltf=plotter.running_mean(mat['model'].t,mat['model'].saltf,1.)

ax1=ax[1]
ax1.plot(mat['obs'].t,mat['obs'].saltf,'k-',label='NH10')
ax1.plot(mat['obs'].t,mat['model'].saltf,'-',color='dodgerblue',label='Long Free')
ax1.plot([tLim[1],tLim[1]],[-10,40],color='k',linewidth=1,linestyle='--')
ax1.plot([tLim[2],tLim[2]],[-10,40],color='k',linewidth=1,linestyle='--')
ax1.set_ylim(21.5,34)
ax1.set_ylabel('Salinity')
leg=ax1.legend(loc='lower left',bbox_to_anchor=(.71,.01))
leg.get_frame().set_alpha(1)

#%%


fig.subplots_adjust(left=.15,bottom=.1,top=.98,right=.95,hspace=.08)

#fig.savefig('wind.png',figsize=(5,2.5),dpi=500)
fig.savefig('wind.pdf')