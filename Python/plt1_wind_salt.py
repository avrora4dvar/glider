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
import matplotlib
import seawater as sw
imp.reload(plotter)

plt.rcParams.update({'font.size':10})
matplotlib.rcParams['font.family']='Arial'

plt.close('all')
#fig=plt.figure(figsize=(5,3))
fig=plt.figure(figsize=(3.66,8.92*.5))

#%%

mat=sio.loadmat('V:\ipasmans\wind_stress.mat',squeeze_me=True,struct_as_record=False)
wind=mat['wind']


wind.t=wind.t-366.
dateRef=datetime(2005,1,1).toordinal()
tLim=np.array([2380,2392,2410,2434])+dateRef

#%% Plot

ax=[]
labels=['a)','b)','c)']
for iax in np.arange(0,2):
    ax1=fig.add_subplot(2,1,iax+1)
    ax1.grid(color=(.7,.7,.7),linestyle='--',linewidth=.5)
    ax1.xaxis.set_major_locator(ticker.FixedLocator(np.arange(wind.t[0],wind.t[-1]+.01,12)))
    if iax==1:
        ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:plotter.num2time(x).strftime('%m/%d')))
        ax1.set_xlabel('2011')
    else:
        ax1.xaxis.set_major_formatter(ticker.NullFormatter())
    ax1.text(.02,.93,s=labels[iax],transform=ax1.transAxes)
    ax1.set_xlim(tLim[0],tLim[-1])
    ax.append(ax1)

#%% Plot wind
iax=0
ax[iax].plot(wind.t,wind.v,linewidth=1,color='b')
ax[iax].plot([tLim[0],tLim[-1]],[0,0],color='k',linewidth=1)
ax[iax].plot([tLim[1],tLim[1]],[-10,10],color='k',linewidth=1,linestyle='--')
ax[iax].plot([tLim[2],tLim[2]],[-10,10],color='k',linewidth=1,linestyle='--')

ax[iax].text(tLim[0]+.5,8,'{0}\n{1}'.format('down-','welling'))
ax[iax].text(tLim[0]+.5,-9,'upwelling')

ax[iax].yaxis.set_major_locator(ticker.MultipleLocator(.05))
ax[iax].set_ylim(-.15,.15)
ax[iax].set_ylabel('{0} {1}'.format(r'Meridional wind stress',r'[$\mathrm{Nm^{-2}}$]'))

#%% Plot No DA, observations salinity NH10

mat=sio.loadmat('V:\ipasmans\\NH10_ana.mat',squeeze_me=True,struct_as_record=False)
mat['obs'].t=mat['obs'].t-366.

#Apply running mean
mat['obs'].saltf=plotter.running_mean(mat['obs'].t,mat['obs'].salt[0,:],1.)
for mod1 in mat['model']:
    mod1.saltf=plotter.running_mean(mat['obs'].t,mod1.salt[0,:],1.)

ax1=ax[1]
ax1.plot(mat['obs'].t,mat['obs'].saltf,'k-',label='Obs')
ax1.plot(mat['obs'].t,mat['model'][3].saltf,'b-',label='No DA')
ax1.set_ylabel('Salinity [ppt]')
leg=ax1.legend(loc='lower right')
leg.get_frame().set_alpha(1)

#%%


fig.subplots_adjust(left=.2,bottom=.1,top=.98,right=.96,hspace=.08)

#fig.savefig('wind.png',figsize=(5,2.5),dpi=500)
fig.savefig('wind.pdf')