# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:08:06 2018

@author: ipasmans
"""

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
from matplotlib import patches
imp.reload(plotter)

plt.close('all')
fig=plt.figure(figsize=(3.66,8.92*.3))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})

#%% Read grid

grd=plotter.read_grid_rho()

#%% Generate plot

ax=[]
for iax in np.arange(0,3):
    ax1=fig.add_subplot(3,1,iax+1)
    ax1.set_axis_off()

    
    ax1.set_xlim(-124.8,-124)
    ax1.set_ylim(100,-20)
    
    ax1.add_patch(patches.Rectangle((-124.4,0),.25,120,fill=False,linestyle='--',linewidth=1.5,color='gray'))

    ax.append(ax1)
    
#%% Different 
    
lon1=np.arange(-124.8,-124.0001,.01)

#No assimilation
val1=15*np.tanh((lon1+124.1)/.1)
val1[val1>3]=np.NaN
ax[0].plot(lon1,val1,linewidth=1.5,color='r')
val1=10-15*np.tanh((lon1+124.2)/.05)
val1[val1<3]=np.NaN
ax[0].plot(lon1,val1,linewidth=1.5,color='b')
val1=35-15*np.tanh((lon1+124.2)/.1)
val1[val1<3]=np.NaN
ax[0].plot(lon1,val1,linewidth=1.5,color='b')
ax[0].text(-124.78,90.3,'a)')
    
#Glider
w=.5+.5*np.tanh((lon1+124.4)/.02)
val1=7*np.tanh((lon1+124.1)/.1)*w
val1=val1+(1-w)*15*np.tanh((lon1+124.1)/.1)
val1[val1>3]=np.NaN
ax[1].plot(lon1,val1,linewidth=1.5,color='r')
val1=w*(5-7*np.tanh((lon1+124.3)/.05))
val1=val1+(1-w)*(10-15*np.tanh((lon1+124.2)/.05))
val1[val1<3]=np.NaN
ax[1].plot(lon1,val1,linewidth=1.5,color='b')
val1=w*(15-18*np.tanh((lon1+124.2)/.1))
val1=val1+(1-w)*(35-15*np.tanh((lon1+124.2)/.1))
val1[val1<3]=np.NaN
ax[1].plot(lon1,val1,linewidth=1.5,color='b')
ax[1].text(-124.78,90.3,'b)')

    
#Combined
val1=7*np.tanh((lon1+124.1)/.1)
val1[val1>3]=np.NaN
ax[2].plot(lon1,val1,linewidth=1.5,color='r')
val1=5-7*np.tanh((lon1+124.3)/.05)
val1[val1<3]=np.NaN
ax[2].plot(lon1,val1,linewidth=1.5,color='b')
val1=15-18*np.tanh((lon1+124.2)/.1)
val1[val1<3]=np.NaN
ax[2].plot(lon1,val1,linewidth=1.5,color='b')
ax[2].text(-124.78,90.3,'c)')


#%%

for ax1 in ax:
    #Coast   
    h=grd['h'][:,223]; h[np.isnan(h)]=3.
    h=np.concatenate(([10e3],h,[10e3]))
    hlon=np.concatenate(([-135.],grd['lon'][:,0],[-120.]))
    p=patch.Polygon(np.column_stack((hlon,h)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
    p.set_zorder(20)
    ax1.add_patch(p) 

#%% Save

fig.subplots_adjust(bottom=.02,top=.98)
fig.savefig('hypo.pdf')
