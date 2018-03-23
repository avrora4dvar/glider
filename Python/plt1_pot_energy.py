# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:25:07 2018

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
import matplotlib
imp.reload(plotter)

#%% 

t=plotter.time2num(datetime(2011,8,8))

plt.close('all')
fig=plt.figure(figsize=(5.61,8.92*.65))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})

#Obslist
obs=plotter.read_obslist()
#Glider location
tLim1=3.*int((t-2)/3)+np.array([2.,5.])
obsll=[(lon1,lat1,z1) for (lon1,lat1,type1,t1,z1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t'],obs['z']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
lonObs,latObs,zObs=zip(*obsll)

mat=sio.loadmat('V:/ipasmans/Epot_2410.mat',squeeze_me=True,struct_as_record=False)

#%% 

ax=[]
ax.append(plt.subplot2grid((2, 2), (0, 0),rowspan=2))
ax.append(plt.subplot2grid((2, 2), (0, 1),colspan=1))
ax.append(plt.subplot2grid((2, 2), (1, 1),colspan=1))

cm=plt.cm.get_cmap('seismic')

#%% Surface


ax[0]=plotter.add_bathy(ax[0])
ax[0]=plotter.add_coast(ax[0])
ax[0]=plotter.layout_surface(ax[0])

val1=mat['surface'][1]-mat['surface'][0]
print([np.nanpercentile(val1,.5),np.nanpercentile(val1,99.5)])
levels1=[tick1 for tick1 in np.arange(-75e3,75.001e3,12.5e3) if np.abs(tick1)>.001e3]
cplot=ax[0].contourf(mat['grd'].lon,mat['grd'].lat,val1,levels=levels1,cmap=cm,extend='both')
cplot=ax[0].contourf(mat['grd'].lon,mat['grd'].lat,val1,levels=levels1,cmap=cm,extend='both')

ax[0].text(-129.9,41.1,'a)')
ax[0].set_xlabel('Longitude')
ax[0].set_ylabel('Latitude')

gplot1=ax[0].plot(lonObs,latObs,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')

cbar=plt.colorbar(cplot,ax=ax[0],orientation='horizontal')
#cbar.set_ticks(levels1)
cbar.locator=ticker.MultipleLocator(25e3)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x*1e-4))
cbar.set_label(r'Pot. energy density difference [$10^4\:\mathrm{J m^{-2}}$]')
cbar.update_ticks()

#%% Cross Z


val1=mat['crossZ'][1]-mat['crossZ'][0]
print([np.nanpercentile(val1,.5),np.nanpercentile(val1,99.5)])
levels1=[tick1 for tick1 in np.arange(-40e6,40.001e6,5e6) if np.abs(tick1)>.001e6]
lon1=np.reshape(mat['grd'].lon[:,0],(-1,1))+np.zeros(np.shape(val1))
z1=np.reshape(mat['z'],(1,-1))+np.zeros(np.shape(val1))
cplot=ax[1].contourf(lon1,-z1,val1,levels=levels1,cmap=cm,extend='both')
cplot=ax[1].contourf(lon1,-z1,val1,levels=levels1,cmap=cm,extend='both')

ax[1].text(-129.9,385,'b)')
ax[1].set_ylim(400,0)
ax[1].set_xlim(-130,-123.8)
ax[1].set_xlabel('Longitude')
ax[1].set_ylabel('Depth [m]')

gplot1=ax[1].plot([np.min(lonObs),np.min(lonObs),np.max(lonObs),np.max(lonObs)],
    [np.min(zObs),np.max(zObs),np.max(zObs),np.min(zObs)],':',color=(.5,.5,0.5),linewidth=1.4)

cbar=plt.colorbar(cplot,ax=ax[1],orientation='horizontal')
#cbar.set_ticks(levels1)
cbar.locator=ticker.MultipleLocator(20e6)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x*1e-7))
cbar.set_label(r'Pot. energy density difference [$10^7\:\mathrm{J m^{-2}}$]')
cbar.update_ticks()

#%% Cross M


val1=mat['crossM'][1]-mat['crossM'][0]
print([np.nanpercentile(val1,.5),np.nanpercentile(val1,99.5)])
levels1=[tick1 for tick1 in np.arange(-28e6,28.001e6,4e6) if np.abs(tick1)>.001e6]
lat1=np.reshape(mat['grd'].lat[0,:],(1,-1)).T+np.zeros(np.shape(val1))
z1=np.reshape(mat['z'],(1,-1))+np.zeros(np.shape(val1))
cplot=ax[2].contourf(lat1,-z1,val1,levels=levels1,cmap=cm,extend='both')
cplot=ax[2].contourf(lat1,-z1,val1,levels=levels1,cmap=cm,extend='both')

ax[2].text(41.1,385,'c)')
ax[2].set_ylim(400,0)
ax[2].set_xlim(41),50
ax[2].set_xlabel('Latitude')
ax[2].set_ylabel('Depth [m]')

gplot1=ax[2].plot([np.min(latObs),np.min(latObs),np.max(latObs),np.max(latObs)],
    [np.min(zObs),np.max(zObs),np.max(zObs),np.min(zObs)],':',color=(.5,.5,0.5),linewidth=1.4)

cbar=plt.colorbar(cplot,ax=ax[2],orientation='horizontal')
#cbar.set_ticks(levels1)e
cbar.locator=ticker.MultipleLocator(8e6)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x*1e-7))
cbar.set_label(r'Pot. energy density difference [$10^7\:\mathrm{J m^{-2}}$]')
cbar.update_ticks()

#%% Save

fig.subplots_adjust(top=.97,bottom=.07,right=.98)
#fig.savefig('Epot2410.pdf')

