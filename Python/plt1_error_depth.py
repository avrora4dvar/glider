# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:39:11 2018

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



#%% Read data

mat=[]
mat.append(sio.loadmat('V:/ipasmans/ritz_exp35.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:/ipasmans/ritz_exp37.mat',squeeze_me=True,struct_as_record=False))
grd=mat[0]['grd']
t=2395+datetime(2005,1,1).toordinal()

#Read grid
grd=plotter.read_grid_rho()

#Obslist
obs=plotter.read_obslist()

#%% Layout

plt.close('all')
fig=plt.figure(figsize=(5,5.5))

plt.rcParams.update({'font.size': 12})


#Colorbar
cmap=plt.cm.get_cmap('seismic')
cgrey=plt.cm.get_cmap('Greys')



#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,1):
        ax1=fig.add_subplot(1,1,1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        plotter.layout_surface(ax1,xLim=[-129,-123],yLim=[41,49])
        if i==0:
            ax1.set_ylabel(r'Latitude')
        ax1.set_xlabel(r'Longitude')
        if i>0:
            ax1.yaxis.set_major_locator(plt.NullLocator())
        ax1.set_xticks([-128,-124])

        ax.append(ax1)
        
    return ax
    
#Fill
ax=plotAx()
ax=ax[0]


z0=mat[0]['mode_depth'].temp
z1=mat[1]['mode_depth'].temp
d=(z1-z0)/np.abs(z0)
z0=mat[0]['mode_depth'].salt
z1=mat[1]['mode_depth'].salt
d=np.where( np.abs((z1-z0)/np.abs(z0))>np.abs(d),(z1-z0)/np.abs(z0),d)

cplot1=ax.contourf(grd.lon,grd.lat,d*100,[tick1 for tick1 in np.arange(-100.0,100.01,5.) if np.abs(tick1)>1],cmap=cmap,extend='both')
        
           
tLim1=3.*int((t-2)/3)+np.array([2.,5.])
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
lon1,lat1=zip(*obsll)
gplot1=ax.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        
#Colorbar
fig.subplots_adjust(right=.8)
cax=fig.add_axes([.75,.1,.02,.8])
cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional')
cbar.set_ticks([tick1 for tick1 in np.arange(-100.0,100.01,20.)])
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x))
cbar.update_ticks()
cbar.set_label(r'Rel. depth change [$\%$]')

plt.savefig('mode_depth.png',dpi=600)
#plt.tight_layout()
#%%
