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

mat=(sio.loadmat('V:/ipasmans/rho_contours_ana.mat',squeeze_me=True,struct_as_record=False))
mat=mat['model']

#Adjust times to Python
for mat1 in mat:
    mat1.t=mat1.t-366.
dateRef=datetime(2005,1,1).toordinal()

#Read grid
grd=plotter.read_grid_rho()

#Obslist
obs=plotter.read_obslist()

plt.close('all')
fig=plt.figure(figsize=(3.74,4.53))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})

#%% 

mat=sio.loadmat('V:/ipasmans/Epot.mat',squeeze_me=True,struct_as_record=False)
mat['t']=mat['t']-366.

#%% Layout

#Colorbar
cmap=plt.cm.get_cmap('seismic')

#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,1):
        ax1=fig.add_subplot(1,1,i+1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        plotter.layout_surface(ax1,xLim=[-129,-123],yLim=[41,49])

        ax1=plotter.add_bathy(ax1)
        ax1=plotter.add_coast(ax1)
        
        ax1.set_ylabel(r'Latitude')
        ax1.set_xlabel(r'Longitude')
     
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
        ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))

        ax.append(ax1)
        
    return ax
ax=plotAx()

#%% Plot

#
ax1=ax[0]
t=plotter.time2num(datetime(2011,7,21,12,0,0))
it=np.abs(mat['t']-t)<.1
cplot1=ax1.contourf(mat['lon'],mat['lat'],np.squeeze(mat['zpot'][1][:,:,it]-mat['zpot'][0][:,:,it]),levels=[tick1 for tick1 in np.arange(-.1,.101,.01) if np.abs(tick1)>.001],cmap=cmap)
cplot2=ax1.contour(mat['lon'],mat['lat'],np.squeeze(mat['dEpot'][1][:,:,it]),colors='k',linewidths=.7,levels=np.arange(0,2.01e5,.1e5))

#Glider location
tLim1=3.*int((t-2)/3)+np.array([2.,5.])
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
lon1,lat1=zip(*obsll)
gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        
#Colorbar
cax=fig.add_axes([.79,.11,.03,.77])
cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional')
cbar.set_ticks([tick1 for tick1 in np.arange(-.1,.101,.02)])
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.2f}'.format(x))
cbar.update_ticks()
cbar.set_label(r'Depth [$\mathrm{m}$]')
        

