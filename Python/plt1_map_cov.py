# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 11:43:24 2017

Create map of study area

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
import seawater as sw

plt.rcParams.update({'font.size':12})


#%% Figure

#Obslist
obs=plotter.read_obslist()

mat=sio.loadmat('V:\ipasmans\ocean2018\cov_bal.mat',squeeze_me=True,struct_as_record=False)
cov=mat['cov']; grd=mat['grd']; loc=mat['loc']

#%% Map

plt.close('all')
fig=plt.figure(figsize=(3.5,4.5),dpi=300)
plt.rcParams.update({'font.size': 12})


#Colorbar
cmap=plt.cm.get_cmap('seismic')


#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,1):
        ax1=fig.add_subplot(1,1,i+1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        #plotter.layout_surface(ax1,xLim=[-127,-123.5],yLim=[43,47])
        plotter.layout_surface(ax1,xLim=[-125.5,-123.5],yLim=[43.8,45.6])
        ax1.set_ylabel(r'Latitude')
        ax1.set_xlabel(r'Longitude')

        ax1.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
        ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=.2))
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.1f}'))
      
        ax.append(ax1)
    return ax
ax=plotAx()
ax=ax[0]

#%% Create

fig.subplots_adjust(left=.14,right=.7,bottom=.15,top=.95)

iz=-1
levelsf=np.arange(-.16,.001,.01)
levels=np.arange(0.1,1.01,.1)
cmap=plt.cm.get_cmap('YlGnBu_r')

rho1=7.5e-4*cov.salt[:,:,iz]*1025-1.7e-4*1025*cov.temp[:,:,iz]


cplot1=ax.contourf(grd.lon,grd.lat,cov.salt[:,:,iz],levels=levelsf,cmap=cmap,extend='min')
cplot2=ax.contour(grd.lon,grd.lat,cov.temp[:,:,iz],levels=levels,linewidths=1,colors='k')


step=3
u1=.5*cov.u[0:-2,1:-2,iz]+.5*cov.u[1:-1,1:-2,iz]
v1=.5*cov.v[1:-2,0:-2,iz]+.5*cov.v[1:-2,1:-1,iz]
print(np.array([np.max(np.abs(u1)),np.max(np.abs(v1))])/.9)
qplot=ax.quiver(grd.lon[1:-2:step,1:-2:step],grd.lat[1:-2:step,1:-2:step],
          u1[::step,::step],v1[::step,::step],angles='uv',units='width',scale=1/10)

#Point
ax.plot(loc.lon,loc.lat,'ko',markerfacecolor='k',markersize=4)

#Colorbar
cax=fig.add_axes([.71,.17,.02,.75])
cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional')
cbar.set_ticks(levelsf[[0,-1]])
cbar.locator=ticker.MultipleLocator(base=.02)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.2f}'.format(x))
cbar.update_ticks()
cbar.set_label(r'T,S-covariance [$\mathrm{ppt ^\circ C}$]')


#%%

plt.savefig('cross_cov_map.pdf')

