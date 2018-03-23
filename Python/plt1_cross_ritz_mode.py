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
import datetime as datetime
import matplotlib.ticker as ticker

plt.rcParams.update({'font.size':12})


#%% Figure

#Obslist
obs=plotter.read_obslist()
mat=sio.loadmat('V:/ipasmans/ritz_mode.mat',squeeze_me=True,struct_as_record=False)
mat=mat['vec']

#%% obs

tLim1=np.array([2395,2398])+datetime.datetime(2005,1,1).toordinal()
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
obsll=np.array(obsll)

#%% Map

plt.close('all')
fig=plt.figure(figsize=(4,4))
plt.rcParams.update({'font.size': 12})

#Colorbar
cmap=plt.cm.get_cmap('seismic')


#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,2):
        ax1=fig.add_subplot(2,1,i+1)
        #plotter.layout_surface(ax1,xLim=[-127,-123.5],yLim=[43,47])
        ax1.set_ylabel(r'Depth [m]')
        if i==1:
            ax1.set_xlabel(r'Longitude')
        ax1.set_ylim([300.,0.])
        ax1.set_xlim([-126,-124])

        ax1.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
        if i==1:
            ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        else:
            ax1.xaxis.set_major_formatter(ticker.NullFormatter())
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=50))
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
      
        ax.append(ax1)
    return ax
ax=plotAx()


#%% Create

fig.subplots_adjust(left=.2,bottom=.13)

#Interpolate to z-grid
iexp=0
z1=np.column_stack((mat.zr[:,0],mat.zr,mat.zr[:,-1]))
salt2=np.column_stack((mat.salt[:,0,iexp],mat.salt[:,:,iexp],mat.salt[:,-1,iexp]))
temp2=np.column_stack((mat.temp[:,0,iexp],mat.temp[:,:,iexp],mat.temp[:,-1,iexp]))
z2,salt2=plotter.z_interp(z1,salt2,np.array([-300.,0]))
z2,temp2=plotter.z_interp(z1,temp2,np.array([-300.,0]))
salt2=salt2
temp2=temp2

#Average depth
zmode=np.sum(mat.zr*np.squeeze(np.abs(mat.temp[:,:,iexp])),axis=1)/np.sum(np.squeeze(np.abs(mat.temp[:,:,iexp])),axis=1)

#Plot contourf
levelsf=[tick1 for tick1 in np.arange(-.08,.081,.01) if np.abs(tick1)>1e-8]
lon2=np.reshape(mat.lon,(-1,1))+np.zeros(np.shape(z2))
cplot1=ax[0].contourf(lon2,-z2,temp2,levels=levelsf,cmap=cmap,extend='both')
ax[0].plot(mat.lon[~np.isnan(zmode)],-zmode[~np.isnan(zmode)],color='k',linewidth=2)

#Glider extent
ax[0].plot(np.min(obsll[:,0])*np.array([1,1]),np.array([0,300]),color=(.5,.5,.5),linewidth=1.5,linestyle='--')
ax[0].plot(np.max(obsll[:,0])*np.array([1,1]),np.array([0,300]),color=(.5,.5,.5),linewidth=1.5,linestyle='--')

#Coast   
h=np.concatenate(([10e3],mat.h,[10e3]))
hlon=np.concatenate(([-135.],mat.lon,[-120.]))
p=patch.Polygon(np.column_stack((hlon,h)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
ax[0].add_patch(p)

ax[0].set_title('Surface Only')

#%% Combined

#Interpolate to z-grid
iexp=2
z1=np.column_stack((mat.zr[:,0],mat.zr,mat.zr[:,-1]))
salt2=np.column_stack((mat.salt[:,0,iexp],mat.salt[:,:,iexp],mat.salt[:,-1,iexp]))
temp2=np.column_stack((mat.temp[:,0,iexp],mat.temp[:,:,iexp],mat.temp[:,-1,iexp]))
z2,salt2=plotter.z_interp(z1,salt2,np.array([-300.,0]))
z2,temp2=plotter.z_interp(z1,temp2,np.array([-300.,0]))
salt2=salt2
temp2=temp2

#Average depth
zmode=np.sum(mat.zr*np.squeeze(np.abs(mat.temp[:,:,iexp])),axis=1)/np.sum(np.squeeze(np.abs(mat.temp[:,:,iexp])),axis=1)

#Plot contourf
levelsf=[tick1 for tick1 in np.arange(-.4,.41,0.05) if np.abs(tick1)>1e-8]
lon2=np.reshape(mat.lon,(-1,1))+np.zeros(np.shape(z2))
cplot1=ax[1].contourf(lon2,-z2,temp2,levels=levelsf,cmap=cmap,extend='both')
ax[1].plot(mat.lon[~np.isnan(zmode)],-zmode[~np.isnan(zmode)],color='k',linewidth=2)

#Glider extent
ax[1].plot(np.min(obsll[:,0])*np.array([1,1]),np.array([0,300]),color=(.5,.5,.5),linewidth=1.5,linestyle='--')
ax[1].plot(np.max(obsll[:,0])*np.array([1,1]),np.array([0,300]),color=(.5,.5,.5),linewidth=1.5,linestyle='--')

#Coast   
h=np.concatenate(([10e3],mat.h,[10e3]))
hlon=np.concatenate(([-135.],mat.lon,[-120.]))
p=patch.Polygon(np.column_stack((hlon,h)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
ax[1].add_patch(p)

ax[1].set_title('Combined')

#%%

plt.savefig('cross_ritz.png',dpi=600)









