# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 12:50:46 2018

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
from matplotlib import dates

imp.reload(plotter)

#%% Read data

#SST data
matT0=sio.loadmat('V:\ipasmans\\sst_rms_exp40_20110721_20110811.mat',squeeze_me=True,struct_as_record=False)

#HFR data
matUV0=sio.loadmat('V:\ipasmans\\hfr_rms_exp40_20110721_20110811.mat',squeeze_me=True,struct_as_record=False)

#Obslist
obs=plotter.read_obslist()

#%% Create figure

plt.close('all')
fig=plt.figure(figsize=(7,8.24*.4))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})


#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,3):
        ax1=fig.add_subplot(1,3,i+1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        #plotter.layout_surface(ax1,xLim=[-127,-123.5],yLim=[43,47])
        plotter.layout_surface(ax1,xLim=[-129,-123],yLim=[41,49])
        
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
        
        if np.mod(i,3)==0:
            ax1.set_ylabel(r'Latitude')
        else:
            ax1.yaxis.set_major_formatter(ticker.NullFormatter())
        if i>=0:
            ax1.set_xlabel(r'Longitude')
        else:
            ax1.xaxis.set_major_formatter(ticker.NullFormatter())
        #t=ax1.text(-128.6,48.5,mat[i]['name'])
        #t.set_bbox(dict(facecolor='w', alpha=1, edgecolor='w'))

        ax.append(ax1)
    return ax

#%% Calculate bias
    
biasT=matT0['error'].M1/matT0['error'].M0
biasU=matUV0['errorU'].M1/matUV0['errorU'].M0
biasV=matUV0['errorV'].M1/matUV0['errorV'].M0

#%% Plot biases

ax=plotAx()
cmap=plt.cm.get_cmap('seismic')
tLim1=dates.date2num([datetime(2011,7,21),datetime(2011,8,13)])

lon1=np.reshape(matT0['grd'].lon,(-1)); lat1=np.reshape(matT0['grd'].lat,(-1))
val1=np.reshape(biasT,(-1)); in1=~np.isnan(val1)
cmap0=ax[0].scatter(lon1[in1],lat1[in1],c=val1[in1],s=.1,cmap=cmap,vmin=-1.6,vmax=1.6)
plt.colorbar(cmap0,ax=ax[0],orientation='horizontal',extend='both',ticks=np.arange(-1.5,1.501,.5),label=u'Temperature bias [\u00B0C]')



p=matplotlib.patches.Rectangle(xy=(0.01,.92),width=.1,height=.07,transform=ax[0].transAxes,zorder=100,fill=True,color='w')
ax[0].add_patch(p)
ax[0].text(0.02,.93,'a)',transform=ax[0].transAxes,color='k',zorder=150)

lon1=np.reshape(matUV0['grd'].lon,(-1)); lat1=np.reshape(matUV0['grd'].lat,(-1))
val1=np.reshape(biasU,(-1)); in1=~np.isnan(val1)
cmap1=ax[1].scatter(lon1[in1],lat1[in1],c=val1[in1],s=.6,cmap=cmap,vmin=-.3,vmax=.3)
plt.colorbar(cmap1,ax=ax[1],orientation='horizontal',extend='both',ticks=np.arange(-.3,0.301,.1),label=r'Zonal velocity bias [ms$^{-1}$]')

p=matplotlib.patches.Rectangle(xy=(0.01,.92),width=.1,height=.07,transform=ax[1].transAxes,zorder=100,fill=True,color='w')
ax[1].add_patch(p)
ax[1].text(0.02,.93,'b)',transform=ax[1].transAxes,color='k',zorder=150)


lon1=np.reshape(matUV0['grd'].lon,(-1)); lat1=np.reshape(matUV0['grd'].lat,(-1))
val1=np.reshape(biasV,(-1)); in1=~np.isnan(val1)
cmap2=ax[2].scatter(lon1[in1],lat1[in1],c=val1[in1],s=.6,cmap=cmap,vmin=-.3,vmax=.3)
plt.colorbar(cmap2,ax=ax[2],orientation='horizontal',extend='both',ticks=np.arange(-.3,0.31,.1),label=r'Meridional velocity bias [ms$^{-1}$]')

p=matplotlib.patches.Rectangle(xy=(0.01,.92),width=.1,height=.07,transform=ax[2].transAxes,zorder=100,fill=True,color='w')
ax[2].add_patch(p)
ax[2].text(0.02,.93,'c)',transform=ax[2].transAxes,color='k',zorder=150)

    
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
lon1,lat1=zip(*obsll)
gplot1=ax[0].plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=1,fillstyle='full')
gplot1=ax[1].plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=1,fillstyle='full')
gplot1=ax[2].plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=1,fillstyle='full')

#%%

fig.subplots_adjust(left=.01,right=.99,top=.98,bottom=.00)
fig.savefig('bias_sst_hfr.pdf')