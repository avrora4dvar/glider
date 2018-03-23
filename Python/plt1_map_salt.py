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
import datetime

plt.rcParams.update({'font.size':12})

dateRef=datetime.datetime(2005,1,1,0,0,0).toordinal()
tLim=np.array([2392,2395])+dateRef
tLim=np.array([2412.1,2412.9])+dateRef

#%% Figure

mat=[]
mat.append(sio.loadmat('V:/ipasmans/Exp36for_avg.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:/ipasmans/Exp37for_avg.mat',squeeze_me=True,struct_as_record=False))

for mat1 in mat:
    mat1['t']=mat1['t']-366

grd=plotter.read_grid_rho()

#Obslist
obs=plotter.read_obslist()

#%% Map

plt.close('all')
fig=plt.figure(figsize=(3,4.5))
plt.rcParams.update({'font.size': 12})


#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,1):
        ax1=fig.add_subplot(1,1,i+1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        #plotter.layout_surface(ax1,xLim=[-127,-123.5],yLim=[43,47])
        plotter.layout_surface(ax1,xLim=[-130,-122.2],yLim=[40.65,50])
        if i==0:
            ax1.set_ylabel(r'Latitude')
        ax1.set_xlabel(r'Longitude')
        if i>0:
            ax1.yaxis.set_major_locator(plt.NullLocator())
        ax.append(ax1)
    return ax
ax=plotAx()
ax=ax[0]


#%% Add SST

#Colorbar
cmap=plt.cm.get_cmap('plasma')

#Plot
ax.set_position([.15,.15,.6,.8])
levelsf=np.arange(-1,1.01,.1)
ind=np.all((obs['type']==4,obs['t']>=tLim[0],obs['t']<tLim[1]),axis=0)
cplot1=ax.scatter(obs['lon'][ind],obs['lat'][ind],c=obs['val'][ind],s=.2,marker='o',cmap=cmap,
                  vmin=10,vmax=16)

it=mat[0]['t']==2412.5+dateRef
ax.contour(mat[0]['lon'],mat[0]['lat'],np.squeeze(mat[0]['salt'][:,:,it]),levels=[31.5],linestyles=':',colors='k',linewidths=1.5)
ax.contour(mat[1]['lon'],mat[1]['lat'],np.squeeze(mat[1]['salt'][:,:,it]),levels=[31.5],linestyles='-',colors='k',linewidths=1.5)

#Add NH10
#pplot=ax.scatter(-124.3,44.642,s=10,c='k',marker='o')
#ax.text(-124.3-.4,44.642+.1,'NH10')


#Colorbar
cax=fig.add_axes([.75,.22,.02,.56])
cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional',extend='both')
cbar.locator=ticker.MultipleLocator(base=1)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.0f}'.format(x))
cbar.update_ticks()
cbar.set_label(r'Sea-surface temperature [$\mathrm{^\circ C}$]')


fig.subplots_adjust(left=.2,right=.7,bottom=.1)
plt.savefig('map_SST_plume.png',dpi=600)

