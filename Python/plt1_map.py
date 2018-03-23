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
import matplotlib.ticker as ticker

plt.rcParams.update({'font.size':8})
matplotlib.rcParams['font.family']='Arial'
plt.close('all')
fig=plt.figure(figsize=(3.66,8.92*.35))
#fig=plt.figure(figsize=(5.61,8.92*.35))

tLim=np.array([2392,2395])+datetime.datetime(2005,1,1,0,0,0).toordinal()
#tLim=np.array([2412,2413])+datetime.datetime(2005,1,1,0,0,0).toordinal()

#%% Figure

#Obslist
obs=plotter.read_obslist()
grd=plotter.read_grid_rho()

#%% Map

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
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
        if i>0:
            ax1.yaxis.set_major_locator(plt.NullLocator())
        #ax1.spines['top'].set_linewidth(.5)
        #ax1.spines['bottom'].set_linewidth(.5)
        #ax1.spines['left'].set_linewidth(.5)
        #ax1.spines['right'].set_linewidth(.5)
        ax.append(ax1)
    return ax

ax=plotAx()


#%% Add velocity

ax1=ax[0]

#Colorbar
cmap=plt.cm.get_cmap('seismic')

#Plot
#ax1.set_position([.15,.15,.6,.8])
levelsf=np.arange(-1,1.01,.1)
ind=np.all((obs['type']==8,obs['t']>=tLim[0],obs['t']<tLim[1]),axis=0)
cplot1=ax1.scatter(obs['lon'][ind],obs['lat'][ind],c=obs['val'][ind],s=.3,marker='o',cmap=cmap,vmin=-1,vmax=1)

#Add glider
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim[0]<=t1<=tLim[1]]
lon1,lat1=zip(*obsll)
gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=1,fillstyle='full')
   

#Add NH10
#pplot=ax.scatter(-124.3,44.642,s=20,c='k',marker='o')
#ax.text(-124.3-.4,44.642+.1,'NH10')

#Colorbar
ax1.text(-129.8,49.4,r'a)')
#cax=fig.add_axes([.75,.15,.02,.8])
cbar=fig.colorbar(cplot1,ax=ax1,orientation='horizontal',spacing='proportional')
cbar.set_ticks(levelsf[[0,-1]])
cbar.locator=ticker.MultipleLocator(base=.5)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x))
cbar.update_ticks()
cbar.set_label(r'Radial velocity [$\mathrm{m s^{-1}}$]')

#fig.add_axes
#plt.axis('tight')
#plt.savefig('map_uv.png',dpi=600)

#%% Add SST

ax1=ax[1]

#Colorbar
cmap=plt.cm.get_cmap('jet')

#Plot
#ax1.set_position([.15,.15,.6,.8])
levelsf=np.arange(-1,1.01,.1)
ind=np.all((obs['type']==4,obs['t']>=tLim[0],obs['t']<tLim[1]),axis=0)
cplot1=ax1.scatter(obs['lon'][ind],obs['lat'][ind],c=obs['val'][ind],s=.05,marker='o',cmap=cmap,
                  vmin=np.floor(np.min(obs['val'][ind])),vmax=np.ceil(np.max(obs['val'][ind])))

#Add glider
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim[0]<=t1<=tLim[1]]
lon1,lat1=zip(*obsll)
gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=1,fillstyle='full')

#Add NH10
#pplot=ax.scatter(-124.3,44.642,s=10,c='k',marker='o')
#ax.text(-124.3-.4,44.642+.1,'NH10')

#Colorbar
ax1.text(-129.8,49.4,r'b)')
#cax=fig.add_axes([.75,.15,.02,.8])
cbar=fig.colorbar(cplot1,ax=ax1,orientation='horizontal',spacing='proportional')
cbar.set_ticks(levelsf[[0,-1]])
cbar.locator=ticker.MultipleLocator(base=2)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.0f}'.format(x))
cbar.update_ticks()
cbar.set_label(u'Sea-surface temperature [\u00B0C]')

#fig.add_axes
#plt.axis('tight')
#plt.savefig('map_SST.png',dpi=600)


#%% Add SSH

ax1=ax[2]

#Colorbar
cmap=plt.cm.get_cmap('seismic')

#Plot
#ax1.set_position([.15,.15,.6,.8])
levelsf=np.arange(-1,1.01,.1)
ind=np.all((obs['type']>=1000,obs['t']>=tLim[0],obs['t']<tLim[1]),axis=0)
cplot1=ax1.scatter(obs['lon'][ind],obs['lat'][ind],c=100*obs['val'][ind],s=.15,marker='o',vmin=-16,vmax=16,cmap=cmap)

#Add NH10
#pplot=ax.scatter(-124.3,44.642,s=10,c='k')
#ax.text(-124.3-.4,44.642+.1,'NH10')

#Add glider
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim[0]<=t1<=tLim[1]]
lon1,lat1=zip(*obsll)
gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=1,fillstyle='full')
   

#Colorbar
ax1.text(-129.8,49.4,r'c)')
#cax=fig.add_axes([.75,.15,.02,.8])
cbar=fig.colorbar(cplot1,ax=ax1,orientation='horizontal',spacing='proportional')
cbar.set_ticks(levelsf[[0,-1]])
cbar.locator=ticker.MultipleLocator(base=8)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.0f}'.format(x))
cbar.update_ticks()
cbar.set_label(r'SSH - track mean [$\mathrm{cm}$]')

#fig.add_axes
#plt.axis('tight')
#plt.savefig('map_SSH.png',dpi=600)
fig.subplots_adjust(bottom=.02,top=.98)
plt.savefig('map_obs.pdf')


#%% Add velocity

#Colorbar
cmap=plt.cm.get_cmap('Greys')

mat=sio.loadmat('V:\ipasmans\Exp40ana_avg.mat',squeeze_me=True,struct_as_record=False)

#Plot
ax[0].set_position([.15,.15,.6,.8])
mask=np.where(np.all(np.stack((grd['h']<1e3,grd['lon']<-123.6,grd['lat']<48),axis=2),axis=2),1,0)
cplot1=ax[0].contourf(grd['lon'],grd['lat'],mask,alpha=.7,vmin=-.1,vmax=2,cmap=cmap)

ax[0].contour(mat['lon'],mat['lat'],mat['salt'][:,:,1],levels=[31.5],colors='k',linestyles='-',linewidths=.7)
ax[0].contour(mat['lon'],mat['lat'],mat['salt'][:,:,-1],levels=[31.5],colors='k',linestyles=':',linewidths=.7)

#Add glider
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6]
lon1,lat1=zip(*obsll)
gplot1=ax[0].plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=1,fillstyle='full')
 
ax[0].annotate('Columbia River',xy=(-123.8,46.2),xytext=(-126.1,45.2),xycoords='data',textcoords='data',arrowprops=dict(facecolor='black', shrink=0.05,headwidth=1,width=.5,headlength=5))

fig.subplots_adjust(bottom=.12,top=.98)
plt.savefig('map_surface_area.pdf')
 
 


