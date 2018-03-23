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
import matplotlib
imp.reload(plotter)


#%% Read data

mat35=(sio.loadmat('V:/ipasmans/ritz_vectors_exp35.mat',squeeze_me=True,struct_as_record=False))
mat36=(sio.loadmat('V:/ipasmans/ritz_vectors_exp36.mat',squeeze_me=True,struct_as_record=False))
mat37=(sio.loadmat('V:/ipasmans/ritz_vectors_exp37.mat',squeeze_me=True,struct_as_record=False))
mat=mat35
t=plotter.time2num(datetime(2011,7,24))

plt.close('all')
fig=plt.figure(figsize=(7.4,8.2))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})

#Read grid
grd=plotter.read_grid_rho()

#Obslist
obs=plotter.read_obslist()


#%% Create axes

ax=[]
for iax in np.arange(0,24):
    ax1=fig.add_subplot(6,4,iax+1)
    
    if iax<12:
        ax1=plotter.add_coast(ax1)
        ax1=plotter.add_bathy(ax1)
        ax1.set_xlim(-129,-123.5)
        ax1.set_ylim(41,49)
        
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
        
        if np.mod(iax,4)==0:
            ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        else:
            ax1.yaxis.set_major_formatter(ticker.NullFormatter())
        if iax>=20:
            ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        else:
            ax1.xaxis.set_major_formatter(ticker.NullFormatter())
        
        if iax<4:
            ax1.set_title('Mode {:d}'.format(iax+1))
    else:
        ax1.set_xlim(-129,-123.5)
        ax1.set_ylim(300,0)
        
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
        
        if np.mod(iax,4)==0:
            ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        else:
            ax1.yaxis.set_major_formatter(ticker.NullFormatter())
        if iax>=20:
            ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
            ax1.set_xlabel('Longitude')
        else:
            ax1.xaxis.set_major_formatter(ticker.NullFormatter())
        
    ax.append(ax1)
    
#%% Plot exp35

mat=mat35
cm=plt.cm.get_cmap('seismic')
levels=[tick1 for tick1 in np.arange(-1,1.01,.1) if np.abs(tick1)>.01]

ax[0].set_ylabel(r'Surface Only')
maxVal=[]
for iax in np.arange(0,4):
    val1=mat['surface'].temp[:,:,iax-0]
    maxVal1=np.nanmax(np.abs(val1))
    val1=val1/maxVal1; maxVal.append(maxVal1)
    ax[iax].contourf(mat['grd'].lon,mat['grd'].lat,val1,levels=levels,cmap=cm,vmin=-1,vmax=1)
    ax[iax].plot(np.array([-130,120]),mat['cross'].lat[iax]*np.array([1,1]),'k-',linewidth=1)
    
    ax[iax].text(-128.5,48,'{:.0f}%'.format(np.round(100*mat['p_ritzS'][iax])))
    
    #Glider location
    tLim1=3.*int((t-2)/3)+np.array([2.,5.])
    obsll=[(lon1,lat1,z1) for (lon1,lat1,type1,t1,z1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t'],obs['z']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
    lonObs,latObs,zObs=zip(*obsll)
    gplot1=ax[iax].plot(lonObs,latObs,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        

ax[12].set_ylabel(r'Surface Only')
s=np.array([1,1,1,1])
for iax in np.arange(12,16):
    val1=mat['cross'].temp[:,:,iax-12]
    maxVal1=np.nanmax(np.abs(val1))
    val1=s[iax-12]*val1/maxVal1
    lon1=np.reshape(mat['grd'].lon[:,0],(-1,1))+mat['cross'].zr*0
    ax[iax].contourf(lon1,-mat['cross'].zr,val1,levels=levels,cmap=cm,vmin=-1,vmax=1)
    
    ax[iax].text(-128.5,48,'{:.0f}%'.format(np.round(100*mat['p_ritzS'][iax-12])))
    
    lon1=mat['grd'].lon[:,0]
    lon1=np.concatenate([[-140],lon1,[-120]])
    h1=mat['grd'].h[:,mat['cross'].ilat[iax-12]-1]
    h1=np.concatenate([[10e3],h1,[10e3]])
    p=patch.Polygon(np.column_stack((lon1,h1)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
    ax[iax].add_patch(p)
    
    gplot1=ax[iax].plot([np.min(lonObs),np.min(lonObs),np.max(lonObs),np.max(lonObs)],
    [np.min(zObs),np.max(zObs),np.max(zObs),np.min(zObs)],':',color=(.5,.5,.5),linewidth=1.)
    
#%% Plot exp36

mat=mat36
cm=plt.cm.get_cmap('seismic')
levels=[tick1 for tick1 in np.arange(-1,1.01,.1) if np.abs(tick1)>.01]

ax[4].set_ylabel(r'Glider Only')
maxVal=[]
for iax in np.arange(4,8):
    val1=mat['surface'].temp[:,:,iax-4]
    maxVal1=np.nanmax(np.abs(val1)); maxVal.append(maxVal1)
    val1=val1/maxVal1
    ax[iax].contourf(mat['grd'].lon,mat['grd'].lat,val1,levels=levels,cmap=cm,vmin=-1,vmax=1)
    ax[iax].plot(np.array([-130,120]),mat['cross'].lat[iax-4]*np.array([1,1]),'k-',linewidth=1)
    
    ax[iax].text(-128.5,48,'{:.0f}%'.format(np.round(100*mat['p_ritzS'][iax-4])))
    
    #Glider location
    tLim1=3.*int((t-2)/3)+np.array([2.,5.])
    obsll=[(lon1,lat1,z1) for (lon1,lat1,type1,t1,z1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t'],obs['z']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
    lonObs,latObs,zObs=zip(*obsll)
    gplot1=ax[iax].plot(lonObs,latObs,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        

ax[16].set_ylabel(r'Glider Only')
for iax in np.arange(16,20):
    val1=mat['cross'].temp[:,:,iax-16]
    maxVal1=np.nanmax(np.abs(val1))
    val1=val1/maxVal1
    lon1=np.reshape(mat['grd'].lon[:,0],(-1,1))+mat['cross'].zr*0
    ax[iax].contourf(lon1,-mat['cross'].zr,val1,levels=levels,cmap=cm,vmin=-1,vmax=1)
    
    ax[iax].text(-128.5,48,'{:.0f}%'.format(np.round(100*mat['p_ritzS'][iax-16])))
    
    lon1=mat['grd'].lon[:,0]
    lon1=np.concatenate([[-140],lon1,[-120]])
    h1=mat['grd'].h[:,mat['cross'].ilat[iax-16]-1]
    h1=np.concatenate([[10e3],h1,[10e3]])
    p=patch.Polygon(np.column_stack((lon1,h1)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
    ax[iax].add_patch(p)
    
    gplot1=ax[iax].plot([np.min(lonObs),np.min(lonObs),np.max(lonObs),np.max(lonObs)],
    [np.min(zObs),np.max(zObs),np.max(zObs),np.min(zObs)],':',color=(.5,.5,.5),linewidth=1.)
    
#%% Plot exp37

mat=mat37
cm=plt.cm.get_cmap('seismic')
levels=[tick1 for tick1 in np.arange(-1,1.01,.1) if np.abs(tick1)>.01]

ax[8].set_ylabel(r'Combined')
maxVal=[]
s=np.array([1,1,1,1])
for iax in np.arange(8,12):
    val1=mat['surface'].temp[:,:,iax-8]
    maxVal1=np.nanmax(np.abs(val1)); maxVal.append(maxVal1)
    val1=val1/maxVal1
    ax[iax].contourf(mat['grd'].lon,mat['grd'].lat,val1,levels=levels,cmap=cm,vmin=-1,vmax=1)
    ax[iax].plot(np.array([-130,120]),mat['cross'].lat[iax-8]*np.array([1,1]),'k-',linewidth=1)
    
    ax[iax].text(-128.5,48,'{:.0f}%'.format(np.round(100*mat['p_ritzS'][iax-8])))
    
    #Glider location
    tLim1=3.*int((t-2)/3)+np.array([2.,5.])
    obsll=[(lon1,lat1,z1) for (lon1,lat1,type1,t1,z1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t'],obs['z']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
    lonObs,latObs,zObs=zip(*obsll)
    gplot1=ax[iax].plot(lonObs,latObs,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        

ax[20].set_ylabel(r'Combined')
for iax in np.arange(20,24):
    val1=mat['cross'].temp[:,:,iax-20]
    maxVal1=np.nanmax(np.abs(val1))
    val1=val1/maxVal1
    lon1=np.reshape(mat['grd'].lon[:,0],(-1,1))+mat['cross'].zr*0
    cplot1=ax[iax].contourf(lon1,-mat['cross'].zr,val1,levels=levels,cmap=cm,vmin=-1,vmax=1)
    
    ax[iax].text(-128.5,48,'{:.0f}%'.format(np.round(100*mat['p_ritzS'][iax-20])))
    
    lon1=mat['grd'].lon[:,0]
    lon1=np.concatenate([[-140],lon1,[-120]])
    h1=mat['grd'].h[:,mat['cross'].ilat[iax-20]-1]
    h1=np.concatenate([[10e3],h1,[10e3]])
    p=patch.Polygon(np.column_stack((lon1,h1)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
    ax[iax].add_patch(p)
    
    gplot1=ax[iax].plot([np.min(lonObs),np.min(lonObs),np.max(lonObs),np.max(lonObs)],
    [np.min(zObs),np.max(zObs),np.max(zObs),np.min(zObs)],':',color=(.5,.5,.5),linewidth=1.)
    
#%% Save

#Colorbar
cax=fig.add_axes([.12,.04,.78,.015])
cbar=fig.colorbar(cplot1,cax=cax,orientation='horizontal')
cbar.set_ticks([tick1 for tick1 in np.arange(-1,1.01,.2)])
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x))
cbar.update_ticks()
#cbar.set_label(r'Density-$10^3$ [$\mathrm{kg m^{-3}}$]')
    
fig.subplots_adjust(top=.97)
fig.savefig('array_modes.pdf')