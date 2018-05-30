# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:58:27 2018

@author: ipasmans

Contains plotting routines for data

"""

import matplotlib.pyplot as plt
import matplotlib.patches as patch
import numpy as np
import scipy.io as sio
import netCDF4
import os 
from datetime import datetime, timedelta
import re

def num2time(tIn):
    tOut=[]
    try:
        for t1 in tIn:
            tOut.append(datetime.fromordinal(int(t1))+timedelta(days=t1%1))
    except TypeError:
        t1=tIn
        tOut=datetime.fromordinal(int(t1))+timedelta(days=t1%1)
    return tOut

def time2num(tIn):
    tOut=[]
    try:
        for t1 in tIn:
            tOut.append( datetime.toordinal(t1)+t1.hour/24.+t1.minute/(24.*60.)+t1.second/(24.*3600.))
    except TypeError:
        t1=tIn
        tOut=datetime.toordinal(t1)+t1.hour/24.+t1.minute/(24.*60.)+t1.second/(24.*3600.)
    return tOut
        

def read_grid_rho():
    """
    Read rho-grid from netcdf file
    """
    nc=netCDF4.Dataset('V:/ipasmans/grd_ow2km_r13lin_mix.nc','r',format='NETCDF3_classic')
    grd={'lon':np.transpose(np.array(nc['lon_rho'][:,:])),
    'lat':np.transpose(np.array(nc['lat_rho'][:,:])),
    'mask':np.transpose(np.array(nc['mask_rho'])),
    'h':np.transpose(np.array(nc['h'][:,:])),
    'zr':np.array(nc['z0_r'][:,:,:]).transpose((2,1,0)),
    'zw':np.array(nc['z0_w'][:,:,:]).transpose((2,1,0))}
    nc.close()
    return grd

def read_obslist():
    """
    Read obslist files
    """
    obsDir='V:/ipasmans/Exp37/'
    indir=os.listdir(obsDir)
    obs={'t':np.array([]),'lon':np.array([]),'lat':np.array([]),'type':np.array([]),
    'val':np.array([]),'z':np.array([]),'sig':np.array([]),'winNo':np.array([])}
    for fname in indir:
        if not fname.startswith('obslist'):
            continue
        nc=netCDF4.Dataset(obsDir+fname,'r',format='NETCDF3_classic')
        obsno=re.findall('obslist(\d+).nc',fname)
        obsno=float(obsno[0])+731947
        obs['lon']=np.append(obs['lon'],nc['lon'][:])
        obs['lat']=np.append(obs['lat'],nc['lat'][:])
        obs['t']=np.append(obs['t'],nc['time'][:]/3600/24+obsno)
        obs['type']=np.append(obs['type'],nc['type'][:])
        obs['z']=np.append(obs['z'],nc['z'][:])
        obs['val']=np.append(obs['val'],nc['obs'][:])
        obs['sig']=np.append(obs['sig'],nc['sig_d'][:])
    return obs

def add_coast(ax):
    """
    Add coastline
    """
    coast=sio.loadmat('V:/ipasmans/coastLine.mat',squeeze_me=True,struct_as_record=False)
    coast=coast['coastLine']
    for c1 in coast:
        p=patch.Polygon(np.transpose([c1.Lon,c1.Lat]),closed=True,facecolor=(.5,.6,.5),edgecolor=None)
        ax.add_patch(p)
    return ax
        
def add_bathy(ax,levels=[200,1000,2000]):
    """
    Add batymetry
    """
    levels=np.array(levels)
    grd=read_grid_rho()
    ax.contour(grd['lon'],grd['lat'],grd['h'],levels=levels,colors=(.7,.7,.7),linewidths=.6)
    return ax
    
def add_glider2d(ax,tLim=[-np.Inf,np.Inf]):
    """
    Add glider locations on surface
    """
    obs=read_obslist()
    ll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and np.min(tLim)<=t1<np.max(tLim) ]
    lon,lat=zip(*ll)
    p=ax.plot(lon,lat,'k,')
    return p

def layout_surface(ax,xLim=[-130,-123.5],yLim=[41,50]):
    ax.set_xlim(xLim)
    ax.set_ylim(yLim)
    ax.set_aspect(1./(.5*np.sqrt(2)))
    return ax
    
def add_uv2d(ax,lon,lat,u,v):
    nx=np.diff(ax.get_xlim())/(35*np.mean(np.diff(lon[:,0])))
    ny=np.diff(ax.get_ylim())/(35*np.mean(np.diff(lat[0,:])))
    nx=int(np.round(nx)); ny=int(np.round(ny))
    
    lon=lon[::nx,::ny]; lat=lat[::nx,::ny]
    u=u[::nx,::ny]; v=v[::nx,::ny]
    
    q1=ax.quiver(lon,lat,u,v,scale_units='width',scale=1.*float(nx),angles='uv',color=(0,0,0),width=0.006,headwidth=6,headlength=9)
    q2=ax.quiver(-124.2,42.8,1,0,scale_units='width',scale=1.*float(nx),angles='uv',color=(0,0,0),width=0.006,headwidth=6,headlength=9,zorder=100)
    q=(q1,q2)
    
    return q

def z_interp(zIn,valIn,zLim):
    """
    Interpolate in z-direction to get smooth contourf plot
    """
    zOut=np.arange(zLim[0],zLim[1],4.*np.sign(zLim[1]-zLim[0]))
    zOut=np.reshape(zOut,(1,-1))
    zOut=zOut+np.zeros((np.size(zIn,0),np.size(zOut,1)))
    zOut[:,-1]=zIn[:,-1]
    
    valOut=np.zeros(np.shape(zOut))
    
    for i0 in np.arange(0,np.size(valIn,0)):
        valOut[i0,:]=np.interp(zOut[i0,:],xp=zIn[i0,:],fp=valIn[i0,:])
        
    return zOut,valOut

