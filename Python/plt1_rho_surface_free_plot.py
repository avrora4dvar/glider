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


tOut=[datetime(2011,7,21,12,0,0),datetime(2011,8,10,12,0,0)]
tOut=plotter.time2num(tOut)
tOut=np.arange(tOut[0],tOut[1],2./24.)



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

levelList=[1026.5]

mat0=mat[-1]
mat=mat[[1,0,2]]
titleList=['Glider Only','Surface Only','Combined']
labelList=['a)','b)','c)','d)','e)','f)']

plt.close('all')
fig=plt.figure(figsize=(5.61,8.92*.75))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})

#%% Layout


#Colorbar
cmap=plt.cm.get_cmap('seismic')
cgrey=plt.cm.get_cmap('Greys')



#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,6):
        ax1=fig.add_subplot(2,3,i+1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        plotter.layout_surface(ax1,xLim=[-129,-123],yLim=[41,49])
        
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(2))
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
        if np.mod(i,3)==0:
            ax1.set_ylabel(r'Latitude')
            ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        else:
            ax1.yaxis.set_major_formatter(ticker.NullFormatter())
        if i>=3:
            ax1.set_xlabel(r'Longitude')
            ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        else:
            ax1.xaxis.set_major_formatter(ticker.NullFormatter())
            ax1.set_title(titleList[i])
            
        t=ax1.annotate(labelList[i],xycoords='axes fraction',textcoords='axes fraction',xy=(.06,.94),xytext=(.06,.94))
        t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
        ax.append(ax1)
        
    return ax

fig.clf()
ax=plotAx()
    
#Fill
im=[]
def plotFrame(t,ax):
    cplot=[]; gplot=[]; qplot=[]; bcplot=[]
    level1=levelList[0]
    
    for mat1,ax1,title1 in zip(mat,ax,titleList):
        #Time
        #fig.suptitle(plotter.num2time(t).strftime('%Y-%b-%d %H:%M')
        
        inCon=level1==mat1.contours
        it=np.interp(t,mat1.t,np.arange(0,len(mat1.t)))
        it1=np.minimum(int(it)+1,len(mat1.t)-1)        
        
        z0=mat0.z[:,:,inCon,int(it)]*(1-it%1)+mat0.z[:,:,inCon,it1]*(it%1)
        z1=mat1.z[:,:,inCon,int(it)]*(1-it%1)+mat1.z[:,:,inCon,it1]*(it%1)
        z0=np.squeeze(z0); z1=np.squeeze(z1)
        mask0=np.where(np.isnan(z0),np.ones(np.shape(z0)),np.zeros(np.shape(z0)))
        mask1=np.where(np.isnan(z1-z0),np.ones(np.shape(z1)),np.zeros(np.shape(z1)))
        print([np.nanpercentile(z1-z0,2.5),np.nanpercentile(z1-z0,97.5)])

        cplot1=ax1.contourf(mat1.lon,mat1.lat,z1-z0,[tick1 for tick1 in np.arange(-60.0,60.01,5.0) if np.abs(tick1)>1],cmap=cmap,extend='both')
        cplot1=ax1.contourf(mat1.lon,mat1.lat,z1-z0,[tick1 for tick1 in np.arange(-60.0,60.01,5.0) if np.abs(tick1)>1],cmap=cmap,extend='both')
        
        if level1<=1024.5:
            bcplot1=ax1.contour(mat1.lon,mat1.lat,z1,levels=np.arange(-200,0.1,10.),colors='k',linewidths=.7)
        else:
            bcplot1=ax1.contour(mat1.lon,mat1.lat,z1,levels=np.arange(-2000,0.1,25.),colors='k',linewidths=.7)
        #ax1.set_title(title1)
        bcplot2=ax1.contourf(mat1.lon,mat1.lat,mask1,levels=[-.5,.5,1],cmap=cgrey,vmin=0.5,vmax=1,alpha=.5)

        
        tLim1=3.*int((t-2)/3)+np.array([2.,5.])
        obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
        lon1,lat1=zip(*obsll)
        gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        
        #Colorbar
        #fig.subplots_adjust(top=1.2)
        cax=fig.add_axes([.1,.07,.8,.02])
        cbar=fig.colorbar(cplot1,cax=cax,orientation='horizontal',spacing='proportional')
        cbar.set_ticks([tick1 for tick1 in np.arange(-60.0,60.1,20.)])
        cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x))
        cbar.update_ticks()
        cbar.set_label(r'Depth change [$\mathrm{m}$]')
        
        cplot.append(cplot1)
        gplot.append(gplot1)
        bcplot.append(bcplot1)
        
        #plt.tight_layout()
        
    return cplot,gplot,bcplot



#%%

t=datetime(2011,7,21).toordinal()+.5
plotFrame(t,ax[0:3])
t=datetime(2011,8,10).toordinal()+.49
plotFrame(t,ax[3:6])

fig.subplots_adjust(top=1,hspace=-.1)
fig.savefig('rho_surface_free.pdf')