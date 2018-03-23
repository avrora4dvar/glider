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

#%% Layout

plt.close('all')
fig=plt.figure(figsize=(4.2,5))

plt.rcParams.update({'font.size': 12})


#Colorbar
cmap=plt.cm.get_cmap('seismic')
cgrey=plt.cm.get_cmap('Greys')



#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,1):
        ax1=fig.add_subplot(1,1,i+1)
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
im=[]
def plotFrame(t):
    fig.clf()
    ax=plotAx()
    cplot=[]; gplot=[]; qplot=[]; bcplot=[]
    mat1=mat[2]
    for level1,ax1 in zip(levelList,ax):
        #Time
        fig.suptitle(plotter.num2time(t).strftime('%Y-%b-%d %H:%M'))

        inCon=level1==mat1.contours
        it=np.interp(t,mat1.t,np.arange(0,len(mat1.t)))
        it1=np.minimum(int(it)+1,len(mat1.t)-1)
        print([it,it1])
        
        
        z0=mat[0].z[:,:,inCon,int(it)]*(1-it%1)+mat[0].z[:,:,inCon,it1]*(it%1)
        z1=mat1.z[:,:,inCon,int(it)]*(1-it%1)+mat1.z[:,:,inCon,it1]*(it%1)
        z0=np.squeeze(z0); z1=np.squeeze(z1)
        mask0=np.where(np.isnan(z0),np.ones(np.shape(z0)),np.zeros(np.shape(z0)))
        mask1=np.where(np.isnan(z1-z0),np.ones(np.shape(z1)),np.zeros(np.shape(z1)))
        
        print([np.nanmin(z1),np.nanmax(z1)])
        print([np.nanmin(z1-z0),np.nanmax(z1-z0)])
        cplot1=ax1.contourf(mat1.lon,mat1.lat,z1-z0,[tick1 for tick1 in np.arange(-40.0,40.01,5.0) if np.abs(tick1)>1],cmap=cmap,extend='both')
        
        if level1<=1024.5:
            bcplot1=ax1.contour(mat1.lon,mat1.lat,z1,levels=np.arange(-200,0.1,10.),colors='k',linewidth=.5)
        else:
            bcplot1=ax1.contour(mat1.lon,mat1.lat,z1,levels=np.arange(-2000,0.1,25.),colors='k',linewidth=.5)
        ax1.set_title(r'{:.1f}'.format(float(level1))+r'$\mathrm{kgm^{-3}}$')
        bcplot2=ax1.contourf(mat1.lon,mat1.lat,mask1,levels=[-.5,.5,1],cmap=cgrey,vmin=0.5,vmax=1,alpha=.5)
        
        tLim1=3.*int((t-2)/3)+np.array([2.,5.])
        obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
        lon1,lat1=zip(*obsll)
        gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        
        #Colorbar
        fig.subplots_adjust(right=.8)
        cax=fig.add_axes([.76,.1,.02,.8])
        cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional')
        cbar.set_ticks([tick1 for tick1 in np.arange(-40.0,40.01,10.)])
        cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x))
        cbar.update_ticks()
        cbar.set_label(r'Depth change [$\mathrm{m}$]')
        
        cplot.append(cplot1)
        gplot.append(gplot1)
        bcplot.append(bcplot1)
        
        #plt.tight_layout()
        
    return cplot,gplot,bcplot

#%%

#Animation
# Set up formatting for the movie files
plt.rcParams['animation.ffmpeg_path'] = r'C:\Users\ipasmans\Anaconda3\ffmpeg\bin\ffmpeg.exe'
#FFwriter=animation.FFMpegWriter()
Writer = animation.writers['ffmpeg']
writer = Writer(fps=10, metadata=dict(artist='ivopasmans'), bitrate=1800)
nFrame=len(tOut)

def aniFunc(i):
    t=np.interp(i,[0,nFrame-1],[np.min(tOut),np.max(tOut)])
    return plotFrame(t)


ani= animation.FuncAnimation(fig,aniFunc,frames=nFrame, blit=False,repeat=False)
#ani= animation.ArtistAnimation(fig, im, interval=50, repeat_delay=3000,blit=True)
ani.save('movie.mp4', writer=writer)

#%%
