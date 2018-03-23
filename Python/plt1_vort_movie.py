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

imp.reload(plotter)

tOut=[datetime(2011,7,21,12,0,0),datetime(2011,8,13,12,0,0)]
tOut=plotter.time2num(tOut)
tOut=np.arange(tOut[0],tOut[1],2./24.)



#%% Read data

mat=[]
mat.append(sio.loadmat('V:\ipasmans\\Exp40ana_avg.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:\ipasmans\Exp36for_avg.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:\ipasmans\Exp37for_avg.mat',squeeze_me=True,struct_as_record=False))

#Adjust times to Python
for mat1 in mat:
    mat1['t']=mat1['t']-366.
dateRef=datetime(2005,1,1).toordinal()

#Read grid
grd=plotter.read_grid_rho()

#Obslist
obs=plotter.read_obslist()

mat[0]['name']='No DA'
mat[1]['name']='Glider Only'
mat[2]['name']='Combined'

#%% Layout

plt.close('all')
fig=plt.figure(figsize=(11,4.8),dpi=300)
plt.rcParams.update({'font.size': 12})


#Colorbar
cmap=plt.cm.get_cmap('seismic')



#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,3):
        ax1=fig.add_subplot(1,3,i+1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        #plotter.layout_surface(ax1,xLim=[-127,-123.5],yLim=[43,47])
        plotter.layout_surface(ax1,xLim=[-129,-123],yLim=[41,49])
        if i==0:
            ax1.set_ylabel(r'Latitude')
        ax1.set_xlabel(r'Longitude')
        ax1.set_title(mat[i]['name'])
        if i>0:
            ax1.yaxis.set_major_locator(plt.NullLocator())
        ax.append(ax1)
    return ax
    
#Fill
im=[]
def plotFrame(t):
    fig.clf()
    ax=plotAx()
    cplot=[]; gplot=[]; qplot=[]
    for mat1,ax1 in zip(mat,ax):
        #Time
        fig.suptitle(plotter.num2time(t).strftime('%Y-%b-%d %H:%M'))
        
        it=np.interp(t,mat1['t'],np.arange(0,len(mat1['t'])))
        vort1=mat1['vort'][:,:,int(it)]*(1.-it%1)+mat1['vort'][:,:,int(it)+1]*(it%1)
        cplot1=ax1.contourf(mat1['lon'],mat1['lat'],vort1,[tick1 for tick1 in np.arange(-1.4e-4,1.401e-4,.2e-4) if np.abs(tick1)>1e-8],cmap=cmap,extend='both')
        
        u1=mat1['u'][:,:,int(it)]*(1.-it%1)+mat1['u'][:,:,int(it)+1]*(it%1)
        v1=mat1['v'][:,:,int(it)]*(1.-it%1)+mat1['v'][:,:,int(it)+1]*(it%1)
        qplot1=plotter.add_uv2d(ax1,mat1['lon'],mat1['lat'],u1,v1)
    
        tLim1=3.*int((t-2)/3)+np.array([2.,5.])
        obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
        lon1,lat1=zip(*obsll)
        gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=5,fillstyle='full')
        
        #Colorbar
        fig.subplots_adjust(right=.8)
        cax=fig.add_axes([.82,.15,.02,.7])
        cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional')
        cbar.set_ticks([tick1 for tick1 in np.arange(-1.4e-4,1.401e-4,.2e-4)])
        cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.1f}'.format(x*1e4))
        cbar.update_ticks()
        cbar.set_label(r'Rel. vorticity [$10^{-4}\,\mathrm{s^{-1}}$]')
        
        cplot.append(cplot1)
        gplot.append(gplot1)
        qplot.append(qplot1)
        
        #plt.tight_layout()
        
    return cplot,qplot,gplot

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
    print([i,(t-tOut[0])*24])
    return plotFrame(t)


ani= animation.FuncAnimation(fig,aniFunc,frames=nFrame, blit=False,repeat=False)
#ani= animation.ArtistAnimation(fig, im, interval=50, repeat_delay=3000,blit=True)
ani.save('movie.mp4', writer=writer)

#%%
