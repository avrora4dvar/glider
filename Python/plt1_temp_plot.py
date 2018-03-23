# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:39:11 2018

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

imp.reload(plotter)

tOut=[datetime(2011,7,21,12,0,0),datetime(2011,8,13,12,0,0)]
tOut=plotter.time2num(tOut)
tOut=np.arange(tOut[0],tOut[1],2./24.)



#%% Read data

mat=[]
mat.append(sio.loadmat('V:\ipasmans\\Exp40ana_avg.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:\ipasmans\Exp35for_avg.mat',squeeze_me=True,struct_as_record=False))
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

mat[0]['name']='a) Observations'
mat[1]['name']='b) Surface Only'
mat[2]['name']='c) Glider Only'
mat[3]['name']='d) Combined'

#%% Layout

plt.close('all')
fig=plt.figure(figsize=(7,8.24))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})


#Colorbar
cmap=plt.cm.get_cmap('jet')



#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,4):
        ax1=fig.add_subplot(2,2,i+1)
        plotter.add_bathy(ax1)
        plotter.add_coast(ax1)
        #plotter.layout_surface(ax1,xLim=[-127,-123.5],yLim=[43,47])
        plotter.layout_surface(ax1,xLim=[-129,-123],yLim=[41,49])
        
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(2))
        
        if np.mod(i,2)==0:
            ax1.set_ylabel(r'Latitude')
        else:
            ax1.yaxis.set_major_formatter(ticker.NullFormatter())
        if i>=2:
            ax1.set_xlabel(r'Longitude')
        else:
            ax1.xaxis.set_major_formatter(ticker.NullFormatter())
        #t=ax1.text(-128.6,48.5,mat[i]['name'])
        #t.set_bbox(dict(facecolor='w', alpha=1, edgecolor='w'))

        ax.append(ax1)
    return ax


    
#%%
    
t=datetime(2011,8,4).toordinal()+0.5
tLim1=3.*int((t-2)/3)+np.array([2.,5.])
ax=plotAx()

def plotFrame(t):
    cplot=[]; gplot=[]; qplot=[]
    for mat1,ax1 in zip(mat[1:],ax[1:]):
        #fig.suptitle(plotter.num2time(t).strftime('%Y-%b-%d %H:%M'))
        
        it=np.interp(t,mat1['t'],np.arange(0,len(mat1['t'])))
        temp1=mat1['temp'][:,:,int(it)]*(1.-it%1)+mat1['temp'][:,:,int(it)+1]*(it%1)
        print([np.nanmin(temp1),np.nanmax(temp1)])
        cplot1=ax1.contourf(mat1['lon'],mat1['lat'],temp1,[tick1 for tick1 in np.arange(10,20.1)],cmap=cmap,extend='neither')
        
        salt1=mat1['salt'][:,:,int(it)]*(1.-it%1)+mat1['salt'][:,:,int(it)+1]*(it%1)
        print([np.nanmin(salt1),np.nanmax(salt1)])
        qplot1=ax1.contour(mat1['lon'],mat1['lat'],salt1,levels=[31.5],colors='k',linestyles=['-'])
    
    
        obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
        lon1,lat1=zip(*obsll)
        gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=5,fillstyle='full')
        
        ax1.set_title(mat1['name'])
        
        #Colorbar
        fig.subplots_adjust(right=.8)
        cax=fig.add_axes([.82,.15,.02,.7])
        cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional')
        cbar.set_ticks([tick1 for tick1 in np.arange(10.,20.1)])
        cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.0f}'.format(x))
        cbar.update_ticks()
        cbar.set_label(r'Temperature [$^\circ\mathrm{C}$]')
        
        cplot.append(cplot1)
        gplot.append(gplot1)
        qplot.append(qplot1)
        
        #plt.tight_layout()
        
    return cplot,qplot,gplot

#Observations
ax1=ax[0]
obsll=[(lon1,lat1,val1) for (lon1,lat1,val1,type1,t1) in zip(obs['lon'],obs['lat'],obs['val'],obs['type'],obs['t']) if type1==4 and t-1.5<=t1<=t+1.5]
lon1,lat1,val1=zip(*obsll)
ax1.scatter(lon1,lat1,c=val1,vmin=10,vmax=20,cmap=cmap,s=.02)
obsll=[(lon1,lat1,val1) for (lon1,lat1,val1,type1,t1) in zip(obs['lon'],obs['lat'],obs['val'],obs['type'],obs['t']) if type1==4 and t-.5<=t1<=t+.5]
lon1,lat1,val1=zip(*obsll)
ax1.scatter(lon1,lat1,c=val1,vmin=10,vmax=20,cmap=cmap,s=.02)
obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
lon1,lat1=zip(*obsll)
gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=5,fillstyle='full')
ax1.set_title(mat[0]['name'])

#time
plotFrame(t)


#%%

fig.subplots_adjust(bottom=.05,top=.96)
fig.savefig('temp_20110804.pdf')