# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:39:11 2018

@author: ipasmans
"""

import matplotlib
from matplotlib import dates as dates
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
mat.append(sio.loadmat('V:\ipasmans\\Exp37ana_avg.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:\ipasmans\Exp38ana_avg.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:\ipasmans\Exp41ana_avg.mat',squeeze_me=True,struct_as_record=False))
mat.append(sio.loadmat('V:\ipasmans\Exp42ana_avg.mat',squeeze_me=True,struct_as_record=False))

#Adjust times to Python
for mat1 in mat:
    mat1['t']=mat1['t']-366.
dateRef=datetime(2005,1,1).toordinal()

#Read grid
grd=plotter.read_grid_rho()

#Obslist
obs=plotter.read_obslist()

mat[0]['name']='a) Combined'
mat[1]['name']='b) Glider T'
mat[2]['name']='c) Glider+SST'
mat[3]['name']='d) Glider+HFR'

#%% Layout

plt.close('all')
fig=plt.figure(figsize=(7,8.24))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})


#Colorbar
cmap=plt.cm.get_cmap('seismic')



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
    


def plotFrame(t):
    fig.clf()
    ax=plotAx()
    cplot=[]; gplot=[]; qplot=[]
    for mat1,ax1 in zip(mat,ax):
        #Time
        #fig.suptitle(plotter.num2time(t).strftime('%Y-%b-%d %H:%M'))
        
        it=np.interp(t,mat1['t'],np.arange(0,len(mat1['t'])))
        vort1=mat1['vort'][:,:,int(it)]*(1.-it%1)+mat1['vort'][:,:,int(it)+1]*(it%1)
        cplot1=ax1.contourf(mat1['lon'],mat1['lat'],vort1,[tick1 for tick1 in np.arange(-1.4e-4,1.401e-4,.2e-4) if np.abs(tick1)>1e-8],cmap=cmap,extend='both')
        
        u1=mat1['u'][:,:,int(it)]*(1.-it%1)+mat1['u'][:,:,int(it)+1]*(it%1)
        v1=mat1['v'][:,:,int(it)]*(1.-it%1)+mat1['v'][:,:,int(it)+1]*(it%1)
        qplot1=plotter.add_uv2d(ax1,mat1['lon'],mat1['lat'],u1,v1)
    
        tLim1=np.empty(2)
        tLim1[0]=dates.date2num(datetime(2011,7,21))
        tLim1[1]=tLim1[0]+np.ceil((t-tLim1[0])/3.)*3.
        obsll=[(lon1,lat1) for (lon1,lat1,type1,t1) in zip(obs['lon'],obs['lat'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
        lon1,lat1=zip(*obsll)
        gplot1=ax1.plot(lon1,lat1,'.',color=(1.,.55,0.),markersize=2,fillstyle='full')
        
        ax1.text(-124.25,43.0,r'1 ms$^{-1}$')
        
        ax1.set_title(mat1['name'])
        
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

#time
t=datetime(2011,8,4).toordinal()+0.5
plotFrame(t)


#%%

fig.subplots_adjust(bottom=.05,top=.96)
fig.savefig('vort_exp38_20110804.pdf')
