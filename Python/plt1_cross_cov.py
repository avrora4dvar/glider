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
import matplotlib.ticker as ticker

plt.rcParams.update({'font.size':12})


#%% Figure

#Obslist
obs=plotter.read_obslist()

mat=sio.loadmat('V:\ipasmans\ocean2018\cov_bal.mat',squeeze_me=True,struct_as_record=False)
cov=mat['cov']; grd=mat['grd']; loc=mat['loc']

#%% Map

plt.close('all')
fig=plt.figure(figsize=(4.5,2.5),dpi=300)
plt.rcParams.update({'font.size': 12})


#Colorbar
cmap=plt.cm.get_cmap('YlGnBu_r')


#Plots
def plotAx():
    ax=[]
    for i in np.arange(0,1):
        ax1=fig.add_subplot(1,1,i+1)
        #plotter.layout_surface(ax1,xLim=[-127,-123.5],yLim=[43,47])
        ax1.set_ylabel(r'Depth [m]')
        ax1.set_xlabel(r'Longitude')
        ax1.set_ylim([100.,0.])
        ax1.set_xlim([-125.5,-124])

        ax1.xaxis.set_major_locator(ticker.MultipleLocator(base=1))
        ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=50))
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
      
        ax.append(ax1)
    return ax
ax=plotAx()
ax=ax[0]

#%% Create

fig.subplots_adjust(left=.14,right=.75,bottom=.2,top=.95)

z1=np.column_stack((grd.zw[:,grd.ilat-1,0],grd.zr[:,grd.ilat-1,:],grd.zw[:,grd.ilat-1,-1]))
salt2=np.squeeze(cov.salt[:,grd.ilat-1,:]); salt2=np.column_stack((salt2[:,0],salt2,salt2[:,-1]))
temp2=np.squeeze(cov.temp[:,grd.ilat-1,:]); temp2=np.column_stack((temp2[:,0],temp2,temp2[:,-1]))
z2,salt2=plotter.z_interp(z1,salt2,np.array([-300.,0]))
z2,temp2=plotter.z_interp(z1,temp2,np.array([-300.,0]))

levelsf=np.arange(-.16,.001,.01)
levels=np.arange(0.1,1.01,.1)


lon2=np.reshape(grd.lon[:,1],(-1,1))+np.zeros(np.shape(z2))
cplot1=ax.contourf(lon2,-z2,salt2,levels=levelsf,cmap=cmap,extend='min')
cplot2=ax.contour(lon2,-z2,temp2,levels=levels,linewidths=1,colors='k')

#Point
ax.plot(loc.lon,0,'ko',markerfacecolor='k',markersize=4)

#Coast   
h=-z1[:,0]; h[np.isnan(h)]=3.
h=np.concatenate(([10e3],h,[10e3]))
hlon=np.concatenate(([-135.],grd.lon[:,0],[-120.]))
p=patch.Polygon(np.column_stack((hlon,h)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
ax.add_patch(p)

#Velocity
if False:
    step=2
    u1=np.squeeze(.5*cov.u[0:-2,grd.ilat-1,:]+.5*cov.u[1:-1,grd.ilat-1,:]); 
    u1=np.column_stack((u1[:,0],u1,u1[:,-1]))
    v1=np.squeeze(.5*cov.v[:,grd.ilat-2,:]+.5*cov.v[:,grd.ilat-1,:])
    v1=np.column_stack((v1[:,0],v1,v1[:,-1]))
    z2,v2=plotter.z_interp(z1,v1,np.array([0,-300]))
    v3=v2+0.0; v3[v3<.005]=np.nan
    plt.scatter(lon2[::step,::step],-z2[::step,::step],s=v3[::step,::step]/.001,c='k',marker='o')
    v3=-v2+0.0; v3[v3<.005]=np.nan
    plt.scatter(lon2[::step,::step],-z2[::step,::step],s=v3[::step,::step]/.001,c='k',marker='x')

#Colorbar
cax=fig.add_axes([.76,.2,.02,.75])
cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='proportional')
cbar.set_ticks(levelsf[[0,-1]])
cbar.locator=ticker.MultipleLocator(base=.02)
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.2f}'.format(x))
cbar.update_ticks()
cbar.set_label(r'T,S-covariance [$\mathrm{ppt ^\circ C}$]')

#%%

plt.savefig('cross_cov_transect.pdf')

#qplot=ax.quiver(lon2[::step,::step],-z2[::step,::step],
#          0.*v2[::step,::step],v2[::step,::step],angles='uv',units='width',scale=1/10)








