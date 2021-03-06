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
import seawater as sw
imp.reload(plotter)

tOut=[datetime(2011,7,21,0,0,0),datetime(2011,8,10,13,0,0)]
tOut=plotter.time2num(tOut)
tOut=np.arange(tOut[0],tOut[1],2./24.)



#%% Read data

mat=sio.loadmat('V:\ipasmans\cross_NH10_ana.mat',squeeze_me=True,struct_as_record=False)
model=mat['model']

#Adjust times to Python
for model1 in model:
    model1.t=model1.t-366.
dateRef=datetime(2005,1,1).toordinal()

#Read grid
grd=plotter.read_grid_rho()

#Obslist
obs=plotter.read_obslist()

model[0].name='Surface'
model[1].name='Glider'
model[2].name='Combined'
model[3].name='Glider T'
model[4].name='Free'

model[0].label='a)'
model[1].label='b)'
model[2].label='c)'
model[3].label='Glider T'
model[4].label='Free'

model=model[[0,1,2]]

#%% Layout

plt.close('all')
fig=plt.figure(figsize=(5.61,8.92*.4))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})


#Colorbar
cmap=plt.cm.get_cmap('plasma')


#Plots

def plotAx():
    ax=[]
    for i in np.arange(0,3):
        ax1=fig.add_subplot(3,1,i+1)

        #Limits
        #ax1.set_title(model[i].name)        
        ax1.set_ylim(300,0)
        ax1.set_xlim(-126.,-124.)
        
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(base=.5))
        if i==2:
            ax1.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.1f}'))
            ax1.set_xlabel('Longitude')
        else:
            ax1.xaxis.set_major_formatter(ticker.NullFormatter())
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(base=100))
        ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
        ax1.set_ylabel(r'Depth [m]')
        
        ax.append(ax1)
    return ax
    
#Fill
im=[]
def plotFrame(t):
    fig.clf()
    ax=plotAx()
    cplot=[]; gplot=[];  bcplot=[]
    for model1,ax1 in zip(model,ax):
        #Time
        #fig.suptitle(plotter.num2time(t).strftime('%Y-%b-%d %H:%M'))
        fig.subplots_adjust(bottom=0.15,top=.95,right=.85)
        
        it=np.interp(t,model1.t,np.arange(0,len(model1.t)))
        salt1=model1.salt[:,:,int(it)]*(1.-it%1)+model1.salt[:,:,int(it)+1]*(it%1)
        temp1=model1.temp[:,:,int(it)]*(1.-it%1)+model1.temp[:,:,int(it)+1]*(it%1)
        z1=model1.z[:,:,int(it)]*(1.-it%1)+model1.z[:,:,int(it)+1]*(it%1)
        rho1=sw.pden(salt1,temp1,-z1,0.0*z1 )
        rho1=rho1-1e3
        
        z2,rho2=plotter.z_interp(z1,rho1,np.array([-300.,0]))
        print([np.nanmax(rho2),np.nanmin(rho2)])
        ax1.contourf(np.reshape(model1.lon[:,0],(-1,1))+np.zeros(np.shape(z2)),-z2,rho2,[tick1 for tick1 in np.arange(18.,28.01,.5)],cmap=cmap,extend='neither')
        cplot1=ax1.contourf(np.reshape(model1.lon[:,0],(-1,1))+np.zeros(np.shape(z2)),-z2,rho2,[tick1 for tick1 in np.arange(18.,28.01,.5)],cmap=cmap,extend='neither')

        bcplot1=ax1.contour(model1.lon,-z1,salt1,levels=[31.5],colors='k',linewidths=.5,linestyles='-')
        bcplot.append(bcplot1)
        bcplot2=ax1.contour(model1.lon,-z1,rho1,levels=[26.5],linestyles='--',colors='k',linewidths=.5)
        bcplot.append(bcplot2)
        
        #Coast   
        h=-z1[:,0]; h[np.isnan(h)]=3.
        h=np.concatenate(([10e3],h,[10e3]))
        hlon=np.concatenate(([-135.],model1.lon[:,0],[-120.]))
        p=patch.Polygon(np.column_stack((hlon,h)),facecolor=(.5,.6,.5),edgecolor=None,closed=True)
        ax1.add_patch(p)
    
        tLim1=3.*int((t-2)/3.)+np.array([2.,5.])
        obsll=[(lon1,z1) for (lon1,z1,type1,t1) in zip(obs['lon'],obs['z'],obs['type'],obs['t']) if type1==6 and tLim1[0]<=t1<=tLim1[1]]
        lon1,z1=zip(*obsll)
        gplot1=ax1.plot([np.min(lon1),np.min(lon1),np.max(lon1),np.max(lon1)],
                         [np.min(z1),np.max(z1),np.max(z1),np.min(z1)],':',
                         color=(.5,.5,.5),linewidth=1.)
        
        #Plot label
        tplot=ax1.text(-124.05,280,model1.label)
        #t.set_bbox(dict(facecolor='w', alpha=1, edgecolor='w'))
        
        #Colorbar
        cax=fig.add_axes([.855,.15,.02,.8])
        cbar=fig.colorbar(cplot1,cax=cax,orientation='vertical',spacing='vertical')
        cbar.set_ticks([tick1 for tick1 in np.arange(18.0,28.01,1.)])
        cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.0f}'.format(x))
        cbar.update_ticks()
        cbar.set_label(r'Pot. density [$\mathrm{kg m^{-3}}$]')
        
        cplot.append(cplot1)
        gplot.append(gplot1)
        
        #plt.tight_layout()
        
    return cplot,gplot,bcplot

#%%

t=datetime(2011,8,8,00,0,0).toordinal()+.5
plotFrame(t)

#%%

fig.savefig('cross_rho_ana.pdf',figsize=(7.48,5))