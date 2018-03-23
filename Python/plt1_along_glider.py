# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:08:06 2018

@author: ipasmans
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:34:36 2018

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
from matplotlib import patches
imp.reload(plotter)

plt.close('all')
fig=plt.figure(figsize=(7.4,8.2))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size':8})

#%% Load

mat=sio.loadmat('V:/ipasmans/along_glider.mat',squeeze_me=True,struct_as_record=False)

#Time to Python
mat['obs'].t=mat['obs'].t-366
    
def dateFormatter(x,pos):
    date=plotter.num2time(x).strftime('%m/%d')
    return date

tLim=plotter.time2num([datetime(2011,7,21),datetime(2011,8,11)])

#%% Create figure

ax=[]
alfabet=['a)','b)','c)','d)','e)']
for iax in np.arange(0,5):
    ax1=fig.add_subplot(5,1,iax+1)
    
    #y
    ax1.set_ylim(200,0)
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax1.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
    ax1.set_ylabel('Depth [m]')
    
    #Text
    ax1.text(tLim[-1]-.9,190,alfabet[iax])
    
    #x
    ax1.set_xlim(tLim[0],tLim[1])
    ax1.xaxis.set_major_locator(ticker.FixedLocator(np.arange(tLim[0],tLim[1]+.001,3)))
    if iax==4:
        ax1.xaxis.set_major_formatter(ticker.FuncFormatter(dateFormatter))
        ax1.set_xlabel('2011')
    else:
        ax1.xaxis.set_major_formatter(ticker.NullFormatter())
        
    ax.append(ax1)
    
#%% Plot observations
    
cm=plt.cm.get_cmap('jet')

ax1=ax[0]
t1=np.reshape(mat['obs'].t,(-1,1))+np.zeros(np.shape(mat['obs'].temp))
z1=np.reshape(mat['obs'].depth,(1,-1))+np.zeros(np.shape(mat['obs'].temp))
val1=mat['obs'].temp
cplot1=ax1.contourf(t1,z1,val1,levels=np.arange(5,18.5),cmap=cm,vmin=5,vmax=18,extend='neither')

t1=np.reshape(mat['obs'].t,(-1,1))+np.zeros(np.shape(mat['obs'].temp))
z1=np.reshape(mat['obs'].depth,(1,-1))+np.zeros(np.shape(mat['obs'].temp))
val1=mat['obs'].salt
cplot2=ax1.contour(t1,z1,val1,linestyles=[(0,(1,10)),(0,(1,5)),(0,(1,1)),(0,(5,10)),(0,(5,5)),'-',(0, (3, 1, 1, 1, 1, 1)),(0, (3, 10, 1, 10, 1, 10))],linewidths=.7,
            colors=['k','k','k','k','k','k','k','k'],levels=[27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5])

#%% Plot models

levels=[27.5,28.5,29.5,30.5,31.5,32.5,33.5,34.5]
for (mat1,ax1) in zip(mat['model'],ax[1:]):

    t1=np.reshape(mat['obs'].t,(-1,1))+np.zeros(np.shape(mat['obs'].temp))
    z1=np.reshape(mat['obs'].depth,(1,-1))+np.zeros(np.shape(mat['obs'].temp))
    val1=mat1.temp
    cplot1=ax1.contourf(t1,z1,val1,levels=np.arange(5,18,.5),cmap=cm,vmin=5,vmax=18,extend='neither')

    t1=np.reshape(mat['obs'].t,(-1,1))+np.zeros(np.shape(mat['obs'].temp))
    z1=np.reshape(mat['obs'].depth,(1,-1))+np.zeros(np.shape(mat['obs'].temp))
    val1=mat1.salt
    cplot2=ax1.contour(t1,z1,val1,linestyles=[(0,(1,10)),(0,(1,5)),(0,(1,1)),(0,(5,10)),(0,(5,5)),'-',(0, (3, 1, 1, 1, 1, 1)),(0, (3, 10, 1, 10, 1, 10))],linewidths=.7,
            colors=['k','k','k','k','k','k','k','k'],levels=levels)
    
#%% Colorbar
for col1,level1 in zip(cplot2.collections,levels):
    col1.set_label('{:.1f} ppt'.format(level1))
leg=ax[-1].legend(loc='upper right',bbox_to_anchor=(0.05,.52))
leg.get_frame().set_alpha(1)

cax=fig.add_axes([.2,.05,.7,.015])
cbar=fig.colorbar(cplot1,cax=cax,orientation='horizontal',spacing='uniform')
cbar.set_ticks([tick1 for tick1 in np.arange(5.0,18.01,1.)])
cbar.formatter=ticker.FuncFormatter(lambda x,pos:'{:.0f}'.format(x))
cbar.update_ticks()
cbar.set_label(u'Temperature [\u00B0C]')

#%% Save
fig.subplots_adjust(top=.98)
fig.savefig('along_glider.pdf')
