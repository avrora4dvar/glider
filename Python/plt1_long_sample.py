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
import datetime
import matplotlib.ticker as ticker

fig=plt.figure(figsize=(7.4,8.92*.5))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})

#%% Load data

mat=sio.loadmat('V:\ipasmans/for_long_sample.mat',squeeze_me=True,struct_as_record=False)
mat=mat['samples']

mat.t=mat.t-366.
mat.t=mat.t+.5
dateRef=datetime.datetime(2005,1,1).toordinal()

ax=[]
for i1 in np.arange(0,4):
    ax1=fig.add_subplot(2,2,i1+1)
    ax1.set_xlim(2392+dateRef,2434+dateRef)
    ax1.grid(color=(.7,.7,.7),linestyle=':',linewidth=.5)
    ax1.plot(np.array([2410,2410])+dateRef,[-100,100],color='k',linewidth=1)
    if i1 >=0:
        ax1.set_xlabel('2011')
        ax1.xaxis.set_major_locator(ticker.FixedLocator(np.arange(2392,2434.1,6)+dateRef))
        ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:plotter.num2time(x).strftime('%m/%d')))
    else:
        ax1.xaxis.set_major_formatter(ticker.NullFormatter())
    if i1==4:
        ax1.set_axis_off()
    ax.append(ax1)


#%% Plot

markers=['*','o','x',None]
colors=['m','r','g','b']
labels=['Surface ','Glider','Combined','No DA']
linewidth=[1,1,1,1]

#SST
pplot=[]
for (val1,label1,marker1,color1,linewidth1) in zip(mat.sst.rms,labels,markers,colors,linewidths):
    print(label1)
    print(np.shape(val1))
    pplot1=ax[0].plot(mat.t,val1,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[0].set_ylim(0,3.5)
    ax[0].set_ylabel(u'RMSE sea-surface temperature [\u00B0C]')
    ax[0].annotate('a)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
    

##SSH
#for (val1,label1,marker1,color1,linewidth1) in zip(mat.ssh.rms,labels,markers,colors,linewidths):
#    pplot1=ax[3].plot(mat.t,val1*100,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
#    ax[3].set_ylim(0,25)
#    ax[3].set_ylabel(r'Sea-surface height [$\mathrm{cm}$]')
#    pplot.append(pplot1)

#uv
for (val1,label1,marker1,color1,linewidth1) in zip(mat.uv.rms,labels,markers,colors,linewidths):
    ax[1].plot(mat.t,val1*100,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[1].set_ylim(0,40)
    ax[1].set_ylabel(r'RMSE HFR velocity [$\mathrm{cms^{-1}}$]')
    ax[1].annotate('b)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
#glider T
for (val1,label1,marker1,color1,linewidth1) in zip(mat.gliderT.rms,labels,markers,colors,linewidths):
    ax[2].plot(mat.t,val1,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[2].set_ylim(0,5)
    ax[2].set_ylabel(u'RMSE glider temperature [\u00B0C]')
    ax[2].annotate('c)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
#glider s
for (val1,label1,marker1,color1,linewidth1) in zip(mat.gliderS.rms,labels,markers,colors,linewidths):
    ax[3].plot(mat.t,val1,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[3].set_ylim(0,1.4)
    ax[3].set_ylabel(r'RMSE glider salinity [$\mathrm{ppt}$]')
    ax[3].annotate('d)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
   
ax[2].legend(loc='lower left',bbox_to_anchor=(.01,.4))
#fig.subplots_adjust(left=.06,right=.94)
fig.tight_layout()
fig.savefig('rms_for.pdf')

