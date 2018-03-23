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

plt.close('all')
fig=plt.figure(figsize=(7.4,8.92*.5))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})

#%% Load data

mat=sio.loadmat('V:\ipasmans/glider_variation_ana.mat',squeeze_me=True,struct_as_record=False)
mat=mat['samples']

mat.t=mat.t-366.
mat.t=mat.t
dateRef=datetime.datetime(2005,1,1).toordinal()

ax=[]
for i1 in np.arange(0,4):
    ax1=fig.add_subplot(2,2,i1+1)
    ax1.set_xlim(2392+dateRef,2413+dateRef)
    ax1.grid(color=(.7,.7,.7),linestyle=':',linewidth=.5)
    #ax1.set_ylim(0,6)
    if i1 >=0:
        ax1.set_xlabel('2011')
        ax1.xaxis.set_major_locator(ticker.FixedLocator(np.arange(2392,2434.1,3)+dateRef))
        ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:plotter.num2time(x).strftime('%m/%d')))
    else:
        ax1.xaxis.set_major_formatter(ticker.NullFormatter())
    if i1==4:
        ax1.set_axis_off()
    ax.append(ax1)
    
#%% Differences of the whole period

def rms(x):
    return np.sqrt(np.nanmean(x*x))

print('sst')
for val1 in mat.sst.rms:
    print(rms(val1)/rms(mat.sst.rms[0,:]))
print('uv')
for val1 in mat.uv.rms:
    print(rms(val1)/rms(mat.uv.rms[0,:]))
print('gliderT')
for val1 in mat.gliderT.rms:
    print(rms(val1)/rms(mat.gliderT.rms[0,:]))
print('gliderS')
for val1 in mat.gliderS.rms:
    print(rms(val1)/rms(mat.gliderS.rms[0,:]))




#%% Plot

markers=['o','*','x',None]
colors=['r','b','g','k']
labels=['Glider Only','Glider T','No HFR,SSH','No SST,SSH']
linewidths=[1,1,1,1]

#SST
pplot=[]
mat.sst.rms=mat.sst.rms[[1,2,4,5],:]/np.reshape(mat.sst.rms[0,:],(1,-1))
for (val1,label1,marker1,color1,linewidth1) in zip(mat.sst.rms,labels,markers,colors,linewidths):
    print(label1)
    print(np.shape(val1))
    pplot1=ax[0].plot(mat.t,val1,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[0].set_ylim(0,7)
    ax[0].set_ylabel(u'Sea-surface temperature')
    ax[0].annotate('a)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
    

##SSH
#for (val1,label1,marker1,color1,linewidth1) in zip(mat.ssh.rms,labels,markers,colors,linewidths):
#    pplot1=ax[3].plot(mat.t,val1*100,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
#    ax[3].set_ylim(0,25)
#    ax[3].set_ylabel(r'Sea-surface height [$\mathrm{cm}$]')
#    pplot.append(pplot1)

#uv
mat.uv.rms=mat.uv.rms[[1,2,4,5],:]/np.reshape(mat.uv.rms[0,:],(1,-1))
for (val1,label1,marker1,color1,linewidth1) in zip(mat.uv.rms,labels,markers,colors,linewidths):
    ax[1].plot(mat.t,val1,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[1].set_ylabel(r'HFR velocity')
    ax[1].set_ylim(0,4)
    ax[1].annotate('b)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
#glider T
mat.gliderT.rms=mat.gliderT.rms[[1,2,4,5],:]/np.reshape(mat.gliderT.rms[0,:],(1,-1))
for (val1,label1,marker1,color1,linewidth1) in zip(mat.gliderT.rms,labels,markers,colors,linewidths):
    ax[2].plot(mat.t,val1,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[2].set_ylim(0,4)
    ax[2].set_ylabel(u'Glider temperature')
    ax[2].annotate('c)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
#glider s
mat.gliderS.rms=mat.gliderS.rms[[1,2,4,5],:]/np.reshape(mat.gliderS.rms[0,:],(1,-1))
for (val1,label1,marker1,color1,linewidth1) in zip(mat.gliderS.rms,labels,markers,colors,linewidths):
    ax[3].plot(mat.t,val1,color=color1,marker=marker1,label=label1,markevery=3,linewidth=linewidth1,markersize=6)
    ax[3].set_ylim(0,4)
    ax[3].set_ylabel(r'Glider salinity ')
    ax[3].annotate('d)',xy=(.02,.93),xytext=(.02,.93),textcoords='axes fraction',xycoords='axes fraction')
    
    
   
ax[3].legend(loc='upper left',bbox_to_anchor=(.05,.99),ncol=2,columnspacing=1)
#fig.subplots_adjust(left=.06,right=.94)
fig.tight_layout()
fig.savefig('rms_ratio_ana.pdf')

