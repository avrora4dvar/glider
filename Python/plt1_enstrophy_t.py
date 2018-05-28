# -*- coding: utf-8 -*-
"""
Created on Wed May  2 12:02:58 2018

Enstrophy plot

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
fig=plt.figure(figsize=(5.61,8.92*.35))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8})

#%% Load data

mat=sio.loadmat('V:\ipasmans\enstrophy_ana.mat',squeeze_me=True,struct_as_record=False)
mat=mat['model']

for mat1 in mat:
    mat1.t=mat1.t-366
    mat1.slope.t=mat1.slope.t-366

ttick=np.arange(2392,2413.01,3)+plotter.time2num(datetime.datetime(2005,1,1))

#%% Figure

ax=fig.add_subplot(1,1,1)
#x
ax.set_xlim(ttick[[0,-1]])
ax.xaxis.set_major_locator(ticker.FixedLocator(ttick))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: plotter.num2time(x).strftime('%m/%d')))
ax.set_xlabel('2011')
#y
ax.set_ylim(0,225)
ax.yaxis.set_major_locator(ticker.MultipleLocator(25))
ax.set_ylabel(r'Enstrophy [$\mathrm{m^2 s^{-2}}$]')
#Grid
ax.grid(color=(.7,.7,.7),linestyle='--')

#%% 

markers=['*','o','x',None]
colors=['m','r','g','b']
labels=['Surface','Glider','Combined','No DA']
linewidths=[1,1,1,1]

for mat1,colors1,markers1,titles1 in zip(mat,colors,markers,labels):
    if titles1=='Glider Fixed':
        continue
    ax.plot(mat1.t,mat1.enstrophy,color=colors1,label=titles1,marker=markers1,markevery=24)

#Legend
leg=ax.legend(loc='upper left')
leg.get_frame().set_alpha(1)

#%% Save

#fig.savefig('enstrophy.png',dpi=600)
fig.savefig('enstrophy.pdf')
