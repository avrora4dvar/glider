# -*- coding: utf-8 -*-
"""
Created on Tue May 22 13:53:35 2018

Plot

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
fig=plt.figure(figsize=(5.61,8.92*.25))
matplotlib.rcParams['font.family']='Arial'
plt.rcParams.update({'font.size': 8,'axes.linewidth':1})

#%% Times

tList=[]
tList.append(plotter.time2num(datetime.datetime(2011,7,9)))
tList.append(plotter.time2num(datetime.datetime(2011,7,21)))
tList.append(plotter.time2num(datetime.datetime(2011,8,8)))
tList.append(plotter.time2num(datetime.datetime(2011,8,11)))
tList.append(plotter.time2num(datetime.datetime(2011,9,1)))

tTick=np.arange(tList[0],tList[-1]+.001,3)

#%% 

ax=fig.add_subplot(111)

ax.set_xlim(tList[0],tList[-1])
ax.xaxis.set_major_locator(ticker.FixedLocator(np.take(tList,[0,1,2,4])))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos:plotter.num2time(x).strftime('%m/%d')))
ax.yaxis.set_major_locator(ticker.NullLocator())
ax.yaxis.set_major_formatter(ticker.NullFormatter())
ax.set_xlabel('2011')

ax.grid(linestyle='--',color='k',linewidth=1)

#%% Arrows

#Colors
cT=(1.,.55,0.)
cF='k'

#Surface
p1=ax.arrow(tList[0],.8,tList[1]-tList[0],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cT,ec=cT)
ax.arrow(tList[1],.8,tList[3]-tList[1],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cT,ec=cT)
p2=ax.arrow(tList[3],.8,tList[-1]-tList[3],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cF,ec=cF)
#Glider
ax.arrow(tList[1],.6,tList[3]-tList[1],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cT,ec=cT)
ax.arrow(tList[3],.6,tList[-1]-tList[3],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cF,ec=cF)
#Combined
ax.arrow(tList[1],.4,tList[3]-tList[1],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cT,ec=cT)
ax.arrow(tList[3],.4,tList[-1]-tList[3],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cF,ec=cF)
#No DA
ax.arrow(tList[1],.2,tList[-1]-tList[1],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cF,ec=cF)
ax.arrow(tList[3],.2,tList[-1]-tList[3],0,linewidth=1.5,head_width=.05,head_length=2.5,length_includes_head=True,fc=cF,ec=cF)

#%% Text

ax.text(tList[1]+1,.82,'Surface Only')
ax.text(tList[1]+1,.62,'Glider Only')
ax.text(tList[1]+1,.42,'Combined')
ax.text(tList[1]+1,.22,'No DA')

ax.text(tList[3],1.01,'Last DA correction',horizontalalignment='center')
ax.text(tList[1]-1.4,.6,'{}\n{}'.format('Same initial conditions','all cases 7/21'),rotation=90,horizontalalignment='center')

leg=ax.legend([p1,p2],['Analysis+Forecast','Forecast Only'],loc='lower right',ncol=2)
leg.get_frame().set_alpha(1)

#%% Print

fig.savefig('timeline.pdf')