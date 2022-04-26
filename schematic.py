# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 17:26:14 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.path as mplPath

pl.close('all')
fig, ax = pl.subplots(1,2)

ax1 = pl.subplot(121)

ax1.axhline(y=0.5,ls='--',color='k')
pl.ylim(-90,90); pl.yticks([-90,-60,-30,0,30,60,90],fontsize=16)

IZ1 = Polygon([(0.1,-5),(0.1,5),(0.9,5),(0.9,-5)],fill=True,fc='b',
               closed=True)
pl.gca().add_patch(IZ1)

ST1 = Polygon([(0.1,5),(0.1,30),(0.9,30),(0.9,5)],fill=True,fc='r',
               closed=True)
ST2 = Polygon([(0.1,-5),(0.1,-30),(0.9,-30),(0.9,-5)],fill=True,fc='r',
               closed=True)
#
pl.gca().add_patch(ST1); pl.gca().add_patch(ST2)
#
SP1 = Polygon([(0.1,30),(0.1,60),(0.9,60),(0.9,30)],fill=True,fc='b',
               closed=True)
SP2 = Polygon([(0.1,-30),(0.1,-60),(0.9,-60),(0.9,-30)],fill=True,fc='b',
               closed=True)

pl.gca().add_patch(SP1); pl.gca().add_patch(SP2)
pl.title('(a)',loc='center',fontsize=17)
pl.ylabel('latitude',fontsize=20)
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.annotate('$P-E=0$',(0.3,-80),fontsize=20)

ax2 = pl.subplot(122)

ax2.axhline(y=0.5,ls='--',color='k')

ax2.axhline(y=0.5,ls='--',color='k')
pl.ylim(-90,90); pl.yticks([-90,-60,-30,0,30,60,90])

IZ2 = Polygon([(0.1,-5),(0.1,5),(0.9,5),(0.9,-5)],fill=True,fc='b',
               closed=True)
pl.gca().add_patch(IZ2)

ST3 = Polygon([(0.1,5),(0.1,30),(0.9,30),(0.9,5)],fill=True,fc='r',
               closed=True)
ST4 = Polygon([(0.1,-5),(0.1,-30),(0.9,-30),(0.9,-5)],fill=True,fc='r',
               closed=True)
#
pl.gca().add_patch(ST3); pl.gca().add_patch(ST4)
#
SP3 = Polygon([(0.1,30),(0.1,60),(0.9,60),(0.9,30)],fill=True,fc='b',
               closed=True)
SP4 = Polygon([(0.1,-30),(0.1,-60),(0.9,-60),(0.9,-30)],fill=True,fc='aliceblue',
               closed=True)
pl.gca().add_patch(SP3); pl.gca().add_patch(SP4)
pl.title('(b)',loc='center',fontsize=17)
ax2.yaxis.set_major_formatter(pl.NullFormatter())
ax2.xaxis.set_major_formatter(pl.NullFormatter())
pl.annotate('$P-E<0$',(0.3,-80),fontsize=20)