# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 17:38:19 2018

@author: np838619
"""

from __future__ import division
import pylab as pl

def Xlabels(rlspts,ax,s,l):
    """
    """
    a = pl.arange(0,rlspts.shape[0],s); b = pl.around(rlspts[a,l],decimals=0)
    ax.set_xticks(a); ax.set_xticklabels(b.astype('int'),fontsize=14)

pl.close('all')
sheddir = '/home/np838619/Watershed/'

shed = 'eaa'; c = 'g'
shed_sns = pl.genfromtxt(sheddir+'variability/seasonal_profs_'+shed+'.csv')
rlspts = pl.genfromtxt(sheddir+'shed_defs/EAA_traj_release_new.txt',skip_header=5)

fig,ax = pl.subplots(2,2,figsize=(10,8))

ax1 = pl.subplot(221)
ax1.plot(shed_sns[:,0],color=c,lw=2)
pl.ylim(-350,100); pl.xlim(0,rlspts.shape[0])
#ax1.set_yticks(pl.linspace(-100,100,9))
ax1.grid(axis='y'); pl.ylabel('kg/m/s',fontsize=16)
pl.text(10,-90,'(a) DJF',fontsize=20); pl.yticks(fontsize=13)
ax1.xaxis.set_major_formatter(pl.NullFormatter())

ax2 = pl.subplot(222)
ax2.plot(shed_sns[:,1],color=c,lw=2)
pl.ylim(-350,100); pl.xlim(0,rlspts.shape[0])
#ax2.set_yticks(pl.linspace(-100,100,9))
ax2.yaxis.set_major_formatter(pl.NullFormatter())
ax2.xaxis.set_major_formatter(pl.NullFormatter())
ax2.grid(axis='y')
pl.text(10,-90,'(b) MAM',fontsize=20)

ax3 = pl.subplot(223)
ax3.plot(shed_sns[:,2],color=c,lw=2)
pl.ylim(-350,100); pl.xlim(0,rlspts.shape[0])
ax3.grid(axis='y')#; ax3.set_yticks(pl.linspace(-100,100,9))
pl.text(10,-90,'(c) JJA',fontsize=20); pl.yticks(fontsize=13)
Xlabels(rlspts,ax3,30,0); ax3.tick_params(axis='x',pad=7)
pl.xlabel('latitude',fontsize=16); pl.ylabel('kg/m/s',fontsize=16)

ax4 = pl.subplot(224)
ax4.plot(shed_sns[:,3],color=c,lw=2)
pl.ylim(-350,100); pl.xlim(0,rlspts.shape[0])
ax4.yaxis.set_major_formatter(pl.NullFormatter())
ax4.grid(axis='y')#; ax4.set_yticks(pl.linspace(-100,100,9))
pl.text(10,-90,'(d) SON',fontsize=20)
Xlabels(rlspts,ax4,30,0); ax4.tick_params(axis='x',pad=7)
pl.xlabel('latitude',fontsize=16)

pl.tight_layout()
#pl.savefig(sheddir+'variability/seasonal_profs_'+shed+'.png')