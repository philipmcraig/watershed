# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 13:40:59 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans

def Lines(m):
    """
    """
    line1 = pl.genfromtxt(sheddir+'shed_defs/'+'Am_traj_release_new.txt',skip_header=5)
    line2 = pl.genfromtxt(sheddir+'shed_defs/'+'AfMe_traj_release_new.txt',skip_header=5)
    line3 = pl.genfromtxt(sheddir+'shed_defs/'+'EAA_traj_release_new.txt',skip_header=5)
    line4 = pl.genfromtxt(sheddir+'shed_defs/'+'SO_traj_release_new.txt',skip_header=5)
    line5 = pl.genfromtxt(sheddir+'shed_defs/'+'Ar_traj_release_new.txt',skip_header=5)
    
    lw = 1.5
    m.plot(line1[:,0],line1[:,1],color='k',lw=lw,latlon=True)
    m.plot(line2[:,0],line2[:,1],color='k',lw=lw,latlon=True)
    m.plot(line3[:,0],line3[:,1],color='k',lw=lw,latlon=True)
    
    F = pl.zeros_like(line4); F[:,0], F[:,1] = m(line4[:,0],line4[:,1])
    m.plot(F[:246,0], F[:246,1],color='k',lw=lw)
    m.plot(F[246:555,0], F[246:555,1],color='k',lw=lw)
    m.plot(F[555:,0],F[555:,1],color='k',lw=lw)
    
    m.plot(line5[:157,0],line5[:157,1],latlon=True,lw=lw,color='k')
    m.plot(line5[157:,0],line5[157:,1],latlon=True,lw=lw,color='k')
    

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
exec(open('/home/np838619/Trajectory/trajfuncs.py').read())

pl.close('all')

trajdir = '/home/np838619/Trajectory/'
eradir = '/panfs/jasmin/era/era-in/netc/'
sheddir = '/home/np838619/Watershed/'

nc1 = Dataset('/home/np838619/Downloads/etopo05.nc','r')
topo = nc1.variables['ROSE'][:]
toplon = nc1.variables['ETOPO05_X'][:]
toplat = nc1.variables['ETOPO05_Y'][:]
nc1.close()
toplat = pl.flipud(toplat); topo = pl.flipud(topo)


nc2 = Dataset(sheddir+'ggis198901010000.nc','r')
eralon = nc2.variables['longitude'][:]
eralat = nc2.variables['latitude'][:]
geop = nc2.variables['Z'][:]
lsm = nc2.variables['LSM'][:]
nc2.close()
geop = pl.squeeze(geop); lsm = pl.squeeze(lsm)

# calculate orography field by dividing geopotential by gravity:
orog = geop/9.8066

t2 = topo.copy()
t = pl.where(topo<0); t2[t[0],t[1]]=0.


line1 = pl.genfromtxt(sheddir+'shed_defs/'+'Am_traj_release_new.txt',skip_header=5)
line1 = Add360(line1)
line2 = pl.genfromtxt(sheddir+'shed_defs/'+'AfMe_traj_release_new.txt',skip_header=5)
line3 = pl.genfromtxt(sheddir+'shed_defs/'+'EAA_traj_release_new.txt',skip_header=5)
line4 = pl.genfromtxt(sheddir+'shed_defs/'+'SO_traj_release_new.txt',skip_header=5)
line4 = Add360(line4)
line5 = pl.genfromtxt(sheddir+'shed_defs/'+'Ar_traj_release_new.txt',skip_header=5)
line5 = Add360(line5)

l1_top = pl.zeros([line1.shape[0]]); l1_oro = pl.zeros_like(l1_top)
l2_top = pl.zeros([line2.shape[0]]); l2_oro = pl.zeros_like(l2_top)
l3_top = pl.zeros([line3.shape[0]]); l3_oro = pl.zeros_like(l3_top)
l4_top = pl.zeros([line4.shape[0]]); l4_oro = pl.zeros_like(l4_top)
l5_top = pl.zeros([line5.shape[0]]); l5_oro = pl.zeros_like(l5_top)

for pt in range(len(line1)):
    l1_top[pt] = BilinInterp(line1[pt],toplon,toplat,t2)
    l1_oro[pt] = BilinInterp(line1[pt],eralon,eralat,orog)

for pt in range(len(line2)):
    l2_top[pt] = BilinInterp(line2[pt],toplon,toplat,t2)
    l2_oro[pt] = BilinInterp(line2[pt],eralon,eralat,orog)

for pt in range(len(line3)):
    l3_top[pt] = BilinInterp(line3[pt],toplon,toplat,t2)
    l3_oro[pt] = BilinInterp(line3[pt],eralon,eralat,orog)

for pt in range(len(line4)):
    l4_top[pt] = BilinInterp(line4[pt],toplon,toplat,t2)
    l4_oro[pt] = BilinInterp(line4[pt],eralon,eralat,orog)

for pt in range(len(line5)):
    l5_top[pt] = BilinInterp(line5[pt],toplon,toplat,t2)
    l5_oro[pt] = BilinInterp(line5[pt],eralon,eralat,orog)


fig,ax = pl.subplots(5,1,figsize=(10,10))
ax1 = pl.subplot(511)
pl.plot(l1_top,label='ETOP05',lw=2.,color='b')
pl.plot(l1_oro,label='ERA-Interim',lw=2.,color='b',ls='--')
pl.xlim(0,line1.shape[0]); pl.xlabel('latitude',fontsize=16,labelpad=-1)
pl.ylim(0,8000); pl.ylabel('metres',fontsize=16)
a = pl.arange(0,155,20); b = pl.around(line1[a,1],decimals=0)
ax1.set_xticks(a); ax1.set_xticklabels(b.astype('int'),fontsize=13)
pl.yticks(pl.linspace(0,8000,5)); pl.grid(axis='y')
#ax1.text(92,2500,'Panama',fontsize=15)
#pl.arrow(103,2100,0,-1300,head_width=1.5, head_length=350,fc='k')
#ax1.text(80,4200,'Nicaragua',fontsize=15)
#pl.arrow(90,3600,0,-2500,head_width=1.5, head_length=350,fc='k')
#ax1.text(65,3000,'Chivela',fontsize=15)
#pl.arrow(73,2700,1.2,-1800,head_width=1.5, head_length=350,fc='k')
#ax1.text(128,6500,'Andes',fontsize=15)
#pl.annotate(s='', xy=(113,5500), xytext=(155,5500), arrowprops=dict(arrowstyle='<->'))
props = dict(boxstyle='square', facecolor='white', alpha=1.)
ax1.text(5,6000,'(a) Americas',fontsize=18,bbox=props)

ax2 = pl.subplot(512)
pl.plot(l2_top,label='ETOP05',lw=2.,color='k')
pl.plot(l2_oro,label='ERA-Interim',lw=2.,color='k',ls='--')
pl.xlim(0,line2.shape[0]); pl.xlabel('latitude',fontsize=16,labelpad=-1)
pl.ylim(0,8000); pl.ylabel('metres',fontsize=16)
a = pl.arange(0,214,30); b = pl.around(line2[a,1],decimals=0)
ax2.set_xticks(a); ax2.set_xticklabels(b.astype('int'),fontsize=13)
pl.yticks(pl.linspace(0,8000,5)); pl.grid(axis='y')
#ax2.text(76,4010,'Turkana',fontsize=15)
#pl.arrow(86,3700,0,-2300,head_width=1.5, head_length=350,fc='k')
#ax2.text(4,4010,'Al Jaboul',fontsize=15)
#pl.arrow(13,3700,0,-2300,head_width=1.5, head_length=350,fc='k')
ax2.text(11,6000,'(b) Africa',fontsize=18,bbox=props)

ax3 = pl.subplot(513)
pl.plot(l3_top,label='ETOP05',lw=2.,color='g')
pl.plot(l3_oro,label='ERA-Interim',lw=2.,color='g',ls='--')
pl.xlim(0,line3.shape[0]); pl.xlabel('latitude',fontsize=16,labelpad=-1)
pl.ylim(0,8000); pl.legend(loc=0)
a = pl.arange(0,165,20); b = pl.around(line3[a,1],decimals=0)
ax3.set_xticks(a); ax3.set_xticklabels(b.astype('int'),fontsize=13)
pl.yticks(pl.linspace(0,8000,5)); pl.grid(axis='y'); pl.ylabel('metres',fontsize=16)
#ax3.text(145,2010,'Australia',fontsize=15)
#pl.annotate(s='', xy=(146,1300), xytext=(165,1300), arrowprops=dict(arrowstyle='<->'))
#ax3.text(66,2300,'Sumatra',fontsize=15)
#pl.annotate(s='', xy=(67,1500), xytext=(83,1500), arrowprops=dict(arrowstyle='<->'))
ax3.text(11,6000,'(c) South-East Asia',fontsize=18,bbox=props)

z=pl.where(line5[:,0]>180)
line5[z[0],0] = line5[z[0],0] - 360

ax5 = pl.subplot(514)
pl.plot(l5_top,label='ETOP05',lw=2.,color='r')
pl.plot(l5_oro,label='ERA-Interim',lw=2.,color='r',ls='--')
pl.xlim(0,line5.shape[0]); pl.xlabel('longitude',fontsize=16,labelpad=-1)
pl.ylim(0,8000); pl.ylabel('metres',fontsize=16)
pl.yticks(pl.linspace(0,8000,5)); pl.grid(axis='y')
a = pl.arange(0,518,70); b = pl.around(line5[a,0],decimals=0)
ax5.set_xticks(a); ax5.set_xticklabels(b.astype('int'),fontsize=13)
#ax5.text(300,6405,'Himalayas & Tibet',fontsize=15)
#pl.annotate(s='', xy=(320,6100), xytext=(386,6100), arrowprops=dict(arrowstyle='<->'))
#ax5.text(110,4005,'Greenland',fontsize=15)
#pl.annotate(s='', xy=(125,3300), xytext=(158,3300), arrowprops=dict(arrowstyle='<->'))
ax5.text(11,6000,'(d) Arctic',fontsize=18,bbox=props)

z=pl.where(line4[:,0]>180)
line4[z[0],0] = line4[z[0],0] - 360

ax4 = pl.subplot(515)
pl.plot(l4_top,label='ETOP05',lw=2.,color='darkgoldenrod')
pl.plot(l4_oro,label='ERA-Interim',lw=2.,color='darkgoldenrod',ls='--')
pl.xlim(0,line4.shape[0]); pl.xlabel('longitude',fontsize=16)
pl.ylim(0,8000); pl.ylabel('metres',fontsize=16)
a = pl.arange(0,579,70); b = pl.around(line4[a,0],decimals=0)
ax4.set_xticks(a); ax4.set_xticklabels(b.astype('int'),fontsize=13)
pl.yticks(pl.linspace(0,8000,5)); pl.grid(axis='y')
#ax4.text(370,6005,'South America',fontsize=15)
#pl.annotate(s='', xy=(360,5400), xytext=(490,5400), arrowprops=dict(arrowstyle='<->'))
#ax4.text(130,2200,'Australia',fontsize=15)
#pl.annotate(s='', xy=(110,1900), xytext=(214,1900), arrowprops=dict(arrowstyle='<->'))
ax4.text(11,6000,'(e) Southern Ocean',fontsize=18,bbox=props)

pl.tight_layout()
pl.subplots_adjust(hspace=0.35)

#t = pl.where(topo<=0); topo[t[0],t[1]]=pl.float64('nan')
#orog = orog*lsm
#p = pl.where(orog<=0); orog[p[0],p[1]] = pl.float64('nan')
#
#m = Basemap(projection='robin',lon_0=180.)
#levels = [0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000]
#cmap = 'hot_r'
#
#
#fig, ax = pl.subplots(2,1,figsize=(10,8))
#
#ax1 = pl.subplot(211)
#m.drawcoastlines()
#lons, lats = pl.meshgrid(toplon,toplat)
#X1, Y1 = m(lons,lats)
#c1 = m.contourf(X1,Y1,topo,norm=pl.Normalize(0,6000,clip=True),extend='max',levels=levels,
#                cmap=cmap)
#pl.title('ETOPO05',fontsize=16)
#Lines(m)
#
#ax2 = pl.subplot(212)
#m.drawcoastlines()
#lons, lats = pl.meshgrid(eralon,eralat)
#X2, Y2 = m(lons,lats)
#c2 = m.contourf(X2,Y2,orog,norm=pl.Normalize(0,6000,clip=True),extend='max',levels=levels,
#                cmap=cmap)
#pl.title('ERA-Interim',fontsize=16)
#Lines(m)
#
#f = pl.gcf()
#colax = f.add_axes([0.9,0.25,0.015,0.5])
#clb = pl.colorbar(c2,cax=colax,extend='min')
#clb.set_label('metres',fontsize=16)
#pl.tight_layout()