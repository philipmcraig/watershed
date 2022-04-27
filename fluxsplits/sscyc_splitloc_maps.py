# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:41:21 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

pl.close('all')
sheddir = '/home/np838619/Watershed/shed_defs/'

amr = pl.genfromtxt(sheddir+'Am_traj_release_new.txt',skip_header=5)
afr = pl.genfromtxt(sheddir+'AfMe_traj_release_new.txt',skip_header=5)
eaa = pl.genfromtxt(sheddir+'EAA_traj_release_new.txt',skip_header=5)
ari = pl.genfromtxt(sheddir+'ArI_traj_release_new.txt',skip_header=5)

amr_inds = pl.array([70,71,41,61])
afr_inds = pl.array([[50,54,63,45],
                    [185,185,188,188]])
eaa_inds = pl.array([72,44,64,23])
ari_inds = pl.array([[42,40,20,33],
                     [62,65,70,67]])

proj = ccrs.PlateCarree()
lat_formatter = LatitudeFormatter()
lon_formatter = LongitudeFormatter(zero_direction_label=True)

fig, ax = pl.subplots(2,2,figsize=(10,10))

ax1 = pl.subplot(221,projection=proj); ax1.set_extent([-130, -60, -20, 60])
ax1.coastlines(); ax1.gridlines()
ax1.set_yticks([-15, 0, 15, 30, 45, 60], crs=ccrs.PlateCarree())
ax1.set_xticks([-120,-105,-90,-75,-60], crs=ccrs.PlateCarree())
ax1.xaxis.set_major_formatter(lon_formatter)
ax1.yaxis.set_major_formatter(lat_formatter)
ax1.plot(amr[:,0],amr[:,1],transform=ccrs.Geodetic(),color='b',lw=2)
for i in range(4):
    ax1.plot(amr[amr_inds[i],0],amr[amr_inds[i],1],marker='_',ms=10,mew=4)

ax2 = pl.subplot(222,projection=proj); ax2.set_extent([-15, 55, -40, 40])
ax2.coastlines(); ax2.gridlines()
ax2.set_yticks([-40,-20, 0, 20, 40], crs=ccrs.PlateCarree())
ax2.set_xticks([-15,0,15,30,45], crs=ccrs.PlateCarree())
ax2.xaxis.set_major_formatter(lon_formatter)
ax2.yaxis.set_major_formatter(lat_formatter)
ax2.plot(afr[:,0],afr[:,1],transform=ccrs.Geodetic(),color='k',lw=2)
for i in range(4):
    ax2.plot(afr[afr_inds[0,i],0],afr[afr_inds[0,i],1],marker='_',ms=10,mew=4)
    ax2.plot(afr[afr_inds[1,i],0],afr[afr_inds[1,i],1],marker='_',ms=10,mew=4)

ax3 = pl.subplot(223,projection=proj); ax3.set_extent([85, 155, -30, 40])
ax3.coastlines(); ax3.gridlines()
ax3.set_yticks([-30,-15, 0, 15, 30,45], crs=ccrs.PlateCarree())
ax3.set_xticks([90,105,120,135,150], crs=ccrs.PlateCarree())
ax3.xaxis.set_major_formatter(lon_formatter)
ax3.yaxis.set_major_formatter(lat_formatter)
ax3.plot(eaa[:,0],eaa[:,1],transform=ccrs.Geodetic(),color='g',lw=2)
for i in range(4):
    ax3.plot(eaa[eaa_inds[i],0],eaa[eaa_inds[i],1],marker='_',ms=10,mew=4)

ax4 = pl.subplot(224,projection=proj); ax4.set_extent([30, 100, 15, 55])
ax4.coastlines(); ax4.gridlines(ylocs=[15,30,45,60],xlocs=[30,45,60,75,90,105])
ax4.set_yticks([15,30,45,60], crs=ccrs.PlateCarree())
ax4.set_xticks([30,45,60,75,90], crs=ccrs.PlateCarree())
ax4.xaxis.set_major_formatter(lon_formatter)
ax4.yaxis.set_major_formatter(lat_formatter)
ax4.plot(ari[:,0],ari[:,1],transform=ccrs.Geodetic(),color='r',lw=2)
for i in range(4):
    ax4.plot(ari[ari_inds[0,i],0],ari[ari_inds[0,i],1],marker='|',ms=10,mew=4)
    ax4.plot(ari[ari_inds[1,i],0],ari[ari_inds[1,i],1],marker='|',ms=10,mew=4)

pl.tight_layout()