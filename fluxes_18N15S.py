# -*- coding: utf-8 -*-
"""
Created on Thu May 26 17:16:09 2016

@author: np838619
"""

import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset

pl.close('all')
#==============================================================================
sheddir = '/home/np838619/Watershed/'
S16_flx = pl.genfromtxt(sheddir+'S16_fluxes.txt')
TJ_flx = pl.genfromtxt(sheddir+'flux_18N15S.txt')

AD_line = pl.genfromtxt(sheddir+'AD_wshed.txt')
S16_wshed = pl.zeros_like(AD_line)
S16_wshed[:,0] = pl.flipud(AD_line[:,1]); S16_wshed[:,1] = pl.flipud(AD_line[:,0])
relpts1 = pl.genfromtxt(sheddir+'NCA_traj_release_new.txt',skip_header=5)
relpts2 = pl.genfromtxt(sheddir+'SA_traj_release.txt',skip_header=5)
relpts = pl.zeros([len(relpts1)+len(relpts2),2])
relpts[:len(relpts1)] = relpts1[:]; relpts[len(relpts1):] = relpts2[:]

ax, fig = pl.subplots()
ax1 = pl.subplot(1,1,1)
pl.plot(pl.flipud(S16_wshed[54:91:,1]),pl.flipud(S16_flx),
                                    linewidth=2.,label='Singh et al. (2016)')
pl.plot(pl.flipud(relpts[127:174,1]),pl.flipud(TJ_flx),color='r',
                                linewidth=2.,label='Watershed')
pl.axhline(y=0,color='k',ls='--')
pl.ylabel('Sv',fontsize=22); pl.xlabel('Latitude',fontsize=22)
ax1.tick_params(axis='both',labelsize=20)
pl.legend(loc=3,fontsize=20)
pl.subplots_adjust(left=0.19,right=0.96,bottom=0.11,top=0.92)
#==============================================================================

#==============================================================================


ncfile = Dataset(sheddir+'ggap200707211200.nc','r')
lon = ncfile.variables['longitude'][:]
lat = ncfile.variables['latitude'][:]
ncfile.close()


ax,fig = pl.subplots() #pl.figure(2)
ax2 = pl.subplot(1,1,1)
m = Basemap(projection='cyl',resolution='l',llcrnrlat=-25,urcrnrlat=25,\
        llcrnrlon=240.,urcrnrlon=305.,lat_ts=10)
m.drawcoastlines()
lons, lats = pl.meshgrid(lon,lat)
X, Y = m(lons,lats)

a,b = m(S16_wshed[:,0],S16_wshed[:,1])
m.plot(a,b,color='b',label='Singh et al. (2016)',linewidth=2.)
c,d = m(relpts[:,0]+360.,relpts[:,1])
m.plot(c,d,color='r',label='Watershed',linewidth=2.)
pl.legend(loc=3,fontsize=19)
m.drawparallels([18,12,0,-15],labels=[True,True],fontsize=20)
pl.show()