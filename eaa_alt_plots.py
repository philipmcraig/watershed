# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:31:33 2020

@author: qx911590
"""

from __future__ import division
import pylab as pl
import xarray as xr
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

pl.close('all')
sheddir = '/home/users/qx911590/np838619/Watershed/'

ncfile = xr.open_dataset('/home/users/qx911590/np838619/etopo05.nc')
topo = xr.DataArray(ncfile.variables['ROSE'][:])
lat = xr.DataArray(ncfile.variables['ETOPO05_Y'][:])
lon = xr.DataArray(ncfile.variables['ETOPO05_X'][:])
ncfile.close()

lat = lat.values
lon = lon.values
topo = topo.values
#topo_mask = pl.ma.masked_where(topo < 0.,topo)

rod11 = pl.genfromtxt(sheddir+'shed_defs/Rod11full_clicks.txt',skip_header=5)
ls15 = pl.genfromtxt(sheddir+'shed_defs/LS15full_clicks.txt',skip_header=5)


ax = pl.axes(projection=ccrs.PlateCarree(),extent=[90,150,-25,33])
ax.coastlines(resolution='50m',linewidth=0.5,color='grey')

ax.pcolormesh(lon,lat,topo,transform=ccrs.PlateCarree(),cmap='Greys',
              norm=pl.Normalize(0,8000))

ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['water'])
ax.add_feature(ocean_50m)

ax.plot(rod11[:,0],rod11[:,1],marker='.',ms=10,lw=1.5,color='r',
        label='Rodriguez et al. (2011)',transform=ccrs.PlateCarree())
ax.plot(ls15[:,0],ls15[:,1],marker='x',ms=6,mew=2,lw=1,color='k',
        label='Levang & Schmitt (2015)',transform=ccrs.PlateCarree())
    
ax.legend()