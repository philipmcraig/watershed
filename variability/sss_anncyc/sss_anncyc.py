# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:38:25 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs

def SalinityCalc(salinity,areas,basin_mask):
    """
    """
    S = NaN2Zero(salinity)
    M = NaN2Zero(basin_mask)
    basin = areas*M*S
    meansal = pl.nansum(basin)/pl.nansum(areas*M)
    
    return meansal

pl.close('all')
exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
#exec(open('/home/np838619/Trajectory/trajfuncs.py').read())

clusdir = '/glusterfs/scenario/users/np838619/'
eccodir = clusdir + 'ECCO/release3/'
yu11dir = clusdir + 'Yu11_regions/'

ncfile = Dataset(eccodir+'SALT.0001.nc','r')
salt = ncfile.variables['SALT'][:]
lat = ncfile.variables['lat'][:]
lon = ncfile.variables['lon'][:]
ncfile.close()

regfile = Dataset(yu11dir+'yu11_regions_masks.nc')
aiz_msk = regfile.variables['Atlantic ITCZ'][:]
piz_msk = regfile.variables['Pacific ITCZ'][:]
wtnp_msk = regfile.variables['NWT Pacific'][:]
wtsp_msk = regfile.variables['SWT Pacific'][:]
eio_msk = regfile.variables['East Indian'][:]
wio_msk = regfile.variables['West Indian'][:]
regfile.close()

sss = salt[:,0]
lat = lat[:,0]; lat2=pl.flipud(lat); lon = lon[0]
#
sss = pl.reshape(sss,(sss.shape[0]/12,12,lat.size,lon.size))
sss_ac = pl.mean(sss,axis=0); del sss, salt

lat_rad = pl.radians(lat2[:])
lon_rad = pl.radians(lon[:])
lat_half = HalfGrid(lat_rad)
nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6)
# loop over latitude and longitude
for i in range(lat_half.shape[0]-1): # loops over 256
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(lon.shape[0]): # loops over 512
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)

masks = pl.array([aiz_msk,piz_msk,wtnp_msk,wtsp_msk,eio_msk,wio_msk])
saltcyc = pl.zeros([12,6])

for i in range(12):
    for j in range(6):
        saltcyc[i,j] = SalinityCalc(sss_ac[i],areas,masks[j])

lw = 2
fig,ax = pl.subplots()
ax.plot(saltcyc[:,0],label='AIZ',lw=lw)
ax.plot(saltcyc[:,1],label='PIZ',lw=lw)
ax.plot(saltcyc[:,2],label='WTNP',lw=lw)
ax.plot(saltcyc[:,3],label='WTSP',lw=lw)
ax.plot(saltcyc[:,4],label='EIO',lw=lw)
ax.plot(saltcyc[:,5],label='WIO',lw=lw)
pl.xlim(0,11); pl.ylim(33,37)
pl.ylabel('psu',fontsize=18); pl.yticks(fontsize=15)
ax.set_xticklabels(['Jan','Mar','May','Jul','Sep','Nov'])
pl.xticks(fontsize=15); ax.tick_params(axis='x',pad=10)
ax.grid(axis='y')
ax.legend(loc=3,ncol=3,columnspacing=0.4,fontsize=16)
#pl.savefig(yu11dir+'sss_regions_anncyc_ecco.png')