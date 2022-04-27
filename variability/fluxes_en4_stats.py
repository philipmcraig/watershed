# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 20:19:55 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from scipy.stats import pearsonr, linregress
import matplotlib.path as mpath
import cartopy
import cartopy.crs as ccrs

def SalinityCalc(salinity,areas,basin_mask,lsm):
    """
    """
    S = NaN2Zero(salinity)
    M = NaN2Zero(basin_mask)
    basin = areas*M*S
    meansal = pl.nansum(basin)/pl.nansum(areas*M*(1-lsm))
    
    return meansal

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
sheddir = '/home/np838619/Watershed/'
clusdir = '/glusterfs/scenario/users/np838619/'
en4dir = clusdir + 'EN4/'

shedfluxes = pl.genfromtxt(sheddir+'variability/shedfluxes_7914_var.csv',
                                                          skip_header=1)
divQ = pl.genfromtxt(sheddir+'variability/divQ_7914_var.csv',skip_header=1)

ncfile = Dataset(en4dir+'en4_yearlymeans.nc','r')
en4lat = ncfile.variables['lat'][:]
en4lon = ncfile.variables['lon'][:]
atl_sss = ncfile.variables['Atlantic'][:]
ind_sss = ncfile.variables['Indian'][:]
pac_sss = ncfile.variables['Pacific'][:]
arc_sss = ncfile.variables['Arctic'][:]
sou_sss = ncfile.variables['Southern'][:]
ncfile.close()
#elon2 = eccolon + 180.;
#en4lat = pl.flipud(en4lat)

maskfile = Dataset(en4dir+'en4_basin_masks.nc','r')
#eccolat = ncfile.variables['lat'][:]
#eccolon = ncfile.variables['lon'][:]
atlmask = maskfile.variables['Atlantic mask'][:]
indmask = maskfile.variables['Indian mask'][:]
pacmask = maskfile.variables['Pacific mask'][:]
arcmask = maskfile.variables['Arctic mask'][:]
soumask = maskfile.variables['Southern mask'][:]
lsm = maskfile.variables['Full mask'][:]
maskfile.close()

lat_rad = pl.radians(en4lat[:])
lon_rad = pl.radians(en4lon[:])
lat_half = HalfGrid(lat_rad)
nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6)
# loop over latitude and longitude
for i in range(lat_half.shape[0]-1): # loops over 256
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(en4lon.shape[0]): # loops over 512
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)

#atl_sss = pl.reshape(atl_sss,(23,12,eccolat.size,eccolon.size))
#ind_sss = pl.reshape(ind_sss,(23,12,eccolat.size,eccolon.size))
#pac_sss = pl.reshape(pac_sss,(23,12,eccolat.size,eccolon.size))
#arc_sss = pl.reshape(arc_sss,(23,12,eccolat.size,eccolon.size))
#sou_sss = pl.reshape(sou_sss,(23,12,eccolat.size,eccolon.size))

#atl_yrs = pl.nanmean(atl_sss,axis=1); ind_yrs = pl.nanmean(ind_sss,axis=1)
#pac_yrs = pl.nanmean(pac_sss,axis=1); arc_yrs = pl.nanmean(arc_sss,axis=1)
#sou_yrs = pl.nanmean(sou_sss,axis=1)

atl_int = pl.zeros([36]); ind_int = pl.zeros_like(atl_int)
pac_int = pl.zeros_like(atl_int); arc_int = pl.zeros_like(atl_int)
sou_int = pl.zeros_like(atl_int)
for yr in range(36):
    atl_int[yr] = SalinityCalc(atl_sss[yr],areas,atlmask,lsm)
    ind_int[yr] = SalinityCalc(ind_sss[yr],areas,indmask,lsm)
    pac_int[yr] = SalinityCalc(pac_sss[yr],areas,pacmask,lsm)
    arc_int[yr] = SalinityCalc(arc_sss[yr],areas,arcmask,lsm)
    sou_int[yr] = SalinityCalc(sou_sss[yr],areas,soumask,lsm)

sss_ints = pl.array([atl_int,ind_int,pac_int,arc_int,sou_int])

#print pl.mean(atl_int)
#print pl.mean(ind_int)
#print pl.mean(pac_int)
#print pl.mean(arc_int)
#print pl.mean(sou_int)

elat2 = pl.zeros([180]); elat2[7:] = en4lat; elat2[:7] = pl.linspace(-90,-84,7)

AS = pl.zeros([36,elat2.size,en4lon.size]); AM = pl.zeros_like(AS[0])
AS[:,:7,:] = 0; AS[:,7:,:] = sou_sss
AM[:7,:] = 0; AM[7:,:] = soumask
sou_sss = AS; soumask = AM
del AS, AM

x = pl.zeros_like(sou_sss[0]); y = pl.zeros_like(x)#; z = pl.zeros_like(x)
for i in range(elat2.size):
    for j in range(en4lon.size):
        r = pearsonr(sou_sss[13:,i,j],divQ[13:,4])
        x[i,j] = r[0]
        if r[1] <= 0.1:
            y[i,j] = r[1]
        else:
            y[i,j] = pl.float64('nan')
y[173:,:] = pl.float32('nan')

norm = pl.Normalize(-0.7,0.7,clip=False)
levels = [-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7]
#
en4lon[0] = 0


lons,lats = pl.meshgrid(en4lon,elat2)
ax = pl.subplot(projection=ccrs.SouthPolarStereo())
ax.set_extent([180,-180,-90,-30],ccrs.PlateCarree())
theta = pl.linspace(0, 2*pl.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
ax.set_boundary(circle, transform=ax.transAxes)
#
ax.coastlines(); ax.gridlines()
cs = ax.contourf(lons,lats,x*soumask,cmap='seismic',norm=norm,levels=levels,
                 transform=ccrs.PlateCarree(),extend='both')
ax.contourf(lons,lats,y,hatches='/',colors='none',transform=ccrs.PlateCarree())
pl.title('EN4 SSS and Southern Ocean $P-E$ correlation',
         fontsize=18)
f = pl.gcf()
colax = f.add_axes([0.13,0.04,0.76,0.05])
pl.colorbar(cs,cax=colax,orientation='horizontal')
#
pl.tight_layout(); pl.subplots_adjust(bottom=0.1,top=0.95)
#pl.savefig(sheddir+'variability/sss_corr_maps/EN4/short/sou_sss_sou_corr.png')