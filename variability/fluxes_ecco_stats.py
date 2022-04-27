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
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.gridspec as gridspec

def SalinityCalc(salinity,areas,basin_mask,lsm):
    """
    """
    S = NaN2Zero(salinity)
    M = NaN2Zero(basin_mask)
    basin = areas*M*S
    meansal = pl.nansum(basin)/pl.nansum(areas*M*(1-lsm))
    
    return meansal

def Regression(a,b,lat,lon):
    """
    """
    x = pl.zeros_like(b[0]); y = pl.zeros_like(x)
    for i in range(lat.size):
        for j in range(lon.size):
            r = linregress(a,b[:,i,j])
            x[i,j] = r[0]
            if r[3] <= 0.1:
                y[i,j] = r[3]
            else:
                y[i,j] = pl.float64('nan')
    
    return x,y

def PlotStats(axx,field,sigs,lons,lats,proj,ext,levels,norm,polar):
    """
    """
    trans = ccrs.PlateCarree()
    axx.set_extent(ext,trans)
#    if polar == True:
#        theta = pl.linspace(0, 2*pl.pi, 100)
#        center, radius = [0.5, 0.5], 0.5
#        verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#        circle = mpath.Path(verts * radius + center)
#        axx.set_boundary(circle, transform=axx.transAxes)
    axx.coastlines()#; axx.gridlines(draw_labels=True)
    cs = axx.contourf(lons,lats,field,cmap='seismic',norm=norm,levels=levels,
                 transform=trans,extend='both')
    axx.contourf(lons,lats,sigs,hatches='/',colors='none',transform=trans)
    
    return cs

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
sheddir = '/home/np838619/Watershed/'
clusdir = '/glusterfs/scenario/users/np838619/'
eccodir = clusdir + 'ECCO/'

shedfluxes = pl.genfromtxt(sheddir+'variability/shedfluxes_7914_var.csv',
                                                          skip_header=1)
divQ = pl.genfromtxt(sheddir+'variability/divQ_7914_var.csv',skip_header=1)

ncfile = Dataset(eccodir+'release3/eccor3_sss_basins.nc','r')
eccolat = ncfile.variables['lat'][:]
eccolon = ncfile.variables['lon'][:]
atl_sss = ncfile.variables['Atlantic'][:]
ind_sss = ncfile.variables['Indian'][:]
pac_sss = ncfile.variables['Pacific'][:]
arc_sss = ncfile.variables['Arctic'][:]
sou_sss = ncfile.variables['Southern'][:]
ncfile.close()
#elon2 = eccolon + 180.;
eccolat = pl.flipud(eccolat)

maskfile = Dataset(eccodir+'ecco_basin_masks.nc','r')
#eccolat = ncfile.variables['lat'][:]
#eccolon = ncfile.variables['lon'][:]
atlmask = maskfile.variables['Atlantic mask'][:]
indmask = maskfile.variables['Indian mask'][:]
pacmask = maskfile.variables['Pacific mask'][:]
arcmask = maskfile.variables['Arctic mask'][:]
soumask = maskfile.variables['Southern mask'][:]
lsm = maskfile.variables['Full mask'][:]
maskfile.close()

lat_rad = pl.radians(eccolat[:])
lon_rad = pl.radians(eccolon[:])
lat_half = HalfGrid(lat_rad)
nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6)
# loop over latitude and longitude
for i in range(lat_half.shape[0]-1): # loops over 256
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(eccolon.shape[0]): # loops over 512
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)

atl_sss = pl.reshape(atl_sss,(23,12,eccolat.size,eccolon.size))
ind_sss = pl.reshape(ind_sss,(23,12,eccolat.size,eccolon.size))
pac_sss = pl.reshape(pac_sss,(23,12,eccolat.size,eccolon.size))
arc_sss = pl.reshape(arc_sss,(23,12,eccolat.size,eccolon.size))
sou_sss = pl.reshape(sou_sss,(23,12,eccolat.size,eccolon.size))

atl_yrs = pl.nanmean(atl_sss,axis=1); ind_yrs = pl.nanmean(ind_sss,axis=1)
pac_yrs = pl.nanmean(pac_sss,axis=1); arc_yrs = pl.nanmean(arc_sss,axis=1)
sou_yrs = pl.nanmean(sou_sss,axis=1)

atl_int = pl.zeros([23]); ind_int = pl.zeros_like(atl_int)
pac_int = pl.zeros_like(atl_int); arc_int = pl.zeros_like(atl_int)
sou_int = pl.zeros_like(atl_int)
for yr in range(23):
#    atl_int[yr] = SalinityCalc(atl_yrs[yr],areas,atlmask,lsm)
#    ind_int[yr] = SalinityCalc(ind_yrs[yr],areas,indmask,lsm)
#    pac_int[yr] = SalinityCalc(pac_yrs[yr],areas,pacmask,lsm)
    arc_int[yr] = SalinityCalc(arc_yrs[yr],areas,arcmask,lsm)
#    sou_int[yr] = SalinityCalc(sou_yrs[yr],areas,soumask,lsm)

sss_ints = pl.array([atl_int,ind_int,pac_int,sou_int])

#print pl.mean(atl_int)
#print pl.mean(ind_int)
#print pl.mean(pac_int)
#print pl.mean(arc_int)
#print pl.mean(sou_int)

elat2=eccolat#elat2 = pl.zeros([360]); elat2[:330] = eccolat; elat2[330:] = pl.linspace(-75.25,-89.75,30)
#
#AS = pl.zeros([23,elat2.size,eccolon.size]); AM = pl.zeros_like(AS[0])
#AS[:,:30,:] = 0; AS[:,30:,:] = sou_sss
#AM[:30,:] = 0; AM[30:,:] = soumask
#sou_sss = AS; soumask = AM
#del AS, AM

#x = pl.zeros_like(atl_yrs[0]); y = pl.zeros_like(x)#; z = pl.zeros_like(x)
#for i in range(elat2.size):
#    for j in range(eccolon.size):
#        r = pearsonr(sou_yrs[:,i,j],divQ[13:,4])
#        x[i,j] = r[0]
#        if r[1] <= 0.1:
#            y[i,j] = r[1]
#        else:
#            y[i,j] = pl.float64('nan')
#y[330:,:] = pl.float32('nan')

norm = pl.Normalize(-4,4,clip=False)
levels = [-4,-2,-1,-0.5,0.5,1,2,4]#[-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7]
#
eccolon[-1] = 180

#
#lons,lats = pl.meshgrid(eccolon,elat2)
#ax = pl.subplot(projection=ccrs.SouthPolarStereo())
#ax.set_extent([180,-180,-90,-30],ccrs.PlateCarree())
#theta = pl.linspace(0, 2*pl.pi, 100)
#center, radius = [0.5, 0.5], 0.5
#verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#circle = mpath.Path(verts * radius + center)
#ax.set_boundary(circle, transform=ax.transAxes)
##
#ax.coastlines(); ax.gridlines()
#cs = ax.contourf(lons,lats,pl.flipud(x*soumask),cmap='seismic',norm=norm,levels=levels,
#                 transform=ccrs.PlateCarree(),extend='both')
#ax.contourf(lons,lats,pl.flipud(y),hatches='/',colors='none',transform=ccrs.PlateCarree())
#pl.title('ECCO SSS and Southern Ocean $P-E$ correlation',
#         fontsize=18)
#f = pl.gcf()
#colax = f.add_axes([0.13,0.04,0.76,0.05])
#pl.colorbar(cs,cax=colax,orientation='horizontal')
##
#pl.tight_layout(); pl.subplots_adjust(bottom=0.1,top=0.95)
#pl.savefig(sheddir+'variability/sss_corr_maps/ecco_r3/sou_sss_sou_corr.png')

#F = pl.array([shedfluxes[13:,3],shedfluxes[13:,4],shedfluxes[13:,5],divQ[13:,3]])
#T = ['(a) Atlantic sector $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
#     '(b) Indian sector $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
#    '(c) Pacific sector $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
#   '(d) Arctic $P-E$']
#proj = ccrs.NorthPolarStereo(); ext = [-180,180,60,90]
#lons,lats = pl.meshgrid(eccolon,eccolat)
#fig,ax = pl.subplots(2,2,figsize=(12,12))
#for i in range(4):
#    x,y = Regression(F[i],arc_yrs,eccolat,eccolon)
#    axx = pl.subplot(2,2,i+1,projection=proj)
#    cs = PlotStats(axx,x,y,lons,lats,proj,ext,levels,norm,True)
#    pl.title(T[i],fontsize=18,loc='center')
#    theta = pl.linspace(0, 2*pl.pi, 100)
#    center, radius = [0.5, 0.5], 0.5
#    verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#    circle = mpath.Path(verts * radius + center)
#    axx.set_boundary(circle, transform=axx.transAxes)
#    gx = axx.gridlines()
#    if i == 0 or i == 2:
#        gx.ylabels_right = False; gx.xlabels_top = False
#    elif i == 1 or i == 3:
#        gx.ylabels_left = False; gx.xlabels_top = False
#    gx.xlocator = mticker.FixedLocator([-100,-80,-60,-40,-20,0,20,40])
#    gx.ylocator = mticker.FixedLocator([-40,-20,0,20,40,60,80])
#    gx.xformatter = LONGITUDE_FORMATTER; gx.yformatter = LATITUDE_FORMATTER
#    gx.xlabel_style = {'size': 13}; gx.ylabel_style = {'size': 13}

#pl.figure(figsize=(12,8))
#gs = gridspec.GridSpec(2, 4)
#ig = [gs[0,:2],gs[0,2:],gs[1,1:3]]#[0,0,1]; iy = [:2,2:,1:3]
#for i in range(3):
#    x,y = Regression(F[i],pac_yrs,eccolat,eccolon)
#    axx = pl.subplot(ig[i],projection=proj)
#    cs = PlotStats(axx,x,y,lons,lats,proj,ext,levels,norm,False)
#    pl.title(T[i],fontsize=18,loc='left')
#    gx = axx.gridlines(draw_labels=True)
#    if i == 0:
#        gx.ylabels_right = False; gx.xlabels_top = False
#    elif i == 1:
#        gx.ylabels_left = False; gx.xlabels_top = False
#    elif i == 2:
#        gx.xlabels_top = False
#    gx.xlocator = mticker.FixedLocator([-150,-90,-60,90,150,210])
#    gx.ylocator = mticker.FixedLocator([-60,-40,-20,0,20,40,60,80])
#    gx.xformatter = LONGITUDE_FORMATTER; gx.yformatter = LATITUDE_FORMATTER
#    gx.xlabel_style = {'size': 13}; gx.ylabel_style = {'size': 13}
#
#f = pl.gcf()
#colax = f.add_axes([0.13,0.08,0.76,0.03])
#cb=pl.colorbar(cs,cax=colax,orientation='horizontal')
#cb.set_label('psu/Sv',fontsize=18)
#cb.ax.tick_params(labelsize=14)
#pl.subplots_adjust(top=0.95,bottom=0.15,left=0.05,right=0.95,wspace=0.0)
#pl.savefig(sheddir+'variability/sss_reg_maps/ecco_r3/arc_sss_panel.png')