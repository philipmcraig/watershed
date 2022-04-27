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
    cs = axx.contourf(lons,lats,field,cmap='seismic_r',norm=norm,levels=levels,
                 transform=trans,extend='both')
    axx.contourf(lons,lats,sigs,hatches='/',colors='none',transform=trans)
    
    return cs

exec(open('/home/users/np838619/PminusE_data/ERA_Int/functions.py').read())

years = pl.linspace(1979,2014,36)

pl.close('all')
sheddir = '/home/users/np838619/Watershed/'
empdir = '/home/users/np838619/PminusE_data/ERA_Int/'
clusdir = '/glusterfs/scenario/users/np838619/'
eccodir = clusdir + 'ECCO/'
sodadir = clusdir+'SODA/'
en4dir = clusdir + 'EN4/'
eradir = clusdir + 'ERA/'

shedfluxes = pl.genfromtxt(sheddir+'variability/shedfluxes_7914_var.csv',
                                                          skip_header=1)
divQ = pl.genfromtxt(sheddir+'variability/divQ_7914_var.csv',skip_header=1)

ncfile = Dataset(empdir+'tcdq_36years_means.nc','r')
eralat = ncfile.variables['lat'][:]
eralon = ncfile.variables['lon'][:]
tcdq_yrs = ncfile.variables['tcdq'][:]
ncfile.close()

#tcdq_yrs = -1*(365*86400)*tcdq_yrs

maskfile = Dataset(eradir+'era_basin_masks.nc','r')
atlmask = maskfile.variables['Atlantic mask'][:]
indmask = maskfile.variables['Indian mask'][:]
pacmask = maskfile.variables['Pacific mask'][:]
arcmask = maskfile.variables['Arctic mask'][:]
soumask = maskfile.variables['Southern mask'][:]
lsm = maskfile.variables['Land-sea mask'][:]
maskfile.close()

lat_rad = pl.radians(eralat[:])
lon_rad = pl.radians(eralon[:])
lat_half = HalfGrid(lat_rad)
nlon = lon_rad.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon
areas = pl.zeros([lat_rad.shape[0],lon_rad.shape[0]])
radius = 6.37*(10**6)
# loop over latitude and longitude
for i in range(lat_half.shape[0]-1): # loops over 256
    latpair = (lat_half[i],lat_half[i+1])
    for j in range(eralon.shape[0]): # loops over 512
        areas[i,j] = AreaFullGaussianGrid(radius,delta_lambda,latpair)

atl_tcdq = tcdq_yrs*atlmask#*(1-lsm); ind_tcdq = tcdq_yrs*indmask#*(1-lsm)
#pac_tcdq = tcdq_yrs*pacmask*(1-lsm); arc_tcdq = tcdq_yrs*arcmask#*(1-lsm)
#sou_tcdq = tcdq_yrs*soumask#*(1-lsm)

#atl_int = pl.zeros([36]); ind_int = pl.zeros_like(atl_int)
#pac_int = pl.zeros_like(atl_int); arc_int = pl.zeros_like(atl_int)
#sou_int = pl.zeros_like(atl_int)
rho = 10**3
#for yr in range(36):
#    atl_int[yr] = NetEmPCalc(areas,atl_tcdq[yr],rho)
#    ind_int[yr] = NetEmPCalc(areas,ind_tcdq[yr],rho)
#    pac_int[yr] = NetEmPCalc(areas,pac_tcdq[yr],rho)
#    arc_int[yr] = NetEmPCalc(areas,arc_tcdq[yr],rho)
#    sou_int[yr] = NetEmPCalc(areas,sou_tcdq[yr],rho)

#print pl.mean(atl_int)
#print pl.mean(ind_int)
#print pl.mean(pac_int)
#print pl.mean(arc_int)
#print pl.mean(sou_int)

#elat2 = pl.zeros([360]); elat2[:330] = eccolat; elat2[330:] = pl.linspace(-75.25,-89.75,30)
#
#AS = pl.zeros([35,elat2.size,eccolon.size]); AM = pl.zeros_like(AS[0])
#AS[:,:30,:] = 0; AS[:,30:,:] = sou_sss
#AM[:30,:] = 0; AM[30:,:] = soumask
#sou_sss = AS; soumask = AM
#del AS, AM

#x = pl.zeros_like(atl_tcdq[0]); y = pl.zeros_like(x)#; z = pl.zeros_like(x)
#for i in range(eralat.size):
#    for j in range(eralon.size):
#        r = linregress(shedfluxes[:,0]*10,atl_tcdq[:,i,j])
#        x[i,j] = r[0]
#        if r[3] <= 0.1:
#            y[i,j] = r[3]
#        else:
#            y[i,j] = pl.float64('nan')

norm = pl.Normalize(-4,4)#pl.Normalize(-0.7,0.7,clip=False)#
levels = [-4,-2,-1,-0.5,-0.25,0.25,0.5,1,2,4]#[-50,-30,-10,-5,5,10,30,50]#[-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7]#
#
eralon[-1] = 360


#lons,lats = pl.meshgrid(eralon,eralat)
#ax = pl.subplot(projection=ccrs.PlateCarree(central_longitude=-30))
#ax.set_extent([-100,40,-40,80],ccrs.PlateCarree())
#theta = pl.linspace(0, 2*pl.pi, 100)
#center, radius = [0.5, 0.5], 0.5
#verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#circle = mpath.Path(verts * radius + center)
#ax.set_boundary(circle, transform=ax.transAxes)
#
#ax.coastlines(); ax.gridlines()
#cs = ax.contourf(lons,lats,x,cmap='seismic',norm=norm,levels=levels,
#                 transform=ccrs.PlateCarree(),extend='both')
#ax.contourf(lons,lats,y,hatches='/',colors='none',transform=ccrs.PlateCarree())
#pl.title('ERA-Interim $p-e$ and Arctic Pacific $P-E$ correlation',
#         fontsize=18)
#f = pl.gcf()
#colax = f.add_axes([0.13,0.1,0.76,0.05])
#cb=pl.colorbar(cs,cax=colax,orientation='horizontal')
#cb.set_label('mm/day/dSv',fontsize=18)
#pl.tight_layout(); pl.subplots_adjust(bottom=0.16,top=0.95)
#pl.savefig(sheddir+'variability/sss_corr_maps/ERAI_tcdq/sou_tcdq_sou_cor .png')

#lat_formatter = LatitudeFormatter()
#lon_formatter = LongitudeFormatter(zero_direction_label=True)

F = pl.array([shedfluxes[:,0],-shedfluxes[:,1],divQ[:,0]])
T = ['(a) Americas $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
     '(b) Africa $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
#    '(c) Pacific sector $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
   '(c) Atlantic $P-E$']
proj = ccrs.PlateCarree(central_longitude=-30); ext = [-100,40,-40,80]
lons,lats = pl.meshgrid(eralon,eralat)
#fig,ax = pl.subplots(2,2,figsize=(12,12))
#for i in range(4):
#fig = pl.figure(figsize=(12,6))
#ax = fig.add_subplot(1, 1, 1, projection=proj)
#x,y = Regression(years,tcdq_yrs,eralat,eralon)
#    axx = pl.subplot(2,2,i+1,projection=proj)
#cs = PlotStats(ax,x,y,lons,lats,proj,ext,levels,norm,False)
#    pl.title(T[i],fontsize=18,loc='center')
#    theta = pl.linspace(0, 2*pl.pi, 100)
#    center, radius = [0.5, 0.5], 0.5
#    verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#    circle = mpath.Path(verts * radius + center)
#    axx.set_boundary(circle, transform=axx.transAxes)
#    axx.gridlines()
#gx = ax.gridlines(draw_labels=True)
#gx.xlabels_bottom = False
#    if i == 0:
#        gx.ylabels_right = False; gx.xlabels_top = False
#    elif i == 1:
#        gx.ylabels_left = False; gx.xlabels_top = False
#    elif i == 2:
#        gx.ylabels_right = False; gx.xlabels_top = False
#    elif i == 3:
#        gx.ylabels_left = False; gx.xlabels_top = False
#    gx.xlocator = mticker.FixedLocator([-100,-80,-60,-40,-20,0,20,40])
#    gx.ylocator = mticker.FixedLocator([-40,-20,0,20,40,60,80])
#gx.xformatter = LONGITUDE_FORMATTER; gx.yformatter = LATITUDE_FORMATTER
#gx.xlabel_style = {'size': 13}; gx.ylabel_style = {'size': 13}

pl.figure(figsize=(10,10))
gs = gridspec.GridSpec(2, 4)
ig = [gs[0,:2],gs[0,2:],gs[1,1:3]]#[0,0,1]; iy = [:2,2:,1:3]
for i in range(3):
    x,y = Regression(F[i]*10,atl_tcdq*86400,eralat,eralon)
    axx = pl.subplot(ig[i],projection=proj)
    cs = PlotStats(axx,x,y,lons,lats,proj,ext,levels,norm,False)
    pl.title(T[i],fontsize=18,loc='left')
    gx = axx.gridlines(draw_labels=True)
    if i == 0:
        gx.ylabels_right = False; gx.xlabels_top = False
    elif i == 1:
        gx.ylabels_left = False; gx.xlabels_top = False
    elif i == 2:
        gx.xlabels_top = False
    gx.xlocator = mticker.FixedLocator([-100,-80,-60,-40,-20,0,20,40])
    gx.ylocator = mticker.FixedLocator([-40,-20,0,20,40,60,80])
    gx.xformatter = LONGITUDE_FORMATTER; gx.yformatter = LATITUDE_FORMATTER
    gx.xlabel_style = {'size': 13}; gx.ylabel_style = {'size': 13}
##
f = pl.gcf()
colax = f.add_axes([0.13,0.07,0.76,0.03])
cb=pl.colorbar(cs,cax=colax,orientation='horizontal')
cb.set_label('mm day$^{-1}$ dSv$^{-1}$',fontsize=18)
cb.ax.tick_params(labelsize=14)
pl.subplots_adjust(top=0.97,bottom=0.15,left=0.1,right=0.92)#,wspace=0.0)
#pl.subplots_adjust(wspace=-0.2)
pl.savefig(sheddir+'variability/sss_reg_maps/ERA_tcdq/atl_tcdq_panel_new.png')
#pl.show()

#print pl.nanmax(x), pl.nanmin(x)
