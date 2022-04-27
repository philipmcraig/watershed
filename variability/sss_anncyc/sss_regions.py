# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 16:16:49 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.path as mplPath
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs

def ModLon(region):
    """
    """
    a = pl.where(region[:,0]<-180)
    region[a[0],0] = region[a[0],0] + 360
    
    return region

def CatchPoly(boundary,m,colour):
    """
    """
    #if line == True: 
    #    lw = 0.1
    #elif line == False:
    #    lw = 0.
    bnd_map = pl.zeros_like(boundary)
    bnd_map[:,0], bnd_map[:,1] = m(boundary[:,0],boundary[:,1])
    
    l = list(bnd_map)
    poly = Polygon(l,fill=False,fc='none',edgecolor=colour,linewidth=1,closed=True)
    
    return poly

def RegionMask(lat,lon,region,m):
    """
    """
    elat = lat.copy(); elon = pl.linspace(-331.25,28.25,lon.size)
    latmax = NearestIndex(lat,region[:,1].max())+1
    latmin = NearestIndex(lat,region[:,1].min())-1
    lonmax = NearestIndex(elon,region[:,0].max())+1
    lonmin = NearestIndex(elon,region[:,0].min())-1

    mask = pl.zeros([elat.size,elon.size])
    
    for i in range(latmin,latmax+1):
        for j in range(lonmin,lonmax+1):
            #if yn == 'y':
            #    elon[j] = elon[j] + 360
            a,b = m(elon[j],elat[i],inverse=False)
            pp = mplPath.Path(list(region))
            X = pp.contains_point((a,b))
            if X == True:
                mask[i,j] = 1

    return mask

def MaskShift(mask,lon):
    """
    """
    elon = pl.linspace(-331.25,28.25,lon.size)
    a = pl.where(elon==lon[0]); a = a[0]
    b = pl.where(lon==elon[-1]); b = b[0]
    
    MS = pl.zeros_like(mask)
    MS[:,:b+1] = mask[:,a:]
    MS[:,b+1:] = mask[:,:a]
    
    return MS

pl.close('all')
exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
#exec(open('/home/np838619/Trajectory/trajfuncs.py').read())

clusdir = '/glusterfs/scenario/users/np838619/'
eccodir = clusdir + 'ECCO/release3/'
yu11dir = clusdir + 'Yu11_regions/'

ncfile = Dataset(eccodir+'SALT.0001.nc','r')
#salt = ncfile.variables['SALT'][:]
lat = ncfile.variables['lat'][:]
lon = ncfile.variables['lon'][:]
ncfile.close()

#sss = salt[:,0]
lat = lat[:,0]; lon = lon[0]
#
#sss = pl.reshape(sss,(sss.shape[0]/12,12,lat.size,lon.size))
#sss_ac = pl.mean(sss,axis=0)

AIZ = pl.genfromtxt(yu11dir+'AIZ_pts.txt')
PIZ = pl.genfromtxt(yu11dir+'PIZ_pts.txt')#; PIZ = ModLon(PIZ)
WTNP = pl.genfromtxt(yu11dir+'WTNP_pts.txt')#; WTNP = ModLon(WTNP)
WTSP = pl.genfromtxt(yu11dir+'WTSP_pts.txt')#; WTSP = ModLon(WTSP)
EIO = pl.genfromtxt(yu11dir+'EIO_pts.txt')#; EIO = ModLon(EIO)
WIO = pl.genfromtxt(yu11dir+'WIO_pts.txt')#; WIO = ModLon(WIO)

pl.figure(figsize=(12,6))
m = Basemap(projection='cyl',lon_0=209.,resolution='l',llcrnrlat=-66,
                                    urcrnrlat=66,llcrnrlon=-331,urcrnrlon=29)
#m2 = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80)

aiz_ply = CatchPoly(AIZ,m,'b'); piz_ply = CatchPoly(PIZ,m,'b')
wtnp_ply = CatchPoly(WTNP,m,'b'); wtsp_ply = CatchPoly(WTSP,m,'b')
eio_ply = CatchPoly(EIO,m,'b'); wio_ply = CatchPoly(WIO,m,'b')

m.drawcoastlines()
pl.gca().add_patch(aiz_ply); pl.gca().add_patch(piz_ply)
pl.gca().add_patch(wtnp_ply); pl.gca().add_patch(wtsp_ply)
pl.gca().add_patch(eio_ply); pl.gca().add_patch(wio_ply)

pl.annotate('AIZ',(-42,18),fontsize=18); pl.annotate('PIZ',(-154,20),fontsize=18)
pl.annotate('WTNP',(-229,17),fontsize=18); pl.annotate('WTSP',(-154,-30),fontsize=18)
pl.annotate('EIO',(-268,-26),fontsize=18); pl.annotate('WIO',(-308,-32),fontsize=18)

m.drawmeridians([0,60,120,180,240,300],linewidth=0,labels=[1,0,0,1],fontsize=16)
m.drawparallels([-60,-30,0,30,60],linewidth=0,labels=[1,0,0,1],fontsize=16)
pl.subplots_adjust(left=0.06,right=0.99)
#pl.savefig(yu11dir+'regions.png')

#aiz_msk = RegionMask(lat,lon,AIZ,m); piz_msk = RegionMask(lat,lon,PIZ,m)
#wtnp_msk = RegionMask(lat,lon,WTNP,m); wtsp_msk = RegionMask(lat,lon,WTSP,m)
#eio_msk = RegionMask(lat,lon,EIO,m); wio_msk = RegionMask(lat,lon,WIO,m)

#x = aiz_msk+piz_msk+wtnp_msk+wtsp_msk+eio_msk+wio_msk

#y = MaskShift(x,lon)

#aiz_msk = MaskShift(aiz_msk,lon); piz_msk = MaskShift(piz_msk,lon)
#wtnp_msk = MaskShift(wtnp_msk,lon); wtsp_msk = MaskShift(wtsp_msk,lon)
#eio_msk = MaskShift(eio_msk,lon); wio_msk = MaskShift(wio_msk,lon)

#write to netcdf
#newnc = Dataset(yu11dir+'yu11_regions_masks.nc','w')
##
#lat_dim = newnc.createDimension('lat',lat.size)
#lon_dim = newnc.createDimension('lon',lon.size)
#lat_in = newnc.createVariable('lat',pl.float64,('lat',))
#lat_in.units = 'degrees_north'
#lat_in.long_name = 'latitude'
#lon_in = newnc.createVariable('lon',pl.float64,('lon',))
#lon_in.units = 'degrees_east'
#lon_in.long_name = 'longitude'
#lat_in = lat
#lon_in = lon
##
#AZ_in = newnc.createVariable('Atlantic ITCZ',pl.float64,('lat','lon'))
#AZ_in.standard_name = 'Atlantic Ocean ITCZ mask'
#AZ_in[:,:] = aiz_msk[:,:]
##
#PZ_in = newnc.createVariable('Pacific ITCZ',pl.float64,('lat','lon'))
#PZ_in.standard_name = 'Pacific Ocean IYCZ mask'
#PZ_in[:,:] = piz_msk[:,:]
##
#NP_in = newnc.createVariable('NWT Pacific',pl.float64,('lat','lon'))
#NP_in.standard_name = 'Western Tropical North Pacific mask'
#NP_in[:,:] = wtnp_msk[:,:]
##
#SP_in = newnc.createVariable('SWT Pacific',pl.float64,('lat','lon'))
#SP_in.standard_name = 'Western Tropical South Pacific'
#SP_in[:,:] = wtsp_msk[:,:]
##
#EI_in = newnc.createVariable('East Indian',pl.float64,('lat','lon'))
#EI_in.standard_name = 'East Indian Ocean mask'
#EI_in[:,:] = eio_msk[:,:]
##
#WI_in = newnc.createVariable('West Indian',pl.float64,('lat','lon'))
#WI_in.standard_name = 'West Indian Ocean mask'
#WI_in[:,:] = wio_msk[:,:]
##
#newnc.close()