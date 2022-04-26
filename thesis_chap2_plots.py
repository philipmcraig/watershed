# -*- coding: utf-8 -*-
"""
Created on Wed May 16 13:49:43 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
sheddir = '/home/np838619/Watershed/'

# read in topography data
ncfile = Dataset('/home/np838619/Downloads/etopo05.nc','r')
topo = ncfile.variables['ROSE'][:] #perhaps extract only values > 0?
lat = ncfile.variables['ETOPO05_Y'][:]
lon = ncfile.variables['ETOPO05_X'][:]
ncfile.close()

#reader = shpreader.Reader(sheddir+'na_bas_30s_beta')

topo, lon2 = shiftgrid(180.0, topo, lon, start=False)
lons,lats = pl.meshgrid(lon2,lat)

#ax = pl.subplot(projection=ccrs.PlateCarree(central_longitude=0))
#ax.coastlines()
#ax.contourf(lons,lats,topo,norm=pl.Normalize(0,8000),levels=pl.linspace(0,8000,9))

pl.figure(figsize=(18,10))
#m = Basemap(projection='cyl',llcrnrlon=-170,llcrnrlat=40,urcrnrlon=-50,urcrnrlat=80)
#m = Basemap(width=2500000,height=3000000,
#            resolution='l',projection='laea',\
#            lat_ts=70,lat_0=70,lon_0=-40.)
m = Basemap(projection='npstere',boundinglat=25,lon_0=0.,round=True)
X, Y = m(lons,lats)
cs=m.contourf(X,Y,topo,norm=pl.Normalize(0,8000),levels=pl.linspace(0,8000,9),
           cmap='YlOrRd',extend='max')
cb = m.colorbar(cs,location='bottom',pad=0.38)
cb.ax.tick_params(labelsize=16)
cb.set_label('m',size=20)


m.drawmeridians([0,90,180,270],linewidth=0.5,dashes=[1,2],
                labels=[1,1,1,1],size=16)
m.drawparallels([20,40,60,80],linewidth=0.5,dashes=[1,2],labels=[1,1,0,0],
                size=16)

lw= 0.3
#na_shp = m.readshapefile(sheddir+'na_bas_30s_beta',"na",linewidth=lw) # North America
#ca_shp = m.readshapefile(sheddir+'ca_bas_30s_beta_ud',"ca",linewidth=lw) # Central America
#sa_shp = m.readshapefile(sheddir+'sa_bas_30s_beta',"sa",linewidth=lw)
#af_shp = m.readshapefile(sheddir+'af_bas_30s_beta',"af",linewidth=lw) # Africa
eu_shp = m.readshapefile(sheddir+'eu_bas_30s_beta',"eu",linewidth=lw)
as_shp = m.readshapefile(sheddir+'as_bas_30s_beta',"as",linewidth=lw) # rest of Asia
#oz_shp = m.readshapefile(sheddir+'Oz_sheds/oz_sheds',"oz",linewidth=0.5)
#au_shp = m.readshapefile(sheddir+'au_bas_30s_beta',"au",linewidth=lw) # Australasia
cd_shp = m.readshapefile(sheddir+'canadaoda_p_1m_v6-0_shp/canadnaoda_p',"cd",linewidth=lw)
gl_shp = m.readshapefile(sheddir+'GRE_Basins_IMBIE2_v13/GRE_Basins_IMBIE2_v13',"gl",linewidth=lw)

amr_pts = pl.genfromtxt(sheddir+'shed_defs/Am_traj_release_new.txt',skip_header=5)
afr_pts = pl.genfromtxt(sheddir+'shed_defs/AfMe_traj_release_new.txt',skip_header=5)
eaa_pts = pl.genfromtxt(sheddir+'shed_defs/EAA_traj_release_new.txt',skip_header=5)
ar_pts = pl.genfromtxt(sheddir+'shed_defs/Ar_traj_release_new.txt',skip_header=5)
#so_pts = pl.genfromtxt(sheddir+'shed_defs/SO_traj_release_new.txt',skip_header=5)
nas_pts = pl.genfromtxt(sheddir+'shed_defs/NAs_traj_release_new.txt',skip_header=5)
#
m.plot(amr_pts[:,0],amr_pts[:,1],lw=0,marker='o',color='g',latlon=True)
m.plot(afr_pts[:,0],afr_pts[:,1],lw=0,marker='o',color='g',latlon=True)
m.plot(eaa_pts[:,0],eaa_pts[:,1],lw=0,marker='o',color='g',latlon=True)
m.plot(ar_pts[:,0],ar_pts[:,1],lw=lw,marker='o',color='g',latlon=True)
#m.plot(so_pts[:,0],so_pts[:,1],lw=0,marker='o',color='g',latlon=True)
#
m.plot(amr_pts[0,0],amr_pts[0,1],marker='o',color='r',latlon=True,ms=10)
m.plot(afr_pts[0,0],afr_pts[0,1],marker='o',color='r',latlon=True,ms=10)
m.plot(eaa_pts[0,0],eaa_pts[0,1],marker='o',color='r',latlon=True,ms=10)
m.plot(nas_pts[:,0],nas_pts[:,1],lw=3,color='b',latlon=True)

m.plot(-113.52,48.57,marker='^',color='b',latlon=True,ms=13)


pl.subplots_adjust(left=0.04,right=0.96,top=0.96,bottom=0.08)#pl.tight_layout()