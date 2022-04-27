# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:58:06 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import os
import cartopy
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
panfs = '/panfs/jasmin/era/era-in/netc/monthly_means/'
sheddir = '/home/np838619/Watershed/'
empdir = '/home/np838619/PminusE_data/ERA_Int/'

ncfile = Dataset(empdir+'evapprec_36yrs.nc','r')
#prec = ncfile.variables['prec'][:]
evap = ncfile.variables['evap'][:]
lon = ncfile.variables['lon'][:]
lat = ncfile.variables['lat'][:]
ncfile.close()

# list of years
years = pl.linspace(1979,2014,36)
# need next part to get rid of decimal point
year_input = [str(int(i)) for i in years]

# empty array for filenames
filenames = pl.zeros([len(years),96],dtype='S17') # Only need every 4th file!

#loop over years:
for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
    #filenames[year] = PrintFiles(path,type)
    filenames[year] = PrintFiles(path,'hafs')

lessfiles = pl.zeros([filenames.shape[0],filenames.shape[1]/4],dtype='S17')
for year in range(filenames.shape[0]):
    filelist = []
    for name in range(filenames.shape[1]):
        if '12.nc' in filenames[year,name]:
            filelist.append(filenames[year,name])
    for i in range(len(filelist)):
        lessfiles[year,i] = filelist[i]

#empty array for evaporation
evap = pl.zeros([len(years),lessfiles.shape[1],1,1,256,512])
#empty array for precipitation
prec = pl.zeros_like(evap)

lessfiles = pl.sort(lessfiles)

#loop over years:
for year in range(len(years)):
    #loop over filenames:
    for name in range(lessfiles.shape[1]):
        #load ncfile
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                str(year_input[year]) + '/' + str(lessfiles[year,name]),'r')
        #extract E & TP data
        evap[year,name] = ncfile.variables['E'][:]
        prec[year,name] = ncfile.variables['TP'][:]
        ncfile.close()
 
prec = pl.squeeze(prec); evap = pl.squeeze(evap)       
# why are there -ve values of TP? Nothing to worry about.

evap_tot = pl.zeros([evap.shape[0],evap.shape[1]/2,evap.shape[2],evap.shape[3]])
prec_tot = pl.zeros_like(evap_tot)
# loop over number of years:
for year in range(evap.shape[0]):
    # loop over number of months:
    for month in range(int(evap.shape[1]/2)):
        evap_tot[year,month] = -1*(evap[year][2*month] + evap[year][2*month+1])
        prec_tot[year,month] = (prec[year][2*month] + prec[year][2*month+1])

evap_years_mean = pl.mean(evap_tot,axis=1)
prec_years_mean = pl.mean(prec_tot,axis=1)

# calculate the climatological monthly means:
evap_mnths_mean = pl.mean(evap_tot,axis=0)
prec_mnths_mean = pl.mean(prec_tot,axis=0)

prec_JJA = pl.mean(prec_tot[:,5:8],axis=1)

N = NearestIndex(lat,15)
S = NearestIndex(lat,-15)


#zm = pl.mean(pt,axis=1)

lat_rad = pl.radians(lat[:])
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


#pt_av = pl.sum(pt*areas[N:S+1,:])/pl.sum(areas[N:S+1,:])
fig = pl.figure(figsize=(14,4))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.coastlines(); ax.set_extent([-100.01,20.01,-5.01,15.01],ccrs.PlateCarree())

lat_seg = lat[N:S+1]
itcz_pos = pl.zeros([prec_JJA.shape[0],lon.size])
for yr in range(prec_JJA.shape[0]):
    #for mt in range(prec_tot.shape[1]):
    pt = prec_JJA[yr,N:S+1,:]*1000
    for i in range(lon.size):
        P = pt[:,i]
        zmax = P.max()
        half = zmax/2
        over = pl.where(P>half)
        
        #itcz_pos[yr,i,0] = lat_seg[pl.argmax(P)]
        itcz_pos[yr,i] = pl.sum(P*pl.array(lat_seg))/pl.sum(P)
        #itcz_pos[yr,i,2] = pl.sum(P[over]*pl.array(lat_seg[over]))/pl.sum(P[over])

lon[-1] = 360
ax.plot(lon,itcz_pos[16,:],transform=ccrs.Geodetic(),color='b',lw=3,label='1995')
ax.plot(lon,itcz_pos[26,:],transform=ccrs.Geodetic(),color='r',lw=3,label='2005')
ax.plot(lon,itcz_pos[31,:],transform=ccrs.Geodetic(),color='purple',lw=3,label='2010')
ax.plot(lon,pl.mean(itcz_pos[:,:],axis=0),transform=ccrs.Geodetic(),color='k',
        lw=3,label='climatology')

ax.legend(loc=4,ncol=2,columnspacing=0.4,fontsize=18)
gl = ax.gridlines(ylocs=[-5,0,5,10,15],draw_labels=True)
gl.xlabels_top = False
gl.xformatter = LONGITUDE_FORMATTER; gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 13}; gl.ylabel_style = {'size': 13}
lons,lats=pl.meshgrid(lon,lat)
#ax.contourf(lons,lats,prec_mn,levels=[1,5,10,15,20,25])
pl.subplots_adjust(top=0.9,bottom=0.1,left=0.05,right=0.95)#pl.tight_layout()