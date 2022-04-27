# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:54:00 2018

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
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import multiprocessing
import timeit

def NaN2Zero(tcdq):
    """Function to convert all the nans in an array into zeros
    """
    # find where the NaNs are in the array:
    nanarray = pl.isnan(tcdq)
    # change the NaNs to zero:
    tcdq[nanarray] = 0.
    
    return tcdq

def SalinityCalc(salinity,areas,basin_mask,lsm):
    """
    """
    S = NaN2Zero(salinity)
    M = NaN2Zero(basin_mask)
    basin = areas*M*S
    meansal = pl.nansum(basin)/pl.nansum(areas*M*(1-lsm))
    
    return meansal

def CellAreas(lat,lon):
    """
    """
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
    
    return areas

def ExtendArray(newlat,lon,old,which):
    """
    """
    break_pt = newlat.size - old.shape[0]
    new = pl.zeros([newlat.size,lon.size])
    if which == 'z':
        new[:,:] = 0.
    elif which == '1':
        new[:,:] = 1.
    elif which == 'n':
        new[:,:] = pl.float32('nan')
    new[break_pt:,:] = old
    
    return new

def ecco_sss(clusdir,years):
    """
    """
    ###########################################################################
    #---------------------------------ECCO SSS---------------------------------
    ###########################################################################
    eccodir = clusdir + 'ECCO/'
    Y = pl.where(years==1992); Y = Y[0]
    
    eccofile = Dataset(eccodir+'release3/eccor3_sss_basins.nc','r')
    atl_ecc = eccofile.variables['Atlantic'][:]
    ind_ecc = eccofile.variables['Indian'][:]
    pac_ecc = eccofile.variables['Pacific'][:]
    arc_ecc = eccofile.variables['Arctic'][:]
    sou_ecc = eccofile.variables['Southern'][:]
    eccolat = eccofile.variables['lat'][:]
    eccolon = eccofile.variables['lon'][:]
    eccofile.close()
    
    # # change shape of arrays to years x months x lat x lon:
    atl_ecc = pl.reshape(atl_ecc,(atl_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    ind_ecc = pl.reshape(ind_ecc,(ind_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    pac_ecc = pl.reshape(pac_ecc,(pac_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    arc_ecc = pl.reshape(arc_ecc,(arc_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    sou_ecc = pl.reshape(sou_ecc,(sou_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    
    # take mean along months axis:
    atl_ecc = pl.mean(atl_ecc,axis=1); ind_ecc = pl.mean(ind_ecc,axis=1)
    pac_ecc = pl.mean(pac_ecc,axis=1); arc_ecc = pl.mean(arc_ecc,axis=1)
    sou_ecc = pl.mean(sou_ecc,axis=1)
    
    eccomask = Dataset(eccodir+'ecco_basin_masks.nc','r')
    atlmsk_ec = eccomask.variables['Atlantic mask'][:]
    indmsk_ec = eccomask.variables['Indian mask'][:]
    pacmsk_ec = eccomask.variables['Pacific mask'][:]
    arcmsk_ec = eccomask.variables['Arctic mask'][:]
    soumsk_ec = eccomask.variables['Southern mask'][:]
    lsm_ec = eccomask.variables['Full mask'][:]
    eccomask.close()
    

    areas_ecc = CellAreas(pl.flipud(eccolat),eccolon)
    
    ecc_sss_yrs = pl.zeros([5,atl_ecc.shape[0]])
    for i in range(atl_ecc.shape[0]):
        ecc_sss_yrs[0,i] = SalinityCalc(atl_ecc[i],areas_ecc,atlmsk_ec,lsm_ec)
        ecc_sss_yrs[1,i] = SalinityCalc(ind_ecc[i],areas_ecc,indmsk_ec,lsm_ec)
        ecc_sss_yrs[2,i] = SalinityCalc(pac_ecc[i],areas_ecc,pacmsk_ec,lsm_ec)
        ecc_sss_yrs[3,i] = SalinityCalc(arc_ecc[i],areas_ecc,arcmsk_ec,lsm_ec)
        ecc_sss_yrs[4,i] = SalinityCalc(sou_ecc[i],areas_ecc,soumsk_ec,lsm_ec)
    
    trends_ecc = pl.zeros([5,2])
    for i in range(5):
        x = linregress(years[Y:],ecc_sss_yrs[i])
        trends_ecc[i,0] = x[0] # slope
        trends_ecc[i,1] = x[1] # intercept
    
    return ecc_sss_yrs, trends_ecc

def soda_sss(clusdir,years,startyear):
    """
    """
    ###########################################################################
    #---------------------------------SODA SSS---------------------------------
    ###########################################################################
    sodadir = clusdir + 'SODA/'
    Y = pl.where(years==startyear); Y = Y[0][0]
    sodayears = pl.linspace(1980,2014,35)
    S = pl.where(sodayears==startyear); S = S[0][0]
    
    sodafile = Dataset(sodadir+'soda_yearlymeans.nc','r')
    atl_sod = sodafile.variables['Atlantic'][:]
    ind_sod = sodafile.variables['Indian'][:]
    pac_sod = sodafile.variables['Pacific'][:]
    arc_sod = sodafile.variables['Arctic'][:]
    sou_sod = sodafile.variables['Southern'][:]
    sodalat = sodafile.variables['lat'][:]
    sodalon = sodafile.variables['lon'][:]
    sodafile.close()
    
    atl_sod[:,256:273,18:61] = pl.float32('nan')
    atl_sod[:,269:285,32:61] = pl.float32('nan')
    atl_sod[:,230:245,54:84] = pl.float32('nan')

    sl2 = pl.linspace(-sodalat[-1],sodalat[-1],sodalon.size/2)
    sodalat = sl2.copy(); del sl2
    
    newatl = []; newind = []; newpac = []; newarc = []; newsou = []
    for i in range(atl_sod.shape[0]):
        newatl.append(ExtendArray(sodalat,sodalon,atl_sod[i],'n'))
        newind.append(ExtendArray(sodalat,sodalon,ind_sod[i],'n'))
        newpac.append(ExtendArray(sodalat,sodalon,pac_sod[i],'n'))
        newarc.append(ExtendArray(sodalat,sodalon,arc_sod[i],'n'))
        newsou.append(ExtendArray(sodalat,sodalon,sou_sod[i],'n'))

    del atl_sod, ind_sod, pac_sod, arc_sod, sou_sod

    atl_sod = pl.asarray(newatl); del newatl
    ind_sod = pl.asarray(newind); del newind
    pac_sod = pl.asarray(newpac); del newpac
    arc_sod = pl.asarray(newarc); del newarc
    sou_sod = pl.asarray(newsou); del newsou
    
    sodamask = Dataset(sodadir+'soda_basin_masks.nc','r')
    atlmsk_sd = sodamask.variables['Atlantic mask'][:]
    indmsk_sd = sodamask.variables['Indian mask'][:]
    pacmsk_sd = sodamask.variables['Pacific mask'][:]
    arcmsk_sd = sodamask.variables['Arctic mask'][:]
    soumsk_sd = sodamask.variables['Southern mask'][:]
    lsm_sd = sodamask.variables['Full mask'][:]
    sodamask.close()
    
    atlmsk_sd[256:273,18:61] = pl.float32('nan')
    atlmsk_sd[269:285,32:61] = pl.float32('nan')
    atlmsk_sd[230:245,54:84] = pl.float32('nan')
    
    atlmsk_sd = ExtendArray(sodalat,sodalon,atlmsk_sd,'n')
    indmsk_sd = ExtendArray(sodalat,sodalon,indmsk_sd,'n')
    pacmsk_sd = ExtendArray(sodalat,sodalon,pacmsk_sd,'n')
    arcmsk_sd = ExtendArray(sodalat,sodalon,arcmsk_sd,'n')
    soumsk_sd = ExtendArray(sodalat,sodalon,soumsk_sd,'1')
    lsm_sd = ExtendArray(sodalat,sodalon,lsm_sd,'1')
    
    areas_sod = CellAreas(pl.flipud(sodalat),sodalon)
    
    sod_sss_yrs = pl.zeros([5,atl_sod.shape[0]])
    for i in range(atl_sod.shape[0]):
        sod_sss_yrs[0,i] = SalinityCalc(atl_sod[i],areas_sod,atlmsk_sd[:],lsm_sd)
        sod_sss_yrs[1,i] = SalinityCalc(ind_sod[i],areas_sod,indmsk_sd[:],lsm_sd)
        sod_sss_yrs[2,i] = SalinityCalc(pac_sod[i],areas_sod,pacmsk_sd[:],lsm_sd)
        sod_sss_yrs[3,i] = SalinityCalc(arc_sod[i],areas_sod,arcmsk_sd[:],lsm_sd)
        sod_sss_yrs[4,i] = SalinityCalc(sou_sod[i],areas_sod,soumsk_sd[:],lsm_sd)
    
    trends_sod = pl.zeros([5,2])
    for i in range(5):
        x = linregress(years[Y:],sod_sss_yrs[i,S:])
        trends_sod[i,0] = x[0] # slope
        trends_sod[i,1] = x[1] # intercept
    
    return sod_sss_yrs, trends_sod

def en4_sss(clusdir,years,startyear):
    """
    """
    ###############################################################################
    #---------------------------------SODA SSS-------------------------------------
    ###############################################################################
    en4dir = clusdir + 'EN4/'
    Y = pl.where(years==startyear); Y = Y[0][0]
    
    en4file = Dataset(en4dir+'en4_yearlymeans.nc','r')
    atl_en4 = en4file.variables['Atlantic'][:]
    ind_en4 = en4file.variables['Indian'][:]
    pac_en4 = en4file.variables['Pacific'][:]
    arc_en4 = en4file.variables['Arctic'][:]
    sou_en4 = en4file.variables['Southern'][:]
    en4lat = en4file.variables['lat'][:]
    en4lon = en4file.variables['lon'][:]
    en4file.close()
    
    atl_en4[:,137:145,8:29] = pl.float32('nan')
    atl_en4[:,144:150,15:29] = pl.float32('nan')
    atl_en4[:,124:131,26:41] = pl.float32('nan')
    
    el2 = pl.arange(-en4lat[-1],en4lat[-1]+1,1)
    en4lat = el2.copy(); del el2
    
    newatl = []; newind = []; newpac = []; newarc = []; newsou = []
    for i in range(atl_en4.shape[0]):
        newatl.append(ExtendArray(en4lat,en4lon,atl_en4[i],'n'))
        newind.append(ExtendArray(en4lat,en4lon,ind_en4[i],'n'))
        newpac.append(ExtendArray(en4lat,en4lon,pac_en4[i],'n'))
        newarc.append(ExtendArray(en4lat,en4lon,arc_en4[i],'n'))
        newsou.append(ExtendArray(en4lat,en4lon,sou_en4[i],'n'))
    
    del atl_en4, ind_en4, pac_en4, arc_en4, sou_en4

    atl_en4 = pl.asarray(newatl); del newatl
    ind_en4 = pl.asarray(newind); del newind
    pac_en4 = pl.asarray(newpac); del newpac
    arc_en4 = pl.asarray(newarc); del newarc
    sou_en4 = pl.asarray(newsou); del newsou
    
    en4mask = Dataset(en4dir+'en4_basin_masks.nc','r')
    atlmsk_e4 = en4mask.variables['Atlantic mask'][:]
    indmsk_e4 = en4mask.variables['Indian mask'][:]
    pacmsk_e4 = en4mask.variables['Pacific mask'][:]
    arcmsk_e4 = en4mask.variables['Arctic mask'][:]
    soumsk_e4 = en4mask.variables['Southern mask'][:]
    lsm_e4 = en4mask.variables['Full mask'][:]
    en4mask.close()
    
    atlmsk_e4[137:145,8:29] = pl.float32('nan')
    atlmsk_e4[144:150,15:29] = pl.float32('nan')
    atlmsk_e4[124:131,26:41] = pl.float32('nan')
    
    atlmsk_e4 = ExtendArray(en4lat,en4lon,atlmsk_e4,'n')
    indmsk_e4 = ExtendArray(en4lat,en4lon,indmsk_e4,'n')
    pacmsk_e4 = ExtendArray(en4lat,en4lon,pacmsk_e4,'n')
    arcmsk_e4 = ExtendArray(en4lat,en4lon,arcmsk_e4,'n')
    soumsk_e4 = ExtendArray(en4lat,en4lon,soumsk_e4,'1')
    lsm_e4 = ExtendArray(en4lat,en4lon,lsm_e4,'1')
    
    areas_en4 = CellAreas(pl.flipud(en4lat),en4lon)
    
    en4_sss_yrs = pl.zeros([5,atl_en4.shape[0]])
    for i in range(atl_en4.shape[0]):
        en4_sss_yrs[0,i] = SalinityCalc(atl_en4[i],areas_en4,atlmsk_e4[:],lsm_e4)
        en4_sss_yrs[1,i] = SalinityCalc(ind_en4[i],areas_en4,indmsk_e4[:],lsm_e4)
        en4_sss_yrs[2,i] = SalinityCalc(pac_en4[i],areas_en4,pacmsk_e4[:],lsm_e4)
        en4_sss_yrs[3,i] = SalinityCalc(arc_en4[i],areas_en4,arcmsk_e4[:],lsm_e4)
        en4_sss_yrs[4,i] = SalinityCalc(sou_en4[i],areas_en4,soumsk_e4[:],lsm_e4)
    
    trends_en4 = pl.zeros([5,2])
    for i in range(5):
        x = linregress(years[Y:],sou_sss_yrs[i,Y:])
        trends_en4[i,0] = x[0] # slope
        trends_en4[i,1] = x[1] # intercept
    
    return en4_sss_yrs, trends_en4

def Ecco_TrendsMap(ax,clusdir,years,ext):
    """
    """
    eccodir = clusdir + 'ECCO/'
    Y = pl.where(years==1992); Y = Y[0]
    
    eccofile = Dataset(eccodir+'release3/eccor3_sss_basins.nc','r')
    atl_ecc = eccofile.variables['Atlantic'][:]
    ind_ecc = eccofile.variables['Indian'][:]
    pac_ecc = eccofile.variables['Pacific'][:]
    arc_ecc = eccofile.variables['Arctic'][:]
    sou_ecc = eccofile.variables['Southern'][:]
    eccolat = eccofile.variables['lat'][:]
    eccolon = eccofile.variables['lon'][:]
    eccofile.close()
    
    # # change shape of arrays to years x months x lat x lon:
    atl_ecc = pl.reshape(atl_ecc,(atl_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    ind_ecc = pl.reshape(ind_ecc,(ind_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    pac_ecc = pl.reshape(pac_ecc,(pac_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    arc_ecc = pl.reshape(arc_ecc,(arc_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    sou_ecc = pl.reshape(sou_ecc,(sou_ecc.shape[0]/12,12,eccolat.size,eccolon.size))
    
    # take mean along months axis:
    atl_ecc = pl.mean(atl_ecc,axis=1); ind_ecc = pl.mean(ind_ecc,axis=1)
    pac_ecc = pl.mean(pac_ecc,axis=1); arc_ecc = pl.mean(arc_ecc,axis=1)
    sou_ecc = pl.mean(sou_ecc,axis=1)
    
    trend = pl.zeros_like(atl_ecc[0]); sigs = pl.zeros_like(trend)
    
    for i in range(eccolat.size):
        for j in range(eccolon.size):
            r = linregress(years,sou_ecc[:,i,j])
            trend[i,j] = r[0]; sigs[i,j] = r[3]
    z = pl.where(sigs>0.1)
    sigs[z[0],z[1]] = pl.float64('nan')
    
    proj = ccrs.PlateCarree(central_longitude=-30)
    #ax = pl.subplot(projection=proj)
    ax.set_extent(ext,crs=ccrs.PlateCarree()); ax.coastlines()
    lons,lats = pl.meshgrid(eccolon,eccolat)
    norm = pl.Normalize(-0.03,0.03,clip=False)
    levels = [-0.03,-0.01,-0.005,0.005,0.01,0.03]
    cs = ax.contourf(lons,lats,trend,transform=ccrs.PlateCarree(),norm=norm,
                     levels=levels,cmap='seismic',extend='both')
    ax.contourf(lons,lats,sigs,hatches='/',colors='none',transform=ccrs.PlateCarree())
    
    return trend,sigs

def Soda_TrendsMap(axx,clusdir,years,ext):
    """
    """
    sodadir = clusdir + 'SODA/'
#    Y = pl.where(years==1992); Y = Y[0][0]
#    sodayears = pl.linspace(1980,2014,35)
#    S = pl.where(sodayears==startyear); S = S[0][0]
    
    sodafile = Dataset(sodadir+'soda_yearlymeans.nc','r')
    atl_sod = sodafile.variables['Atlantic'][:]
    ind_sod = sodafile.variables['Indian'][:]
    pac_sod = sodafile.variables['Pacific'][:]
    arc_sod = sodafile.variables['Arctic'][:]
    sou_sod = sodafile.variables['Southern'][:]
    sodalat = sodafile.variables['lat'][:]
    sodalon = sodafile.variables['lon'][:]
    sodafile.close()
    
    atl_sod = NaN2Zero(atl_sod); ind_sod = NaN2Zero(ind_sod)
    pac_sod = NaN2Zero(pac_sod); arc_sod = NaN2Zero(arc_sod)
    sou_sod = NaN2Zero(sou_sod)
    all_sss = atl_sod+ind_sod+pac_sod+arc_sod+sou_sod
    
    trend = pl.zeros_like(atl_sod[0]); sigs = pl.zeros_like(trend)
    
    for i in range(sodalat.size):
        for j in range(sodalon.size):
            r = linregress(years[:],all_sss[:,i,j])
            trend[i,j] = r[0]; sigs[i,j] = r[3]
    z = pl.where(sigs>0.1)
    sigs[z[0],z[1]] = pl.float64('nan')
    
    proj = ccrs.PlateCarree(central_longitude=-30)
    #ax = pl.subplot(projection=proj)
    axx.set_extent(ext,crs=ccrs.PlateCarree()); axx.coastlines()
    lons,lats = pl.meshgrid(sodalon,sodalat)
    norm = pl.Normalize(-1,1,clip=False)
    levels = [-1,-0.75,-0.5,-0.3,-0.1,-0.05,0.05,0.1,0.3,0.5,0.75,1]
    cs = axx.contourf(lons,lats,trend*10,transform=ccrs.PlateCarree(),norm=norm,
                     levels=levels,cmap='seismic',extend='min')
    axx.contourf(lons,lats,sigs,hatches='/',colors='none',transform=ccrs.PlateCarree())
    
    return cs#trend,sigs

def En4_TrendsMap(axx,clusdir,years,ext):
    """
    """
    en4dir = clusdir + 'EN4/'
#    Y = pl.where(years==startyear); Y = Y[0][0]
    
    en4file = Dataset(en4dir+'en4_yearlymeans.nc','r')
    atl_en4 = en4file.variables['Atlantic'][:]
    ind_en4 = en4file.variables['Indian'][:]
    pac_en4 = en4file.variables['Pacific'][:]
    arc_en4 = en4file.variables['Arctic'][:]
    sou_en4 = en4file.variables['Southern'][:]
    en4lat = en4file.variables['lat'][:]
    en4lon = en4file.variables['lon'][:]
    en4file.close()
    
    trend = pl.zeros_like(atl_en4[0]); sigs = pl.zeros_like(trend)
    
    for i in range(en4lat.size):
        for j in range(en4lon.size):
            r = linregress(years[13:],sou_en4[13:,i,j])
            trend[i,j] = r[0]; sigs[i,j] = r[3]
    z = pl.where(sigs>0.1)
    sigs[z[0],z[1]] = pl.float64('nan')
    
    proj = ccrs.PlateCarree(central_longitude=-30)
    #ax = pl.subplot(projection=proj)
    axx.set_extent(ext,crs=ccrs.PlateCarree()); axx.coastlines()
    lons,lats = pl.meshgrid(en4lon,en4lat)
    norm = pl.Normalize(-0.03,0.03,clip=False)
    levels = [-0.03,-0.01,-0.005,0.005,0.01,0.03]
    cs = axx.contourf(lons,lats,trend,transform=ccrs.PlateCarree(),norm=norm,
                     levels=levels,cmap='seismic',extend='both')
    axx.contourf(lons,lats,sigs,hatches='/',colors='none',transform=ccrs.PlateCarree())
    
    return cs

#def main():

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
clusdir = '/glusterfs/scenario/users/np838619/'

sodadir = clusdir + 'SODA/'
en4dir = clusdir + 'EN4/'

years = pl.linspace(1979,2014,36)
sodayears = pl.linspace(1980,2014,35)
eccoyears = pl.linspace(1992,2014,23)
    
#    ecc_sss_yrs, trends_ecc = ecco_sss(clusdir,years)
#    sod_sss_yrs, trends_sod = soda_sss(clusdir,years,1992)
#    en4_sss_yrs, trends_en4 = en4_sss(clusdir,years,1992)
#    
#    en4_sss_lng, trln_en4 = en4_sss(clusdir,years,1979)
#    sod_sss_lng, trln_sod = soda_sss(clusdir,years,1980)
#    
#    ylims = [(35.4,35.8),(34.6,35.1),(34.4,34.6),(28,33),(34.2,34.4)]
#    T = ['(a) Atlantic','(b) Indian','(c) Pacific','(d) Arctic','(e) Southern']
#    
#    fig,ax = pl.subplots(5,1,figsize=(8,12))
#    
#    for i in range(5):
#        axx = pl.subplot(5,1,i+1)
#        axx.plot(years,en4_sss_yrs[i,:],color='r',lw=2,label='EN4')
#        axx.plot(sodayears,sod_sss_yrs[i,:],color='b',lw=2,label='SODA')
#        axx.plot(eccoyears,ecc_sss_yrs[i,:],color='g',lw=2,label='ECCO')
#        
#        axx.plot(eccoyears,trends_en4[i,0]*eccoyears+trends_en4[i,1],lw=2,color='r',ls='--')
#        axx.plot(eccoyears,trends_sod[i,0]*eccoyears+trends_sod[i,1],lw=2,color='b',ls='--')
#        axx.plot(eccoyears,trends_ecc[i,0]*eccoyears+trends_ecc[i,1],lw=2,color='g',ls='--')
#        
#        axx.plot(years,trln_en4[i,0]*years+trln_en4[i,1],lw=2,color='r',ls='--')
#        axx.plot(sodayears,trln_sod[i,0]*sodayears+trln_sod[i,1],lw=2,color='b',ls='--')
#        
#        pl.ylim(ylims[i][0],ylims[i][1]); pl.xlim(1979,2014); axx.grid()
#        pl.xticks(pl.linspace(1979,2014,8)); pl.yticks(fontsize=14)
#        if i == 0:
#            pl.yticks(pl.linspace(ylims[0][0],ylims[0][1],5))
#            axx.legend(loc=4,columnspacing=0.4,fontsize=14,ncol=3)
#        if i != 4:
#            axx.xaxis.set_major_formatter(pl.NullFormatter())
#        
#        pl.title(T[i],loc='left',fontsize=18)
#        pl.ylabel('psu',fontsize=18)
#    
#    axx.tick_params(axis='x', which='major', pad=10)
#    pl.xticks(fontsize=16)    
#    
#    pl.tight_layout()
    
start_time = timeit.default_timer()

#ax = pl.subplots()#1,3,figsize=(13,5.5)
proj = ccrs.PlateCarree(central_longitude=0)
ext = [-180,180,-80,80]
#xlocs = [20,50,80,110,140,170]; ylocs = [-40,-20,0,20,40]
 
#    ax1 = pl.subplot(131,projection=proj)
#    cs=Ecco_TrendsMap(ax1,clusdir,eccoyears,ext)
#    theta = pl.linspace(0, 2*pl.pi, 100)
#    center, radius = [0.5, 0.5], 0.5
#    verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#    circle = mpath.Path(verts * radius + center)
#    ax1.set_boundary(circle, transform=ax1.transAxes)
#    g1 = ax1.gridlines()
#    g1.ylabels_right = False; g1.xlabels_top = False
#    g1.xlocator = mticker.FixedLocator(xlocs)
#    g1.ylocator = mticker.FixedLocator(ylocs)
#    g1.xformatter = LONGITUDE_FORMATTER; g1.yformatter = LATITUDE_FORMATTER
#    g1.xlabel_style = {'size': 13}; g1.ylabel_style = {'size': 13}
#    pl.title('(a) ECCOv4')

fig = pl.figure(figsize=(12,6))
ax = fig.add_subplot(1, 1, 1, projection=proj)
cs=Soda_TrendsMap(ax,clusdir,sodayears,ext)
#    theta = pl.linspace(0, 2*pl.pi, 100)
#    center, radius = [0.5, 0.5], 0.5
#    verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#    circle = mpath.Path(verts * radius + center)
#    ax2.set_boundary(circle, transform=ax2.transAxes)
g2 = ax.gridlines(draw_labels=True)
g2.xlabels_bottom = False#; g2.ylabels_left = False; g2.xlabels_top = False
#    g2.xlocator = mticker.FixedLocator(xlocs)
#    g2.ylocator = mticker.FixedLocator(ylocs)
g2.xformatter = LONGITUDE_FORMATTER; g2.yformatter = LATITUDE_FORMATTER
g2.xlabel_style = {'size': 13}; g2.ylabel_style = {'size': 13}#
#pl.title('(b) SODA')

#    ax3 = pl.subplot(133,projection=proj)
#    cs=En4_TrendsMap(ax3,clusdir,years,ext)
#    theta = pl.linspace(0, 2*pl.pi, 100)
#    center, radius = [0.5, 0.5], 0.5
#    verts = pl.vstack([pl.sin(theta), pl.cos(theta)]).T
#    circle = mpath.Path(verts * radius + center)
#    ax3.set_boundary(circle, transform=ax3.transAxes)
#    g3 = ax3.gridlines()
#    g3.ylabels_left = False; g3.xlabels_top = False
#    g3.xlocator = mticker.FixedLocator(xlocs)
#    g3.ylocator = mticker.FixedLocator(ylocs)
#    g3.xformatter = LONGITUDE_FORMATTER; g3.yformatter = LATITUDE_FORMATTER
#    g3.xlabel_style = {'size': 13}; g1.ylabel_style = {'size': 13}
#    pl.title('(c) EN4')

f = pl.gcf()
colax = f.add_axes([0.1,0.11,0.8,0.04])
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_label('$\\times 10$ psu/yr',fontsize=18)
cb.set_ticks([-1.0,-0.75,-0.5,-0.3,-0.1,-0.05,0.05,0.1,0.3,0.5,0.75,1.0])
cb.set_ticklabels([-1.0,-0.75,-0.5,-0.3,-0.1,-0.05,0.05,0.1,0.3,0.5,0.75,1.0])
cb.ax.tick_params(labelsize=16)
pl.subplots_adjust(left=0.05,right=0.95,top=0.93,wspace=0.07,bottom=0.18)

elapsed = timeit.default_timer() - start_time
print elapsed/60
pl.show()

#if __name__ == "__main__":
#    main()