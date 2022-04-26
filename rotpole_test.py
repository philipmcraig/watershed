# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:52:13 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans

def MapCoords(coords,basinmap):
    """Function to convert from the matplotlib figure co-ordinates to the actual
    longitude-latitude co-ordinates in degrees.
    
    Args:
        coords (list): x,y co-ordinates of clicked points from matplotlib figure
        basinmap (Basemap object): map used for clicking
    
    Returns:
        lonlat (array): longitude-latitude co-ordinates of clicked points
    """
    # need empty array, len(coords) X 2:
    boundary = pl.zeros([len(coords),len(coords[0])]) 
    boundary[:] = coords[:] # make every row of boundary a pair of co-ordinates from coords
    # Transform co-ordinates from figure x,y to lon-lat:
    ilon, ilat = basinmap(boundary[:,0],boundary[:,1],inverse=True)
    lonlat=pl.array([ilon,ilat]).T # make array of lon, lat values
    
    return lonlat

bounds = pl.array([19,38,194,78])
#m = Basemap(llcrnrlon=bounds[0],llcrnrlat=bounds[1],urcrnrlon=bounds[2],
#                urcrnrlat=bounds[3],lon_0=100.,
#                o_lat_p=90.,o_lon_p=90.,resolution='l',projection='rotpole')
m=Basemap(lat_0=65,lon_0=100.,lat_1=70,lat_2=70,width=10000000,height=6000000,
             resolution='l',projection='aea')
m.drawcoastlines()

coords = pl.ginput(n=0,timeout=0)
lonlat = MapCoords(coords,m)
boundary = pl.zeros_like(lonlat)
boundary[:,0], boundary[:,1] = m(lonlat[:,0],lonlat[:,1],inverse=False)
m.plot(boundary[:,0],boundary[:,1],'b') # add line segments
m.plot(boundary[:,0],boundary[:,1],'ro')


#m2 = Basemap(lat_0=65,lon_0=100.,lat_1=70,lat_2=70,width=10000000,height=6000000,
#             resolution='l',projection='aea')
#m2.drawcoastlines()
#as_shp = m2.readshapefile("as_bas_30s_beta","as",linewidth=0.5)