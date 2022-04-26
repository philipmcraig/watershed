# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 14:42:58 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import itertools

pl.close('all')
sheddir = '/home/np838619/Watershed/'
gldir = sheddir + 'GRE_Basins_IMBIE2_v13/'
nadir = sheddir + 'NA_Watersheds/'

m = Basemap(projection='aea',lat_0=71,lon_0=-42,lat_1=65,lat_2=80,width=4000000,
            height=2700000,resolution='l')
m.drawcoastlines(linewidth=0.1)
gl_shp = m.readshapefile(gldir+"GRE_Basins_IMBIE2_v13","gl",linewidth=0.5)
#na_shp = m.readshapefile(nadir+"watershed_p_v2","na",linewidth=0.5)