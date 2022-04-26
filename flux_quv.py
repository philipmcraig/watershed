# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 14:43:16 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset
import os



exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

# list of years
years = pl.linspace(1979,2014,36)
# need next part to get rid of decimal point
year_input = [str(int(i)) for i in years]

# empty array for filenames
filenames = pl.zeros([len(years),12],dtype='S13')

for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
    #filenames[year] = PrintFiles(path,type)
    filenames[year] = PrintFiles(path,'ggap')
filenames = pl.sort(filenames,axis=1)

ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/'+
                year_input[0]+'/'+filenames[0,0],'r')
pres = ncfile.variables['p'][:]
lon = ncfile.variables['longitude'][:]
lat = ncfile.variables['latitude'][:]
ncfile.close()
q = pl.zeros([len(years),1,len(pres),len(lat),len(lon)])
u = pl.zeros_like(q); v = pl.zeros_like(q)


#loop over years:
for year in range(len(years)):
    #loop over filenames:
    q_yr = pl.zeros([12,1,37,256,512])
    u_yr = pl.zeros_like(q_yr); v_yr = pl.zeros_like(q_yr)
    for name in range(filenames.shape[1]):
        #load ncfile
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        #extract E & TP data
        q_yr[name] = ncfile.variables['Q'][:]
        u_yr[name] = ncfile.variables['U'][:]
        v_yr[name] = ncfile.variables['V'][:]
        ncfile.close()
    q[year] = pl.mean(q_yr,axis=0)
    u[year] = pl.mean(u_yr,axis=0); v[year] = pl.mean(v_yr,axis=0)

# remove 1D axes:
q = pl.squeeze(q); u = pl.squeeze(u); v = pl.squeeze(v)

# annual means:
q_mean = pl.mean(q,axis=0)
u_mean = pl.mean(u,axis=0); v_mean = pl.mean(v,axis=0)

# vertical integrals:
dp = pl.zeros(len(pres)-1)
for i in range(len(dp)):
    dp[i] = pres[i]-pres[i+1]
q_int = pl.sum(q_mean,axis=0)
u_int = pl.sum(u_mean,axis=0); v_int = pl.sum(v_mean,axis=0)

relpts = pl.genfromtxt('/home/np838619/Watershed/NCA_traj_release.txt',skip_header=5)