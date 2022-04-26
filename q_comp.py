# -*- coding: utf-8 -*-
"""
Created on Sun May 29 13:21:45 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

eradir = '/panfs/jasmin/era/era-in/netc/ggap/'
year1 = '2007'; year2 = '2011'
mon_name = 'jul'#; mon_name2 = 'jul'
mon_no1 = '07'; mon_no2 = '07'

datadir1 = eradir  + year1 + '/' + mon_name + year1 + '/'
datadir2 = eradir  + year2 + '/' + mon_name + year2 + '/'

days = pl.linspace(1,31,31)
days_input = [str(int(i)).zfill(2) for i in days]

times = pl.linspace(0,18,4)
times_input = [str(int(i)).zfill(2) for i in times]

filenames1 = PrintFiles(datadir1,'ggap')
filenames2 = PrintFiles(datadir2,'ggap')

q_jul07 = pl.zeros([len(times)*len(days),1,37,256,512])
q_jul11 = pl.zeros_like(q_jul07)

for name in range(len(filenames1)):
    nc1 = Dataset(datadir1+filenames1[name],'r')
    q_jul07[name] = nc1.variables['Q'][:]
    nc1.close()
    nc2 = Dataset(datadir2+filenames2[name],'r')
    q_jul11[name] = nc2.variables['Q'][:]
    nc2.close()

sheddir = '/home/np838619/Watershed/'
trajq = pl.genfromtxt(sheddir+'spechum.txt'); trajq = pl.reshape(trajq,(16,182))

# remove 1D axes:
q_jul07 = pl.squeeze(q_jul07); q_jul11 = pl.squeeze(q_jul11)
levs = pl.linspace(0,len(q_jul11[0]),len(q_jul11[0]))
pl.plot(q_jul07[84,:,34,274],levs); pl.plot(q_jul11[84,:,34,274],levs)
pl.plot(trajq[:,0],pl.linspace(0,37,16))