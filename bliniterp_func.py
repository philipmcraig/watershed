# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 14:48:54 2017

@author: np838619
"""

relpt = rlspts[0]; flux = zon_ann[0]

p1 = pl.zeros([3]); p2 = pl.zeros([3]); p3 = pl.zeros([3]); p4 = pl.zeros([3])
if relpt[0] > lon[-1]:
    a = -1
else:
    a = NearestIndex(lon,relpt[0]) # nearest longitude index
b = NearestIndex(lat,relpt[1]) # nearest latitude index

if relpt[0] > lon[-1]:
    p1[0] = lon[a]; p3[0] = lon[a]; p2[0] = lon[0]; p4[0] = lon[0]
elif lon[a] < relpt[0]: # nearest lon west of relpt
    p1[0] = lon[a]; p3[0] = lon[a];  p2[0] = lon[a+1]; p4[0] = lon[a+1]
elif lon[a] > relpt[0]: # nearest lon east of relpt
    p2[0] = lon[a]; p4[0] = lon[a]; p1[0] = lon[a-1]; p3[0] = lon[a-1]
    
# does not take 0 meridian into account yet

if lat[b] < relpt[1]: # nearest lat south of relpt
    p1[1] = lat[b]; p2[1] = lat[b]; p3[1] = lat[b-1]; p4[1] = lat[b-1]
elif lat[b] > relpt[1]: # nearest lat north of relpt
    p3[1] = lat[b]; p4[1] = lat[b]; p1[1] = lat[b+1]; p2[1] = lat[b+1]

# values of flux at p1,p2,p3,p4:
nrth_lat = pl.where(lat==p3[1]); sth_lat = pl.where(lat==p1[1])
west_lon = pl.where(lon==p1[0]); east_lon = pl.where(lon==p2[0])
p1[2] = flux[sth_lat[0][0],west_lon[0][0]]
p2[2] = flux[sth_lat[0][0],east_lon[0][0]]
p3[2] = flux[nrth_lat[0][0],west_lon[0][0]]
p4[2] = flux[nrth_lat[0][0],east_lon[0][0]]

if relpt[0] > lon[-1]:
    dx = (360. + lon[0]) - lon[-1]
else:
    dx = p2[0] - p1[0]
dy = p3[1] - p2[1]

if relpt[0] > lon[-1]:
    f1 = (((360+p2[0])-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
    f2 = (((360+p2[0])-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
else:
    f1 = ((p2[0]-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
    f2 = ((p2[0]-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]

F = ((p3[1]-relpt[1])/dy)*f1 + ((relpt[1]-p2[1])/dy)*f2