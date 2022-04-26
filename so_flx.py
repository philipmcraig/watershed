# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 14:36:06 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset
import os

def ReadTxtFile(filename):
	"""Function to read a .txt file and make every line of the file an element of
	a list.

	Args:
		filename (string): path to .txt file
	
	Returns:
		flist (list): contents of file
	"""
	f = open(filename,'r')
	flist = []
	for line in f.readlines():
		flist.append(line.split())
	f.close()
	
	return flist

def NormalVector(endpts):
    """Function to find the normal vectors to the line segments along a watershed.

	Args:
		sheddir (string): path to directory with end points file
		loc (string): stem indicating which watershed is being used

	Returns:
		n (array): unit normal vectors to line segments
    """
    endpts = Add360(endpts); endpts = pl.radians(endpts)
    R = 6.37*10**6 # radius of the Earth
    # convert endpts array to local cartesian co-ordinates:
    #loccar = LocalCartesian(endpts)
    #lc = pl.zeros_like(endpts)
    #lc[:,0] = R*endpts[:,0]*pl.cos(endpts[:,1])
    #lc[:,1] = R*pl.cos(endpts[:,1])
    
    
    nhat = pl.zeros([endpts.shape[0]-1,2])
    for point in range(endpts.shape[0]-1):
        if endpts[point+1,1] == endpts[point,1]:
            nhat[point] = pl.array([0,1])
        elif endpts[point+1,0] == endpts[point,0]:
            nhat[point] = pl.array([1,0])
        else:
            dx = R*pl.cos(endpts[point,1])*(endpts[point+1,0]-endpts[point,0])
            dy = R*(endpts[point+1,1]-endpts[point,1])
            n = pl.array([-dy,dx])
            nhat[point] = n/pl.norm(n) # normalize

    return nhat

def LocalCartesian(coords):
    """Function to convert latitude-longitude co-ordinates into local
	Cartesian co-ordinates.

	Args:
		co-ords (array): Nx2 array of co-ordinates in degrees

	Returns:
		loccar (array): Nx2 array of local Cartesian co-ordinates
    """
    coords = pl.radians(coords) # convert from degrees to radians
    loccar = pl.zeros_like(coords) # creat empty array same size as coords
    R = 6.37*10**6 # radius of the Earth
	# x = R*cos(lat)*lon
    loccar[:,0] = R*pl.sin(coords[:,1])*pl.cos(coords[:,0])#coords[:,0]#
    loccar[:,1] = R*pl.sin(coords[:,1])*pl.sin(coords[:,0]) #coords[:,1] y = R*lat
    
    return loccar

def TrajSegLabel(loc):
    """Function to assign each trajectory release point a label referring to which
	segment of the watershed it is from.

	Args:
		loc (string): stem of continental watershed e.g. NCA = North/Central America
	
	Returns:
             seglab (array): segment labels (integers) of each release point
             rlspts (array): trajectory release points
    """
    sheddir = '/home/np838619/Watershed/shed_defs/'
    
    #endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
    #endpts = pl.asarray(endlist[5:],dtype='float')
    
    #rlslist = ReadTxtFile(sheddir + loc + '_traj_release.txt')
    #rlspts = pl.asarray(rlslist[5:],dtype='float')
    rlspts = pl.genfromtxt(sheddir+loc+'_traj_release.txt',skip_header=5)
    #if loc == 'Af' or loc == 'Eu':
    #    rlspts = pl.flipud(rlspts)
    
    seglab = [1] # first release point has to be on the first segment
    count = 1
    
    for rls in range(1,rlspts.shape[0]):
        if rlspts[rls,0] == rlspts[rls-1,0] and rlspts[rls,1] == rlspts[rls-1,1]:
            count = count + 1
            seglab.append(count)
        else:
            count = count
            seglab.append(count)
    seglab = pl.asarray(seglab)
    
    return seglab, rlspts

def BilinInterp(relpt,lon,lat,flux):
    """
    """
    # First find p1,p2,p3,p4: the points forming a rectangle around relpt
    # Start with empty arrays for p1,p2,p3,p4; p = (x,y,flux):
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
    
    return F

def MidPts(rlspts):
    """
    """
    #rlspts[:,0] = rlspts[:,0] + 360.
    #loccar = LocalCartesian(rlspts)
    #m = Basemap(llcrnrlon=-180,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=80.,
                #lat_ts=20.,resolution='l',projection='mill',suppress_ticks=True)
    #pl.zeros_like(rlspts)
    #loccar[:,0], loccar[:,1] = m(rlspts[:,0],rlspts[:,1])
    
    dl = pl.zeros([rlspts.shape[0]-1])
    for pt in range(rlspts.shape[0]-1):
        dl[pt] = Haversine(rlspts[pt],rlspts[pt+1])
    
    midpt = pl.zeros([rlspts.shape[0]])
    
    for pt in range(1,rlspts.shape[0]-1):
        midpt[pt] = 0.5*(dl[pt]+dl[pt-1])
    midpt[0] = dl[0]; midpt[-1] = dl[-1]
    
    return midpt

def RR2(labels,loc):
    rlspts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_traj_release.txt',skip_header=5)
    repeat = []
    for r in range(1,rlspts.shape[0]):
        if rlspts[r,0] == rlspts[r-1,0] and rlspts[r,1] == rlspts[r-1,1]:
            repeat.append(r)
    
    L2 = []#; RP2 = [] 
    for r in range(rlspts.shape[0]):
        if r not in repeat:
        #RP2.append(interp_pts[r])
            L2.append(labels[r])
    labs = pl.asarray(L2)# rlspts = pl.asarray(RP2)
    
    return labs

def MidFlux(rlspts,lon,lat,zon,mer):
    """
    """
    flux_uv = pl.zeros_like(rlspts); 
    for i in range(len(flux_uv)):
        a = NearestIndex(lon,rlspts[i,0]); b = NearestIndex(lat,rlspts[i,1])
        if lon[a] == rlspts[i,0]:
            flux_uv[i,0] = zon[b,a]; flux_uv[i,1] = mer[b,a]
        else:
            flux_uv[i,0] = BilinInterp(rlspts[i],lon,lat,zon)
            flux_uv[i,1] = BilinInterp(rlspts[i],lon,lat,mer)
    
    midflux = pl.zeros_like(rlspts)
    for f in range(1,midflux.shape[0]-1):
        midflux[f,0] = (flux_uv[f,0]+flux_uv[f-1,0])/2
        midflux[f,1] = (flux_uv[f,1]+flux_uv[f-1,1])/2
    midflux[0,:] = flux_uv[0,:]
    midflux[-1,:] = flux_uv[-1,:]
    
    return midflux

def MF2(rlspts,lon,lat,zon,mer):
    """
    """
    flux_uv = pl.zeros_like(rlspts); 
    for i in range(len(flux_uv)):
        a = NearestIndex(lon,rlspts[i,0]); b = NearestIndex(lat,rlspts[i,1])
        if lon[a] == rlspts[i,0]:
            flux_uv[i,0] = zon[b,a]; flux_uv[i,1] = mer[b,a]
        else:
            flux_uv[i,0] = BilinInterp(rlspts[i],lon,lat,zon_mean)
            flux_uv[i,1] = BilinInterp(rlspts[i],lon,lat,mer_mean)
    
    midflux = pl.zeros([rlspts.shape[0]-1,2])
    for f in range(midflux.shape[0]):
        midflux[f,0] = (flux_uv[f+1,0]+flux_uv[f,0])/2
        midflux[f,1] = (flux_uv[f+1,1]+flux_uv[f,1])/2
    #midflux[0,:] = flux_uv[0,:]
    #midflux[-1,:] = flux_uv[-1,:]
    
    return midflux


def NormalFlux(midflux,labs,nhat):
    """
    """
    FdotN = pl.zeros([midflux.shape[0]])
    for i in range(FdotN.shape[0]):
        segno = labs[i,0] - 1
        FdotN[i] = pl.dot(midflux[i],nhat[segno])
        #if segno not in labs:
        #    pass
        #else:
            #print i, segno
            #FdotN[i] = pl.dot(midflux[i],nhat[segno])
    
    return FdotN

def Add360(rlspts):
    """
    """
    for i in range(rlspts.shape[0]):
        if rlspts[i,0] < 0.:
            rlspts[i,0] = rlspts[i,0] + 360.
    
    return rlspts

def Haversine(pt1,pt2):
    """
    """
    pt1 = pl.radians(pt1); pt2 = pl.radians(pt2)
    
    phi1 = pt1[1]; phi2 = pt2[1]; dphi = phi2 - phi1
    lam1 = pt1[0]; lam2 = pt2[0]; dlam = lam2 - lam1
#    
    dsig = 2*pl.arcsin(pl.sqrt((pl.sin(dphi/2)**2 + pl.cos(phi1)*pl.cos(phi2)*pl.sin(dlam/2)**2)))
    
    R = 6.37e6
    
    d = R*dsig
    
    return d

def ShedFluxes(sheddir,loc,lon,lat,zon,mer):
    """
    """
    endpts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_clicks.txt',skip_header=5)
    labs = TrajSegLabel(loc)
    l2 = pl.zeros([labs[0].size,3]); l2[:,0] = labs[0]; l2[:,1:] = labs[1]
    labs = l2.copy()
    nhat = NormalVector(endpts)
    labs = RR2(labs,loc)
    rlspts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_traj_release_new.txt',skip_header=5)
    rlspts = Add360(rlspts)
    midpt = MidPts(rlspts)
    midflux = MidFlux(rlspts,lon,lat,zon,mer)
    FdotN = NormalFlux(midflux,labs,nhat)
    fluxes = FdotN[:]*midpt[:]/(10**9)
    
    return fluxes

def FluxPerLength(rlspts,fluxes):
    """
    """
    distances = pl.zeros([rlspts.shape[0]-1])
    
    for pt in range(1,rlspts.shape[0]-1):
        distances[pt] = Haversine(rlspts[pt-1],rlspts[pt])
    
    print pl.sum(distances)/1000, 'km'
    flux_per_len = pl.sum(fluxes)/pl.sum(distances)
    
    return flux_per_len

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
panfs = '/panfs/jasmin/era/era-in/netc/monthly_means/'
sheddir = '/home/np838619/Watershed/'

ncfile = Dataset(sheddir+'wvfluxes_7914.nc','r')
lon = ncfile.variables['lon'][:]
lat = ncfile.variables['lat'][:]
zon_ann = ncfile.variables['tcuq'][:]
mer_ann = ncfile.variables['tcvq'][:]
ncfile.close()

zon_mean = pl.mean(zon_ann,axis=0); mer_mean = pl.mean(mer_ann,axis=0)

soa_pts = pl.genfromtxt(sheddir+'shed_defs/SOA_traj_release_new.txt',skip_header=5)
soi_pts = pl.genfromtxt(sheddir+'shed_defs/SOI_traj_release_new.txt',skip_header=5)
sop_pts = pl.genfromtxt(sheddir+'shed_defs/SOP_traj_release_new.txt',skip_header=5)

soa_flx = ShedFluxes(sheddir,'SOA',lon,lat,zon_mean,mer_mean)
soi_flx = ShedFluxes(sheddir,'SOI',lon,lat,zon_mean,mer_mean)
sop_flx = ShedFluxes(sheddir,'SOP',lon,lat,zon_mean,mer_mean)

soa_av = FluxPerLength(soa_pts[95:],soa_flx[95:])#; print soa_av*(10**8)
soi_av = FluxPerLength(soi_pts[:113],soi_flx[:113])#; print soi_av*(10**8)
sop_av = FluxPerLength(sop_pts[30:181],sop_flx[30:181])#; print sop_av*(10**8)

print pl.sum(soa_flx[95:]), 'Sv', soa_av, 'Sv/m'
print pl.sum(soi_flx[:113]), 'Sv', soi_av, 'Sv/m'
print pl.sum(sop_flx[30:181]), 'Sv', sop_av, 'Sv/m'