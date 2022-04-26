# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 10:03:52 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset
import os
import glob
import xarray as xr

def NearestIndex(array_in,point_in):
    """Function to the the nearest index to a specified geographic co-ordinate from an
    array of latitude or longitude co-ordinates
    
    Args:
        array_in (array): longitude or latitude array
        point_in (float): longitude or latitude co-ordinate
    
    Returns:
        index (int): index of array which has value closest to point_in
    """
    index = pl.absolute(array_in-point_in).argmin()
    
    return index

def NormalVector(endpts):
    """Function to find the normal vectors to the line segments along a watershed.

    Args:
	endpts (array): longitudes, latitudes of endpoints of catchment boundary segments

    Returns:
	nhat (array): unit normal vectors to line segments
    """
    # add 360 degrees to any longitudes less than zero & convert to radians
    endpts = Add360(endpts); endpts = pl.radians(endpts)
    R = 6.37*10**6 # radius of the Earth

    # convert endpts array to local cartesian co-ordinates:
    #loccar = LocalCartesian(endpts)
    #lc = pl.zeros_like(endpts)
    #lc[:,0] = R*endpts[:,0]*pl.cos(endpts[:,1])
    #lc[:,1] = R*pl.cos(endpts[:,1])
    
    
    nhat = pl.zeros([endpts.shape[0]-1,2]) # one less vector than endpoints, 2 elements
    for point in range(endpts.shape[0]-1):
        if endpts[point+1,1] == endpts[point,1]: # if latitudes of consecutive endpoints are equal...
            nhat[point] = pl.array([0,1]) # unit normal vector pointing north
        elif endpts[point+1,0] == endpts[point,0]: # if longitudes of consecutive endpoints are equal ...
            nhat[point] = pl.array([1,0]) # unit normal vector pointing east
        else:
            dx = R*pl.cos(endpts[point,1])*(endpts[point+1,0]-endpts[point,0])
            dy = R*(endpts[point+1,1]-endpts[point,1])
            n = pl.array([-dy,dx])
            nhat[point] = n/pl.norm(n) # normalize

    return nhat

def LocalCartesian(coords):
    """Function to convert latitude-longitude co-ordinates into local
	Cartesian co-ordinates.
	This function is not used. IGNORE

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
		loc (string): stem of continental watershed e.g. Am = Americas
	
	Returns:
             seglab (array): segment labels (integers) of each release point
             rlspts (array): trajectory release points
    """
    sheddir = '/home/users/qx911590/np838619/Watershed/shed_defs/'
    
    #endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
    #endpts = pl.asarray(endlist[5:],dtype='float')
    
    #rlslist = ReadTxtFile(sheddir + loc + '_traj_release.txt')
    #rlspts = pl.asarray(rlslist[5:],dtype='float')

    # !!!!!!!!must use '_traj_release' files with repeated points!!!!!!!!!
    rlspts = pl.genfromtxt(sheddir+loc+'_traj_release.txt',skip_header=5)
    
    seglab = [1] # first release point has to be on the first segment
    count = 1 # segment number to label points with

    # note that rls in loop below & rlspts refers to "release points"
    # I was thinking of trajectory release points at the time
    # not important for this code
    
    for rls in range(1,rlspts.shape[0]): # start from 1 as 1st label already assigned
        if rlspts[rls,0] == rlspts[rls-1,0] and rlspts[rls,1] == rlspts[rls-1,1]:
	    # if current point is same as previous point ...
            count = count + 1 # ... increase count by 1
            seglab.append(count) # assign this point with new label
        else: # current point is not the same as previous point ...
            count = count # ... do not increase count
            seglab.append(count) # ... assign this point with same label
    seglab = pl.asarray(seglab) # easier to work with an array
    
    return seglab, rlspts

def BilinInterp(relpt,lon,lat,flux):
    """Bilinear Interpolation of a gridded variable to a point along a catchment boundary.
	This function is written specifically to work with the ERA-Interim longitude & 
	latitude arrays. It may require editing to work for other datasets.
	https://en.wikipedia.org/wiki/Bilinear_interpolation

    Args:
	relpt (array): longitude, latitude co-ordinates of point on boundary
	lon (array): longitude array of gridded dataset
	lat (array): latitude array of gridded dataset
	flux (array): 2D array of variable to be interpolated

    Returns:
	F (float): variable interpolated from grid to point on boundary
    """
    # First find p1,p2,p3,p4: the points forming a rectangle around relpt
    # Start with empty arrays for p1,p2,p3,p4; p = (x,y,flux):
    p1 = pl.zeros([3]); p2 = pl.zeros([3]); p3 = pl.zeros([3]); p4 = pl.zeros([3])
    if relpt[0] > lon[-1]: # if longitude is greater than final element on lon array ...
        a = -1 # nearest index is -1 (final element of array in Python)
    else: # otherwise find the nearest index in lon array
        a = NearestIndex(lon,relpt[0]) # nearest longitude index
    b = NearestIndex(lat,relpt[1]) # nearest latitude index
    
    if relpt[0] > lon[-1]: # if longitude is greater than final lon array element
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
    
    if relpt[0] > lon[-1]: # if longitude is greater then final lon array element
        dx = (360. + lon[0]) - lon[-1]
    else:
        dx = p2[0] - p1[0]
    dy = p3[1] - p2[1]
    
    if relpt[0] > lon[-1]: # if longitude is greater then final lon array element
        f1 = (((360+p2[0])-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
        f2 = (((360+p2[0])-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
    else:
        f1 = ((p2[0]-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
        f2 = ((p2[0]-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
    
    F = ((p3[1]-relpt[1])/dy)*f1 + ((relpt[1]-p2[1])/dy)*f2
    
    return F

def MidPts(rlspts):
    """Calculate the distance between the midpoints between points along catchment boundary.
	Use this in the main code to weight the moisture flux normal to the boundary at each
	midpoint - this gives the total flux between two points along the boundary.

    Args:
	rlspts (array): longitude, latitude co-ordinates of points along catchment boundary

    Returns:
	midpt (array): distances between midpoints
    """
    #rlspts[:,0] = rlspts[:,0] + 360.
    #loccar = LocalCartesian(rlspts)
    #m = Basemap(llcrnrlon=-180,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=80.,
                #lat_ts=20.,resolution='l',projection='mill',suppress_ticks=True)
    #pl.zeros_like(rlspts)
    #loccar[:,0], loccar[:,1] = m(rlspts[:,0],rlspts[:,1])
    
    # calculate the great circle distance between points along catchment boundary
    dl = pl.zeros([rlspts.shape[0]-1]) # one less distance between points than number of points
    for pt in range(rlspts.shape[0]-1):
        dl[pt] = Haversine(rlspts[pt],rlspts[pt+1]) # great circle distance between points
    
    midpt = pl.zeros([rlspts.shape[0]]) # empty array for distances between midpoints
    
    for pt in range(1,rlspts.shape[0]-1):
        midpt[pt] = 0.5*(dl[pt]+dl[pt-1]) # average of 2 adjacent great circle distances
    # can't take average of 2 dl for endpoints, so just use first & last dl
    midpt[0] = dl[0]; midpt[-1] = dl[-1]
    
    return midpt

def RR2(interp_pts,labels):
    """Find and remove repeated points along catchment boundaries. This avoids duplicate
	calculations and zero lengths between repeated points.
	RR stands for "Remove Repeated". This is the 2nd version of the function, no idea
	where the first one is.

    Args:
	interp_pts (array): longitude, latitude points of points along catchment boundary
	labels (array): segment labels

    Returns:
	rlspts (array): longitude, latitude points of non-repeated points along boundary
	labs (array): segment labels
    """
    repeat = [] # empty list for repeated points
    for r in range(1,interp_pts.shape[0]): # start from index 1
	# check if longitude point equals previous & if latitude point equals previous
        if interp_pts[r,0] == interp_pts[r-1,0] and interp_pts[r,1] == interp_pts[r-1,1]:
            repeat.append(r) # put any repeated points in the list
    
    RP2 = []; L2 = [] # empty lists for non-repeated boundary points and segment labels
    for r in range(interp_pts.shape[0]):
        if r not in repeat: # if this point is not repeated...
            RP2.append(interp_pts[r]) # add the point to list
            L2.append(labels[r]) # and add the segment label to list
    rlspts = pl.asarray(RP2); labs = pl.asarray(L2) # turn the lists into arrays
    
    return rlspts, labs

def MidFlux(rlspts,lon,lat,zon,mer):
    """Calculate the moisture flux at the midpoints between points along catchment boundary

    Args:
	rlspts (array): longitude, latitude co-ordinates of points along catchment boundary
	lon (array): longitude array of gridded dataset
	lat (array): latitude array of gridded dataset
	zon (array): zonal moisture flux (qu)
	mer (array): meridional moisture flux (qv)

    Returns:
	midflux (array): zonal & meridional moisture fluxes interpolated to the midpoints
    """
    flux_uv = pl.zeros_like(rlspts); 
    for i in range(len(flux_uv)):
	# find nearest indices in lon,lat arrays to the co-ordinates of each point along boundary
        a = NearestIndex(lon,rlspts[i,0]); b = NearestIndex(lat,rlspts[i,1])
        if lon[a] == rlspts[i,0]: # check if nearest longitude is the same as the longitude of boundary
	    # if it is, then just use the values of qu & qv at that point (no interpolation needed)
            flux_uv[i,0] = zon[b,a]; flux_uv[i,1] = mer[b,a]
        else: # otherwise use bilinear interpolation
            flux_uv[i,0] = BilinInterp(rlspts[i],lon,lat,zon)
            flux_uv[i,1] = BilinInterp(rlspts[i],lon,lat,mer)
    
    midflux = pl.zeros_like(rlspts)
    for f in range(1,midflux.shape[0]-1):
	# flux at midpoints is average of flux at 2 adjacent points along boundary
        midflux[f,0] = (flux_uv[f,0]+flux_uv[f-1,0])/2
        midflux[f,1] = (flux_uv[f,1]+flux_uv[f-1,1])/2
    # flux at end points can't be average of 2 points so keep them the same
    midflux[0,:] = flux_uv[0,:]
    midflux[-1,:] = flux_uv[-1,:]
    
    return midflux

def MF2(rlspts,lon,lat,zon,mer):
    """Don't seem to be using this function. IGNORE
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
    """Calculate the moisture flux normal to a catchment boundary at each
	point along the boundary.

    Args:
	midflux (array): qu & qv at the midpoints between points along boundaries
	labs (array): segment labels
	nhat (array): normal vectors for each segment of the catchment boundary
 
    Returns:
	FdotN (array): moisture flux normal to each midpoint along boundary
    """
    FdotN = pl.zeros([midflux.shape[0]])
    for i in range(FdotN.shape[0]):
        segno = labs[i] - 1 # required because Python is zero based
        FdotN[i] = pl.dot(midflux[i],nhat[segno]) # scalar product
        #if segno not in labs:
        #    pass
        #else:
            #print i, segno
            #FdotN[i] = pl.dot(midflux[i],nhat[segno])
    
    return FdotN

def Add360(rlspts):
    """Add 360 degrees to negative longitude values. Required because some catchment boundaries
	were defined with longitude ranging from -180 to 180.

    Args:
	rlspts (array): longitude, latitude co-ordinates of points along boundary

    Returns:
	rlspts (array): longitude, latitude co-ordinates of points along boundary
    """
    for i in range(rlspts.shape[0]):
        if rlspts[i,0] < 0.: # check if longitude is less than zero
            rlspts[i,0] = rlspts[i,0] + 360. # if less than zero, add 360
    
    return rlspts

def Haversine(pt1,pt2):
    """Haversine formula for calculating the great circle distance between 2 points.
	https://en.wikipedia.org/wiki/Haversine_formula

    Args:
	pt1, pt2 (arrays): longitude, latitude co-ordinates of two points
 
    Returns:
	d (float): distance between pt1 & pt2
    """
    pt1 = pl.radians(pt1); pt2 = pl.radians(pt2) # need to be in radians for this
    # phi1,phi2 are the latitudes & lam1,lam2 are the longitudes
    phi1 = pt1[1]; phi2 = pt2[1]; dphi = phi2 - phi1 # difference in longitude
    lam1 = pt1[0]; lam2 = pt2[0]; dlam = lam2 - lam1 # difference in latitude
    
    # can't remember why the next line is dsig, but the source I initially used for the formula
    # must have used a 'd sigma'
    dsig = 2*pl.arcsin(pl.sqrt((pl.sin(dphi/2)**2 + pl.cos(phi1)*pl.cos(phi2)*pl.sin(dlam/2)**2)))
    
    R = 6.37e6 # radius of Earth
    
    d = R*dsig # great circle distance
    
    return d

def FluxPerLength(rlspts,fluxes):
    """
    """
    distances = pl.zeros([rlspts.shape[0]-1])
    
    for pt in range(1,rlspts.shape[0]-1):
        distances[pt] = Haversine(rlspts[pt-1],rlspts[pt])
    
    print pl.sum(distances)
    flux_per_len = pl.sum(fluxes)/pl.sum(distances)
    
    return flux_per_len

#exec(open('/home/users/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
panfs = '/panfs/jasmin/era/era-in/netc/monthly_means/'
sheddir = '/home/users/qx911590/np838619/Watershed/'

# list of years
years = pl.linspace(2010,2014,5)
# need next part to get rid of decimal point
year_input = [str(int(i)) for i in years]

# empty array for filenames
#filenames = pl.zeros([len(years),12],dtype='S13')
filenames = glob.glob(sheddir+'eraint_vertints_mm_20*')

#for year in range(len(years)):
#    #path = path to mm folder + year
#    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + str(year_input[year])
#    #filenames[year] = PrintFiles(path,type)
#    filenames[year] = PrintFiles(path,'ggaw')
#filenames = pl.sort(filenames,axis=1)

#path = panfs + '2007/'

#ncfile = Dataset(sheddir+'wvfluxes_7914.nc','r')
#lon = ncfile.variables['lon'][:]
#lat = ncfile.variables['lat'][:]
#zon_ann = ncfile.variables['tcuq'][:]
#mer_ann = ncfile.variables['tcvq'][:]
#ncfile.close()

# empty arrays for zonal & meridional water vapour fluxes:
#wv_zon = pl.zeros([filenames.shape[0],filenames.shape[1],1,1,256,512]) # zonal component
#wv_mer = pl.zeros_like(wv_zon) # meridional component

tcwv = pl.zeros([len(filenames),12,241,480])
tcuq = tcwv.copy()
tcvq = tcwv.copy()
tcdq = tcwv.copy()

##loop over years:
#for year in range(len(years)):
#loop over filenames:
for name in range(len(filenames)):
    #load ncfile
    ncfile = xr.open_dataset(filenames[name])
    if name == 0:
        lon = xr.DataArray(ncfile.longitude)
        lat = xr.DataArray(ncfile.latitude)
    tcwv[name,:,:,:] = xr.DataArray(ncfile.tcwv)
    tcuq[name,:,:,:] = xr.DataArray(getattr(ncfile,'p71.162')).data
    tcvq[name,:,:,:] = xr.DataArray(getattr(ncfile,'p72.162')).data
    tcdq[name,:,:,:] = xr.DataArray(getattr(ncfile,'p84.162')).data
    ncfile.close()
#
## remove the 1D axes:
#wv_zon = pl.squeeze(wv_zon,); wv_mer = pl.squeeze(wv_mer)

# climatological monthly means:
#zon_mnths = pl.mean(wv_zon,axis=0); mer_mnths = pl.mean(wv_mer,axis=0)

# annual means:
zon_mean = pl.mean(tcuq,axis=(0,1)); mer_mean = pl.mean(tcvq,axis=(0,1))
#zon_mean = wv_zon[28,6]; mer_mean = wv_mer[28,6]
#zon97 = wv_zon[18,5:]; zon98 = wv_zon[19,:5]
#mer98 = wv_mer[18,5:]; mer98 = wv_mer[19,:5]
#zonEN = pl.zeros([12,256,512]); merEN = pl.zeros_like(zonEN)
#zonEN[:7] = zon88[:]; zonEN[7:] = zon89[:]
#merEN[:7] = mer88[:]; merEN[7:] = mer89[:]
#zon_mean = pl.mean(zonEN,axis=0); mer_mean = pl.mean(merEN,axis=0)
#zon_ann = pl.mean(wv_zon,axis=1); mer_ann = pl.mean(wv_mer,axis=1)
#zon_mean = wv_zon[28,6]; mer_mean = wv_mer[28,6]
#W = pl.sqrt(wv_zon**2 + wv_mer**2) # magnitude of water vapour flux

# make quiver plot:

#latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
#                                    year_input[0] + '/' + filenames[0,0],'r')
#lon = latlon.variables['longitude'][:]
#lat = latlon.variables['latitude'][:]
#latlon.close()

# can't remember why these are here, need to check
#lon_0 = lon.mean()
#lat_0 = lat.mean()

#zon_mn, newlons = shiftgrid(180.0, wv_zon, lon, start=False)
#mer_mn, newlons = shiftgrid(180.0, wv_mer, lon, start=False)
#Ws, newlons = shiftgrid(180.0, W, lon, start=False)

"""pl.figure(1)
m = Basemap(projection='cyl',resolution='l',llcrnrlat=-80,urcrnrlat=80,\
        llcrnrlon=-180.,urcrnrlon=179.,lat_ts=10)
m.drawcoastlines()
lons, lats = pl.meshgrid(newlons,lat)
X, Y = m(lons,lats)
cmap = pl.get_cmap('GnBu')
cs = m.pcolormesh(X, Y, Ws, cmap=cmap)#,norm=pl.Normalize(0,0.015))
uproj,vproj,xx,yy = m.transform_vector(pl.flipud(wv_zon),pl.flipud(wv_mer),
                                       newlons,lat[::-1],nx=30,ny=30,returnxy=True)
m.quiver(xx,yy,uproj,vproj,scale=6000)

pl.figure(2)
m2 = Basemap(llcrnrlon=-180.0,llcrnrlat=0.,urcrnrlon=-45.,urcrnrlat=70.,\
                lat_ts=20.,resolution='l',projection='merc',suppress_ticks=True)
m2.drawcoastlines()
m2.drawmeridians(lon)
m2.drawparallels(lat)"""

# ----------------------------------------------------------
# !!!! the endpts arrays come from the '_clicks' files !!!!!
# ----------------------------------------------------------

endpts1 = pl.genfromtxt(sheddir+'shed_defs/'+'Am_clicks.txt',skip_header=5)
endpts2 = pl.genfromtxt(sheddir+'shed_defs/'+'AfMe_clicks.txt',skip_header=5)
endpts3 = pl.genfromtxt(sheddir+'shed_defs/'+'EAA_clicks.txt',skip_header=5)
endpts4 = pl.genfromtxt(sheddir+'shed_defs/'+'SO_clicks.txt',skip_header=5)
endpts5 = pl.genfromtxt(sheddir+'shed_defs/'+'Ar_clicks.txt',skip_header=5)
endpts6 = pl.genfromtxt(sheddir+'shed_defs/'+'LS15full_clicks.txt',skip_header=5)
#endpts = pl.zeros([len(endpts1)+len(endpts2),2])
#endpts[:len(endpts1)] = endpts1[:]; endpts[len(endpts1):] = endpts2[:]
#AD_line = pl.genfromtxt(sheddir + 'AD_wshed.txt')
#S16_wshed = pl.zeros_like(AD_line)
#S16_wshed[:,0] = pl.flipud(AD_line[:,1]); S16_wshed[:,1] = pl.flipud(AD_line[:,0])

# --------------------------------------------------------------------------------
# !!! segment labels and initial rlspts arrays come from '_traj_release_ files !!!
# --------------------------------------------------------------------------------

labs1, rlspts1 = TrajSegLabel('Am'); labs2, rlspts2 = TrajSegLabel('AfMe')
labs3, rlspts3 = TrajSegLabel('EAA'); labs4, rlspts4 = TrajSegLabel('SO')
labs5, rlspts5 = TrajSegLabel('Ar')#; labs6, rlspts6 = TrajSegLabel('LS15full')

#seglab = pl.zeros([len(labs1)+len(labs2)])
#seglab[:len(labs1)] = labs1[:]; seglab[len(labs1):] = labs2[:]
#rlspts = pl.zeros([len(rlspts1)+len(rlspts2),2])
#rlspts[:len(rlspts1)] = rlspts1[:]; rlspts[len(rlspts1):] = rlspts2[:]
#f = S16_wshed; f[:,0] = f[:,0] + 360.


# ----------------------------------------------------------------------
# !!! check if any longitude co-ordinates are less than zero and fix !!!
# ----------------------------------------------------------------------
rlspts1 = Add360(rlspts1); rlspts2 = Add360(rlspts2)
rlspts3 = Add360(rlspts3); rlspts4 = Add360(rlspts4)
rlspts5 = Add360(rlspts5); rlspts6 = Add360(rlspts6)

# -------------------------------------------------------
# !!!!! calculate normal vectors from endpts arrays !!!!!
# -------------------------------------------------------
nhat1 = NormalVector(endpts1); nhat2 = NormalVector(endpts2)
nhat3 = NormalVector(endpts3); nhat4 = NormalVector(endpts4)
nhat5 = NormalVector(endpts5); nhat6 = NormalVector(endpts6)
#S16_normals = NormalVector(S16_wshed)



# -----------------------------------------------------
# !!!!! remove repeated points and segment labels !!!!!
# -----------------------------------------------------
rlspts1, labs1 = RR2(rlspts1,labs1); rlspts2, labs2 = RR2(rlspts2,labs2)
rlspts3, labs3 = RR2(rlspts3,labs3); rlspts4, labs4 = RR2(rlspts4,labs4);
rlspts5, labs5 = RR2(rlspts5,labs5); rlspts6, labs6 = RR2(rlspts6,labs6)

# ---------------------------------------------------
# !!!!!! calculate distances between midpoints !!!!!!
# ---------------------------------------------------
midpt1 = MidPts(rlspts1); midpt2 = MidPts(rlspts2); midpt3 = MidPts(rlspts3)
midpt4 = MidPts(rlspts4); midpt5 = MidPts(rlspts5); midpt6 = MidPts(rlspts6)
#S16_midpt = pl.zeros([S16_normals.shape[0]])
#for i in range(S16_midpt.shape[0]):
#    S16_midpt[i] = Haversine(S16_wshed[i],S16_wshed[i+1])

#flux_uv1 = pl.zeros_like(rlspts1); 
#for i in range(len(flux_uv)):
#    #a = NearestIndex(lon,S16_wshed[i,0]); b = NearestIndex(lat,S16_wshed[i,1])
#    #if lon[a] == S16_wshed[i,0]:
#    #    flux_uv[i,0] = tst_zon[b,a]; flux_uv[i,1] = tst_mer[b,a]
#    #else:
#    flux_uv[i,0] = BilinInterp(rlspts[i],lon,lat,zon_mean)
#    flux_uv[i,1] = BilinInterp(rlspts[i],lon,lat,mer_mean)
#

# -------------------------------------------
# !!!!!! calculate qu, qv at midpoints !!!!!!
# -------------------------------------------
midflux1 = MidFlux(rlspts1,lon,lat,zon_mean,mer_mean)
midflux2 = MidFlux(rlspts2,lon,lat,zon_mean,mer_mean)
midflux3 = MidFlux(rlspts3,lon,lat,zon_mean,mer_mean)
midflux4 = MidFlux(rlspts4,lon,lat,zon_mean,mer_mean)
midflux5 = MidFlux(rlspts5,lon,lat,zon_mean,mer_mean)
midflux6 = MidFlux(rlspts6,lon,lat,zon_mean,mer_mean)
#S16_midflux = MF2(S16_wshed,lon,lat,zon_mean,mer_mean)

# ------------------------------------------------
# !!!!!! calculate flux normal to midpoints !!!!!!
# ------------------------------------------------
FdotN1 = NormalFlux(midflux1,labs1,nhat1); FdotN2 = NormalFlux(midflux2,labs2,nhat2)
FdotN3 = NormalFlux(midflux3,labs3,nhat3); FdotN4 = NormalFlux(midflux4,labs4,nhat4)
FdotN5 = NormalFlux(midflux5,labs5,nhat5); FdotN6 = NormalFlux(midflux6,labs6,nhat6)
#FdotN_S16 = pl.zeros([S16_midflux.shape[0]])
#for i in range(S16_midflux.shape[0]):
#    FdotN_S16[i] = pl.dot(S16_midflux[i],S16_normals[i])

# ------------------------------------------------------------------------------------------
# !!! weight flux normal to midpoints by distances between midpoints and rescale by 10^9 !!!
# ------------------------------------------------------------------------------------------

fluxes1 = FdotN1[:]*midpt1[:]/(10**9); fluxes2 = FdotN2[:]*midpt2[:]/(10**9)
fluxes3 = FdotN3[:]*midpt3[:]/(10**9); fluxes4 = FdotN4[:]*midpt4[:]/(10**9)
fluxes5 = FdotN5[:]*midpt5[:]/(10**9); fluxes6 = FdotN6[:]*midpt6[:]/(10**9)
#S16_fluxes = FdotN_S16[:]*S16_midpt[:]/(10**9)

# -----------------------------------------------------------------------------
# !!! sum fluxes to get net moisture flux across catchment boundaries in Sv !!!
# -----------------------------------------------------------------------------
print 'Americas = ', pl.sum(fluxes1[:]), ' Sv'
print 'Africa & Middle East = ', pl.sum(fluxes2[:]), ' Sv'
print 'East Asia & Australia = ', pl.sum(fluxes3[:]), ' Sv'
print 'Southern Ocean = ', pl.sum(fluxes4[:]), ' Sv'
print 'Arctic = ', pl.sum(fluxes5[:]), ' Sv'
print 'Levang & Schmitt = ', pl.sum(fluxes6[:]), ' Sv'
#print 'Singh et al (2016) = ', pl.sum(S16_fluxes[:]), ' Sv'
#pl.close(); pl.close()


amr_ind1 = NearestIndex(lat,rlspts1[0,1]); amr_ind2 = NearestIndex(lat,rlspts1[-1,1])
afr_ind1 = NearestIndex(lat,rlspts2[0,1]); afr_ind2 = NearestIndex(lat,rlspts2[-1,1])
eaa_ind1 = NearestIndex(lat,rlspts3[0,1]); eaa_ind2 = NearestIndex(lat,rlspts3[-1,1])

amr_zm = pl.mean(zon_mean[amr_ind1:amr_ind2+1],axis=1)
afr_zm = pl.mean(zon_mean[afr_ind1:afr_ind2+1],axis=1)
eaa_zm = pl.mean(zon_mean[eaa_ind1:eaa_ind2+1],axis=1)
ZM = pl.mean(zon_mean,axis=1)

amr_pts = pl.zeros([len(amr_zm),2]); afr_pts = pl.zeros([len(afr_zm),2])
eaa_pts = pl.zeros([len(eaa_zm),2])

amr_pts[:,1] = lat[amr_ind1:amr_ind2+1]; afr_pts[:,1] = lat[afr_ind1:afr_ind2+1]
eaa_pts[:,1] = lat[eaa_ind1:eaa_ind2+1]
amr_mp = MidPts(amr_pts); afr_mp = MidPts(afr_pts); eaa_mp = MidPts(eaa_pts)

amr_zi = pl.sum(amr_zm*amr_mp)/(10**9); afr_zi = pl.sum(afr_zm*afr_mp)/(10**9)
eaa_zi = pl.sum(eaa_zm*eaa_mp)/(10**9)
print amr_zi, ' Sv'
print afr_zi, ' Sv'
print eaa_zi, ' Sv'

fig,ax = pl.subplots(1,3,figsize=(12,8))

ax1 = pl.subplot(131)
ax1.plot(ZM,lat,lw=2,color='dimgray',ls='--')
ax1.plot(FdotN1,rlspts1[:,1],lw=2,color='b')
pl.ylim(-40,60); pl.xlim(-300,200)
ax1.set_yticks(pl.arange(-40,61,10)); pl.yticks(fontsize=16)
ax1.grid(axis='y'); ax1.axvline(x=0,color='k',ls='--')
ax1.set_ylabel('latitude',fontsize=22,labelpad=-12)
pl.xticks(fontsize=14); ax1.set_xlabel('kg/m/s',fontsize=22)
ax1.tick_params(axis='x', which='major', pad=10)
ax1.annotate('(a) Americas',(-275,54),fontsize=18,
             bbox=dict(boxstyle='round',facecolor='w',edgecolor='k',alpha=0.9))
ax1.annotate('-0.23 Sv',(-275,41),fontsize=20,color='b')
ax1.annotate('-0.24 Sv',(-275,31),fontsize=20,color='dimgrey')#,bbox={'facecolor':'w'})
#pl.title('(a)',fontsize=20)

ax2 = pl.subplot(132)
ax2.plot(ZM,lat,lw=2,color='dimgray',ls='--')
ax2.plot(FdotN2,rlspts2[:,1],lw=2,color='k')
pl.ylim(-40,60); pl.xlim(-300,200)
ax2.set_yticks(pl.arange(-40,61,10))
ax2.grid(axis='y'); ax2.axvline(x=0,color='k',ls='--')
ax2.yaxis.set_major_formatter(pl.NullFormatter())
pl.xticks(fontsize=14); ax2.set_xlabel('kg/m/s',fontsize=22)
ax2.tick_params(axis='x', which='major', pad=10)
ax2.annotate('(b) Africa',(-275,54),fontsize=18,
             bbox=dict(boxstyle='round',facecolor='w',edgecolor='k',alpha=0.9))
ax2.annotate('-0.16 Sv',(-275,41),fontsize=20,color='k')
ax2.annotate('-0.36 Sv',(-275,31),fontsize=20,color='dimgrey')
#pl.title('(b)',fontsize=20)

ax3 = pl.subplot(133)
ax3.plot(ZM,lat,lw=2,color='dimgray',ls='--')
ax3.plot(FdotN3,rlspts3[:,1],lw=2,color='g')
pl.ylim(-40,60); pl.xlim(-300,200)
ax3.set_yticks(pl.arange(-40,61,10))
ax3.grid(axis='y'); ax3.axvline(x=0,color='k',ls='--')
ax3.yaxis.set_major_formatter(pl.NullFormatter())
pl.xticks(fontsize=14); ax3.set_xlabel('kg/m/s',fontsize=22)
ax3.tick_params(axis='x', which='major', pad=10)
ax3.annotate('(c) South-East Asia',(-275,54),fontsize=18,
             bbox=dict(boxstyle='round',facecolor='w',edgecolor='k',alpha=0.9))
ax3.annotate('0.16 Sv',(-275,41),fontsize=20,color='g')
ax3.annotate('-0.50 Sv',(-275,31),fontsize=20,color='dimgrey')
#pl.title('(c)',fontsize=20)

pl.subplots_adjust(left=0.07,right=0.96)

#fix, ax = pl.subplots(3,3,figsize=(15,10))
#
#ax1 = pl.subplot(331)
#pl.plot(FdotN1,lw=2.,color='b'); pl.xlim(0,len(fluxes1)); pl.ylim(-300,100)
#pl.text(3,-290,'(a) Americas',fontsize=20); pl.ylabel('kg/m/s',fontsize=22)
#pl.yticks(fontsize=13); pl.grid(axis='y')
#a = pl.arange(0,155,20); b = pl.around(rlspts1[a,1],decimals=0)
#ax1.set_xticks(a); ax1.set_xticklabels(b.astype('int'),fontsize=14)
#ax1.tick_params(axis='y',pad=5); ax1.tick_params(axis='x',pad=7)
#
#ax2 = pl.subplot(332)
#pl.plot(FdotN2,lw=2.,color='k'); pl.xlim(0,len(fluxes2)); pl.ylim(-300,100)
#pl.text(3,-290,'(b) Africa & Middle East',fontsize=20)
#ax2.yaxis.set_major_formatter(pl.NullFormatter())
#pl.grid(axis='y')
#a = pl.arange(0,214,30); b = pl.around(rlspts2[a,1],decimals=0)
#ax2.set_xticks(a); ax2.set_xticklabels(b.astype('int'),fontsize=14)
#pl.xlabel('latitude',fontsize=16); ax2.tick_params(axis='x',pad=7)
#
#ax3 = pl.subplot(333)
#pl.plot(FdotN3,lw=2.,color='g'); pl.xlim(0,len(fluxes3)); pl.ylim(-300,100)
#pl.text(3,-290,'(c) South-East Asia',fontsize=20)
#ax3.yaxis.set_major_formatter(pl.NullFormatter())
#pl.grid(axis='y')
#a = pl.arange(0,165,20); b = pl.around(rlspts3[a,1],decimals=0)
#ax3.set_xticks(a); ax3.set_xticklabels(b.astype('int'),fontsize=14)
#ax3.tick_params(axis='x',pad=7)
#
#ax4 = pl.subplot(334)
#pl.plot(FdotN5[62:259],lw=2.,color='r')
#pl.xlim(0,len(fluxes5[62:259])); pl.ylim(-60,60)
#pl.text(80,-55,'(d) Arctic Atlantic',fontsize=20); pl.ylabel('kg/m/s',fontsize=22)
#pl.grid(axis='y'); pl.yticks(fontsize=13)
#R = rlspts5[62:259,0]; R[:95] = R[:95] - 360
#a = pl.arange(0,197,30); b = pl.around(R[a],decimals=0)
#ax4.set_xticks(a); ax4.set_xticklabels(b.astype('int'),fontsize=14)
#ax4.tick_params(axis='y',pad=5); ax4.tick_params(axis='x',pad=7)
#
#ax5 = pl.subplot(335)
#pl.plot(FdotN5[259:361],lw=2.,color='r')
#pl.xlim(0,len(fluxes5[259:361])); pl.ylim(-60,60)
#pl.text(48,-55,'(e) Arctic Indian',fontsize=20)
#ax5.yaxis.set_major_formatter(pl.NullFormatter())
#pl.grid(axis='y')
#R = rlspts5[259:361,0]
#a = pl.arange(0,102,20); b = pl.around(R[a],decimals=0)
#ax5.set_xticks(a); ax5.set_xticklabels(b.astype('int'),fontsize=14)
#pl.xlabel('longitude',fontsize=16); ax5.tick_params(axis='x',pad=7)
#
#ax6 = pl.subplot(336)
#flx_arp = pl.zeros([len(fluxes5[361:])+len(fluxes5[:62])])
#flx_arp[:(len(rlspts5)-361)] = FdotN5[361:]; flx_arp[len(rlspts5)-361:] = FdotN5[:62]
#pl.plot(flx_arp,lw=2.,color='r')
#pl.xlim(0,len(flx_arp)); pl.ylim(-60,60)
#pl.text(105,-55,'(f) Arctic Pacific',fontsize=20)
#ax6.yaxis.set_major_formatter(pl.NullFormatter()); pl.grid(axis='y')
#R = pl.zeros([len(flx_arp)])
#R[:len(rlspts5)-361] = rlspts5[361:,0]; R[len(rlspts5)-361:] = rlspts5[:62,0]
#R[148:] = R[148:] - 360.
#a = pl.arange(0,220,30); b = pl.around(R[a],decimals=0)
#ax6.set_xticks(a); ax6.set_xticklabels(b.astype('int'),fontsize=14)
#ax6.tick_params(axis='x',pad=7)
#
#ax7 = pl.subplot(337)
#pl.plot(FdotN4[394:],lw=2.,color='darkgoldenrod')
#pl.xlim(0,len(fluxes4[394:])); pl.ylim(-200,100)
#pl.text(55,-180,'(g) Southern Atlantic',fontsize=20); pl.ylabel('kg/m/s',fontsize=22)
#pl.grid(axis='y'); pl.yticks(fontsize=13)
#R = rlspts4[394:,0]; R[:161] = R[:161] - 360.
#
#ax7.tick_params(axis='y',pad=5); ax7.tick_params(axis='x',pad=7)
#
#ax8 = pl.subplot(338)
#pl.plot(FdotN4[:184],lw=2.,color='darkgoldenrod')
#pl.xlim(0,len(fluxes4[:184])); pl.ylim(-200,100)
#pl.text(65,-180,'(h) Southern Indian',fontsize=20)
#ax8.yaxis.set_major_formatter(pl.NullFormatter()); pl.grid(axis='y')
#R = rlspts4[:184,0]
#a = pl.arange(0,184,30); b = pl.around(R[a],decimals=0)
#ax8.set_xticks(a); ax8.set_xticklabels(b.astype('int'),fontsize=14)
#pl.xlabel('longitude',fontsize=16); ax8.tick_params(axis='x',pad=7)
#
#ax9 = pl.subplot(339)
#pl.plot(FdotN4[184:394],lw=2.,color='darkgoldenrod')
#pl.xlim(0,len(fluxes4[184:394])); pl.ylim(-200,100)
#pl.text(78,-180,'(i) Southern Pacific',fontsize=20)
#ax9.yaxis.set_major_formatter(pl.NullFormatter()); pl.grid(axis='y')
#R = rlspts4[184:394,0]; R[62:210] = R[62:210] - 360.
#a = pl.arange(0,210,30); b = pl.around(R[a],decimals=0)
#ax9.set_xticks(a); ax9.set_xticklabels(b.astype('int'),fontsize=14)
#ax9.tick_params(axis='x',pad=7)
#
#pl.subplots_adjust(wspace=0.09)
#pl.tight_layout()

#profiles = pl.zeros([9,len(flx_arp)])
#profiles[0,:len(FdotN1)] = FdotN1; profiles[1,:len(FdotN2)] = FdotN2
#profiles[2,:len(FdotN3)] = FdotN3
#profiles[3,:len(FdotN5[62:259])] = FdotN5[62:259]
#profiles[4,:len(FdotN5[259:361])] = FdotN5[259:361]
#profiles[5,:len(flx_arp)] = flx_arp
#profiles[6,:len(FdotN4[394:])] = FdotN4[394:]
#profiles[7,:len(FdotN4[:184])] = FdotN4[:184]
#profiles[8,:len(FdotN4[184:394])] = FdotN4[184:394]
#
#f = open(sheddir+'profiles_1014.csv','w')
#pl.savetxt(f,profiles.T,delimiter=',')
#f.close()

#soa_ave = FluxPerLength(rlspts4[489:],fluxes4[489:])
#soi_ave = FluxPerLength(rlspts4[:112],fluxes4[:112])
#sop_ave = FluxPerLength(rlspts4[213:362],fluxes4[213:362])

#print soa_ave, ' m**2/s'
#print soi_ave, ' m**2/s'
#print sop_ave, ' m**2/s'

#fig,ax = pl.subplots(2,1)
#
#ax1 = pl.subplot(211)
#m = Basemap(projection='cyl',lon_0=-160.); m.drawcoastlines()
#mer_shft, lons = shiftgrid(20.0, mer_mean, lon, start=False)
#newlon, newlat = pl.meshgrid(lons,lat)
#X, Y = m(newlon,newlat)
#levels = [-150,-100,-50,-25,-10,10,25,50,100,150]
#cs = m.contourf(X,Y,mer_shft,cmap='seismic',norm=pl.Normalize(-150,150),levels=levels,extend='both')
#m.colorbar()
#m.plot(-123.59,-35,latlon=True,color='yellow',marker='x',markersize=12,mew=3)
#m.plot(-25.42,-35,latlon=True,color='yellow',marker='x',markersize=12,mew=3)
#m.plot(rlspts4[489:-1,0],rlspts4[489:-1,1],latlon=True,color='k',lw=2.)
#m.plot(rlspts4[213:364,0],rlspts4[213:364,1],latlon=True,color='k',lw=2.)
#
#ax2 = pl.subplot(212)
#m = Basemap(projection='cyl',lon_0=-160.); m.drawcoastlines()
#zon_shft, lons = shiftgrid(20.0, zon_mean, lon, start=False)
#newlon, newlat = pl.meshgrid(lons,lat)
#X, Y = m(newlon,newlat)
#levels = [-250,-200,-100,-50,-25,-10,10,25,50,100,200,250]
#cs = m.contourf(X,Y,zon_shft,cmap='seismic',norm=pl.Normalize(-250,250),levels=levels,extend='both')
#m.colorbar()
#m.plot(rlspts4[489:-1,0],rlspts4[489:-1,1],latlon=True,color='k',lw=2.)
#m.plot(rlspts4[213:363,0],rlspts4[213:363,1],latlon=True,color='k',lw=2.)
#
#mer_satl = pl.zeros([rlspts4[489:-1].shape[0]])
#zon_satl = pl.zeros_like(mer_satl)
#mer_spac = pl.zeros([rlspts4[213:363].shape[0]])
#zon_spac = pl.zeros_like(mer_spac)
#
#for i in range(mer_satl.shape[0]):
#    mer_satl[i] = BilinInterp(rlspts4[489:-1][i],lon,lat,mer_mean)
#    zon_satl[i] = BilinInterp(rlspts4[489:-1][i],lon,lat,zon_mean)
#
#for j in range(mer_spac.shape[0]):
#    mer_spac[j] = BilinInterp(rlspts4[213:363][j],lon,lat,mer_mean)
#    zon_spac[j] = BilinInterp(rlspts4[213:363][j],lon,lat,zon_mean)
#
#pl.figure(2)
#pl.plot(mer_satl,lw=2.,label='Atl $qv$')
#pl.plot(zon_satl,lw=2.,label='Atl $qu$')
#pl.plot(mer_spac,lw=2.,label='Pac $qv$')
#pl.plot(zon_spac,lw=2.,label='Pac $qu$')
#pl.axhline(y=0,color='k',ls='--')
#pl.legend(ncol=2)
#pl.ylabel('kg m$^{-1}$ s$^{-1}$',fontsize=22,labelpad=-3.5)
