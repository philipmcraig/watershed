"""
Read in etopo05 topographic data and overlay it with HydroSHEDS global watershed data.
Watersheds and topography are used to indicate continental divides between oceans.
Click on resulting plot to output figure co-ordinates of continental divides,
use a Basemap fucntion to convert these to longitude-latitude co-ordinates and
output them as an array (then print to a file, which will be done later).

Last updated: 20/12/16 11:42 AM 12th commit
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans
from scipy.interpolate import interp1d
from netCDF4 import Dataset
import itertools
#import sys
#sys.path.append(r'/home/np838619/PminusE_data/ERA_Int/functions')
#from functions import *

exec(open('/home/users/qx911590/np838619/PminusE_data/ERA_Int/functions.py').read())


def SegmentLength(point1,point2,shedmap):
    """Function to find the length (in metres) of a line segment.
    
    Args:
        point1 (array): x & y co-ordinates of a clicked point in longitude &
                        latitude (degrees)
        point2 (array): x & y co-ordinates of a clicked point in longitude &
                        latitude (degrees)
        shedmap (Basemap object): map used for clicking
    
    Returns:
        length (float): length of a line segment (in metres)
    """
    # convert the longitude/latitude co-ordinates of the line end points into
    # co-ordinates in metres
    a, b = shedmap(point1[0],point1[1]) # end point 1
    c, d = shedmap(point2[0],point2[1]) # end point 2
    # Use Pythagoras Theorem to calculate the length of the line segment
    length = pl.sqrt((c-a)**2 + (d-b)**2)
    
    return length

def Haversine(pt1,pt2):
    """
    """
    pt1 = pl.radians(pt1); pt2 = pl.radians(pt2)
    
    phi1 = pt1[1]; phi2 = pt2[1]; dphi = phi2 - phi1
    lam1 = pt1[0]; lam2 = pt2[0]; dlam = lam2 - lam1
    
    dsig = 2*pl.arcsin(pl.sqrt((pl.sin(dphi/2)**2 + pl.cos(phi1)*pl.cos(phi2)*pl.sin(dlam/2)**2)))
    
    R = 6.37e6
    
    d = R*dsig
    
    return d

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
    loccar[:,0] = R*pl.cos(coords[:,1])*pl.cos(coords[:,0])
    loccar[:,1] = R*pl.cos(coords[:,1])*pl.sin(coords[:,0]) # y = R*lat
    
    return loccar

def SegmentInterp(point1,point2,seg_length,spacing):
    """Function to perform linear interpolation along a line segment on a map.
    
    Args:
        point1 (array): x & y co-ordinates of a clicked point in longitude &
                        latitude (degrees)
        point2 (array): x & y co-ordinates of a clicked point in longitude &
                        latitude (degrees)
        seg_length (float): lenth of line segment
        spacing (float): required distance between interpolated points
    
    Returns:
        R (array): co-ordinates of interpolated points
    """
    #x_points = (point1[0],point2[0]) # get the x co-ords of each point
    #y_points = (point1[1],point2[1]) # get the y co-ords of each point
    x1 = point1[0]; x2 = point2[0]
    y1 = point1[1]; y2 = point2[1]
    #point1 = pl.radians(point1); point2 = pl.radians(point2)
    #x1, y1 = shedmap(point1[0],point1[1])
    #x2, y2 = shedmap(point2[0],point2[1])
    #R = 6.37*10**6
    #x1 = R*pl.cos(point1[1])*pl.cos(point1[0])
    #y1 = R*pl.cos(point1[1])*pl.sin(point1[0])
    #x2 = R*pl.cos(point2[1])*pl.cos(point2[0])
    #y2 = R*pl.cos(point2[1])*pl.sin(point2[0])
    
    # Use scipy interp1d function to get interpolation function:
    interp_func = interp1d((x1,x2),(y1,y2))
    
    points = seg_length/spacing + 1; points = pl.ceil(points)
    
    # get an array of 10 points between the two x co-ordinates:
    if x1 == x2:
        # Use scipy interp1d function to get interpolation function:
        interp_func = interp1d((y1,y2),(x1,x2))
        y_new = pl.linspace(y1,y2,points)
        x_new = interp_func(y_new)
    else:
        # Use scipy interp1d function to get interpolation function:
        interp_func = interp1d((x1,x2),(y1,y2))
        x_new = pl.linspace(x1,x2,points)
        # use the interpolation function to calculate the corresponding new y co-ords
        y_new = interp_func(x_new)
    
    R = pl.zeros([len(x_new),2]) # make empty array for figure co-ordinates
    R[:,0] = x_new[:]; R[:,1] = y_new[:]
    
    return R

def RemoveRepeated(interp_pts):
    """Function to remove repeated release points.
    
    Args:
        interp_pts (array): nX2 array of release points, some of which are
                            repeated
    
    Returns:
        rp2 (array): mX2 array of release points, where m < n
    """
    # empty array for lengths between release pts:
    rel_len = pl.zeros([len(interp_pts)-1])
    for l in range(len(interp_pts)-1): 
        rel_len[l] = pl.sqrt((interp_pts[l+1,0] - interp_pts[l,0])**2 + 
                                (interp_pts[l+1,1] - interp_pts[l,1])**2)
    
    # Remove repeated points:
    rp2 = [interp_pts[0]] # include 1st release point
    for l in range(len(rel_len)):
        if rel_len[l] != 0.: # check for nonzero length
            rp2.append(interp_pts[l+1]) # add next point to list
    rp2 = pl.asarray(rp2) # convert list to array
    
    return rp2

def RR2(interp_pts):
    repeat = []
    for r in range(1,interp_pts.shape[0]):
        if interp_pts[r,0] == interp_pts[r-1,0] and interp_pts[r,1] == interp_pts[r-1,1]:
            repeat.append(r)
    
    RP2 = []
    for r in range(interp_pts.shape[0]):
        if r not in repeat:
            RP2.append(interp_pts[r])
    rlspts = pl.asarray(RP2)
    
    return rlspts

def XYtoLL(point_in):
    """
    """
    R = 6.37e6
    lam = pl.arctan2(point_in[1],point_in[0]); lon_out = pl.degrees(lam)
    phi = pl.arccos(point_in[0]/(R*pl.cos(lam))); lat_out = pl.degrees(phi)
    #phi = point_in[1]/R; lat_out = pl.degrees(phi)
    #lam = point_in[0]/(R*pl.cos(phi)); lon_out = pl.degrees(lam)
    
    point_out = pl.array([lon_out,lat_out])
    
    return point_out

pl.close('all')
sheddir = '/home/users/qx911590/np838619/Watershed/'

# read in topography data
ncfile = Dataset('/home/users/qx911590/np838619/etopo05.nc','r')
topo = ncfile.variables['ROSE'][:] #perhaps extract only values > 0?
lat = ncfile.variables['ETOPO05_Y'][:]
lon = ncfile.variables['ETOPO05_X'][:]
ncfile.close()

#ncfile = Dataset(sheddir+'ggis198901010000.nc','r')
#geop = ncfile.variables['Z'][:]
#ncfile.close()
#geop = pl.squeeze(geop)

era = Dataset(sheddir+'ggis198901010000.nc','r')
eralat = era.variables['latitude'][:]
eralon = era.variables['longitude'][:]
era.close()
lonshift = eralon-180.

# get rid of all values below zero, bathymetry unimportant:
#for i in range(topo.shape[0]):
#    for j in range(topo.shape[1]):
#        if topo[i,j] < 0.:
#            topo[i,j] = pl.float64('nan')
#topo = geop/9.81
topo_mask = pl.ma.masked_where(topo < 0.,topo)
#topo_mask = pl.ma.masked_invalid(topo_mask)

lon_0 = lon.mean(); lat_0 = lat.mean()

bounds = pl.array([88,-30,155,20]) # East Asia/Australia, Indian/Pacific divide
isect = pl.array([[106.15417,-6.06455], # west Java
                  [106.15417,-24.9273], # Indian Ocean
                  [114.803,-24.9273], # West Australia 1
                  [115.792,-25.679], # West Australia 2
                  [119.28654,-25.89960], # West Australia 3
                  [137.29319,-18.87625], # Oz_clicks -3
                  [141.40260,-22.16378], # Oz_clicks -2
                  [144.970,-20.480]  # Oz_clicks -1
                    ])
label = 'Rod11'
name = 'Rodriguez et al. (2011) Indian Ocean/Australia'
# select which map to use:
#shed = 2 # shed corresponds to index of 0-axis in bounds array
print name + ' Watershed'

#if shed in range(0,5):
m = Basemap(llcrnrlon=bounds[0],llcrnrlat=bounds[1],urcrnrlon=bounds[2],
            urcrnrlat=bounds[3],lat_ts=20.,resolution='l',projection='cyl',
            suppress_ticks=True)


topo, lons = shiftgrid(180.0, topo_mask, lon, start=False)
newlon, newlat = pl.meshgrid(lons,lat)
X, Y = m(newlon,newlat)
topo_land = maskoceans(newlon,newlat,topo,inlands=False,resolution='l',grid=5)

cmap = pl.get_cmap('YlOrRd')
#m.fillcontinents(color=None,lake_color='white',zorder=1)
cs = m.pcolormesh(X,Y,topo_land,cmap=cmap)
m.drawcoastlines(linewidth=0.1)


#--------------------READ IN ALL THE SHAPEFILES--------------------------------

as_shp = m.readshapefile("as_bas_30s_beta","as",linewidth=0.5) # rest of Asia
au_shp = m.readshapefile("au_bas_30s_beta","au",linewidth=0.5,color='lightgrey') # Australasia
oz_shp = m.readshapefile(sheddir+"Oz_sheds/oz_sheds","oz",linewidth=0.5)
m.plot(isect[:,0],isect[:,1],'ro')#; m.plot(isect[5,0],isect[5,1],'ro')
    #au_shp = m.readshapefile("/home/np838619/Watershed/42343_shp/rbasin_polygon","au",linewidth=0.5)
pl.tight_layout()


#-----------------------CLICK TO ALONG CONTINENTAL DIVIDES-----------------------
print 'Click along watershed. Click middle mouse button to stop.'
coords = pl.ginput(n=6,timeout=0)
pl.show()

# Use MapCoords function to get the longitude/latitude co-ordinates of points:
coords = pl.asarray(coords)
lonlat = pl.zeros([isect.shape[0]+coords.shape[0],2])#pl.asarray(coords)#lonlat = MapCoords(coords,m)
lonlat[:5] = isect[:5]
lonlat[5:coords.shape[0]+5] = coords[:]
lonlat[5+coords.shape[0]:] = isect[5:]
#lonlat[-1] = isect[-1]

#lonlat[-2,0] = lonlat[-1,0]
#print lonlat

# Write clicked points to a file:
f = open(sheddir+'shed_defs/'+label +'_clicks.txt','w')
f.write('Longitude, Latitude co-ordinates of clicked points discretely defining the '
        + name + ' continental divide. \nUsed to release back trajectories. \n \n' 
        + str(len(lonlat)) + '  points defined in this file. \n \n')
pl.savetxt(f,lonlat,fmt='%9.5f')
f.close()


#boundary = pl.zeros_like(lonlat) # set up empty array same size as lonlat
# Get the figure co-ordinates of the clicks:
#boundary[:,0], boundary[:,1] = m(lonlat[:,0],lonlat[:,1],inverse=False)
# PLot the points and line segments:
m.plot(lonlat[:,0],lonlat[:,1],'b') # add line segments
m.plot(lonlat[:,0],lonlat[:,1],'ro') # add markers
m.drawparallels([30,45,60],labels=[1,0,0,1],linewidth=0,ax=None)
m.drawmeridians([-60,-90,-120,-150],labels=[1,0,0,1],linewidth=0,ax=None)
pl.savefig(sheddir+'shed_defs/'+label +'_clicks.png')

length = [] # set up empty list for segment lengths
for i in range(lonlat.shape[0]-1): # loop over number of segments
    #length.append(SegmentLength(lonlat[i],lonlat[i+1],m)) # calculate segment length
    length.append(Haversine(lonlat[i],lonlat[i+1]))
 
len_tot = pl.sum(length) # calculate the total length of watershed
print 'Total length of watershed: ', pl.sum(length), ' metres'
spacing = 75000. # impose spacing of 79km between points
#no_of_points = 200 # number of interpolated points required
# Using np_of_points, what should the spacing between each point be?
#spacing = len_tot/no_of_points
print 'Spacing between trajectory release points: ', spacing, ' metres'


F = [] # empty list for interpolated segments 
for i in range(len(length)): # loop over number of segments
    # interpolate between the end points of each segment & append each segment to list::
    F.append(SegmentInterp(lonlat[i],lonlat[i+1],length[i],spacing))

# List of interpolated points is currently a list of length no. of segments,
# with each element an array of points
# Following commands flatten the list to be length no. of points with each element
# an array with the x & y co-ordinate of each interpolated point
D = itertools.chain(*F)
S = list(D)

interp_pts = pl.zeros([len(S),len(S[0])]) # empty array for interpolated points
interp_pts[:] = S[:] # stick interpolated points into array


f = open(sheddir+'shed_defs/'+label+'_traj_release.txt','w')
f.write('Longitude, Latitude co-ordinates defining the ' + name + 
        ' continental divide. \nUsed to release back trajectories. \n \n' 
        + str(len(interp_pts)) + '  points defined in this file. \n \n')
pl.savetxt(f,interp_pts,fmt='%9.5f')
f.close()

# Remove repeated release points:
rp2 = RR2(interp_pts)#; points = pl.zeros_like(rp2)
#for i in range(rp2.shape[0]):
#    points[i] = XYtoLL(rp2[i])

# Convert interpolated points to lon-lat co-ordinates before writing to file:
#interp_lonlat = pl.zeros_like(rp2)
m.plot(rp2[:,0],rp2[:,1],'go') # plot the interpolated points
pl.savefig(sheddir+'shed_defs/'+label+'_interp_pts.png')
#interp_lonlat[:,0] = rp2[:,0]; interp_lonlat[:,1] = rp2[:,1]


# Write interpolated points to a file:
f = open(sheddir+'shed_defs/'+label+'_traj_release_new.txt','w')
f.write('Longitude, Latitude co-ordinates defining the ' + name + 
        ' continental divide. \nUsed to release back trajectories. \n \n' + 
        str(len(rp2)) + '  points defined in this file. \n \n')
pl.savetxt(f,rp2,fmt='%9.5f')
f.close()


#start_time = timeit.default_timer()
#h = RemoveRepeated(interp_pts)
#elapsed = timeit.default_timer() - start_time
#print elapsed