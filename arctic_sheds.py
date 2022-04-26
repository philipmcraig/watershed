"""
Arctic Ocean drainage divide
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

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())


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
    R = pl.zeros([points,2]) # make empty array for figure co-ordinates
    
    if x1 == x2:
        y_new = pl.linspace(y1,y2,points)
        R[:,0] = x1; R[:,1] = y_new[:]
    else:
        # get an array of points between the two x co-ordinates:
        x_new = pl.linspace(x1,x2,points)
        # use the interpolation function to calculate the corresponding new y co-ords
        y_new = interp_func(x_new)
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
sheddir = '/home/np838619/Watershed/'
gldir = sheddir + 'GRE_Basins_IMBIE2_v13/'
cddir = sheddir + 'canadaoda_p_1m_v6-0_shp/'
#nadir = sheddir + 'NA_Watersheds/'

# read in topography data
ncfile = Dataset('/home/np838619/Downloads/etopo05.nc','r')
topo = ncfile.variables['ROSE'][:] #perhaps extract only values > 0?
lat = ncfile.variables['ETOPO05_Y'][:]
lon = ncfile.variables['ETOPO05_X'][:]
ncfile.close()

#ncfile = Dataset(sheddir+'ggis198901010000.nc','r')
#geop = ncfile.variables['Z'][:]
#ncfile.close()
#geop = pl.squeeze(geop)

era = Dataset(sheddir+'ggap200707211200.nc','r')
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

bounds = pl.array([[88,-26,155,33], # Scandanavia to Bering Strait
                    [-172,50,-50,85], # Bering Strait to Greenland
                    [-77,66,34,85], # Greenland to Scandanavia
                    ])
isect = pl.array([[21.27,69.32], # Northermost point of Scandanavian drainage divide
                    [41.98,40.08], # Arctic/Middle East
                    [92.18,32.71], # Arctic/East Asia
                    [-169.03,65.77], # Bering Strait
                    [-117.45,52.03], # Arctic/Americas, the "Snow Dome"
                    [-72.62,78.08], # Canada/Greenland
                    [-27.09,80.5], # Greenland/Norwegian Sea
                    [21.27,80.5] # Svalbard
                    ])
label = ['ScBS','BSGl','GlSc']
name = ['Scandanavia to Bering Strait','Bering Strait to Greenland','Greenland to Scandanavia']
# select which map to use:
shed = 2 # shed corresponds to index of 0-axis in bounds array
print name[shed]# + ' Watershed'
if shed == 0.: print '65 points needed!'
if shed == 1.: print '22 points needed!'
if shed == 2.: print '7 points needed!'

if shed > 0:
    m = Basemap(llcrnrlon=bounds[shed,0],llcrnrlat=bounds[shed,1],urcrnrlon=bounds[shed,2],
                urcrnrlat=bounds[shed,3],lat_ts=20.,resolution='l',projection='cyl',
                suppress_ticks=True)
elif shed == 0:
    m = Basemap(lat_0=60,lon_0=100.,lat_1=45,lat_2=55,width=10000000,height=7000000,
             resolution='l',projection='aea')
    m.drawcountries(linewidth=0.5,color='lightgray')
    #m = Basemap(llcrnrlon=bounds[shed,0],llcrnrlat=bounds[shed,1],urcrnrlon=bounds[shed,2],
    #            urcrnrlat=bounds[shed,3],lon_0=100.,
    #            o_lat_p=90.,o_lon_p=90.,resolution='l',projection='rotpole')

topo, lons = shiftgrid(180.0, topo_mask, lon, start=False)
newlon, newlat = pl.meshgrid(lons,lat)
X, Y = m(newlon,newlat)
topo_land = maskoceans(newlon,newlat,topo,inlands=False,resolution='l',grid=5)

cmap = pl.get_cmap('YlOrRd')
#m.fillcontinents(color=None,lake_color='white',zorder=1)
cs = m.pcolormesh(X,Y,topo_land,cmap=cmap)
m.drawcoastlines(linewidth=0.1)


#--------------------READ IN ALL THE SHAPEFILES--------------------------------
if shed == 0:
    as_shp = m.readshapefile("as_bas_30s_beta","as",linewidth=0.5) # rest of Asia
    eu_shp = m.readshapefile("eu_bas_30s_beta","eu",linewidth=0.5) # Europe + some of Asia
    m.plot(isect[0,0],isect[0,1],'ro',latlon=True)
    m.plot(isect[1,0],isect[1,1],'ro',latlon=True)
    m.plot(isect[2,0],isect[2,1],'ro',latlon=True)
    m.plot(isect[3,0],isect[3,1],'ro',latlon=True)
elif shed == 1.:
    cd_shp = m.readshapefile(cddir+"canadnaoda_p","cd",linewidth=0.5)
    m.plot(isect[3,0],isect[3,1],'ro'); m.plot(isect[4,0],isect[4,1],'ro')
    m.plot(isect[5,0],isect[5,1],'ro')
elif shed == 2.:
    gl_shp = m.readshapefile(gldir+"GRE_Basins_IMBIE2_v13","gl",linewidth=0.5)
    m.plot(isect[5,0],isect[5,1],'ro'); m.plot(isect[6,0],isect[6,1],'ro')
    m.plot(isect[7,0],isect[7,1],'ro'); m.plot(isect[0,0],isect[0,1],'ro')
    m.drawcountries(linewidth=0.5,color='lightgray')
#elif shed == 3.:
#    sa_shp = m.readshapefile("sa_bas_30s_beta","sa",linewidth=0.5) # Europe + some of Asia
#elif shed == 4.:
#    as_shp = m.readshapefile("as_bas_30s_beta","as",linewidth=0.5) # rest of Asia
#    au_shp = m.readshapefile("au_bas_30s_beta","au",linewidth=0.5) # Australasia
#elif shed == 5.:
#    as_shp = m.readshapefile("as_bas_30s_beta","as",linewidth=0.5) # rest of Asia
#    eu_shp = m.readshapefile("eu_bas_30s_beta","eu",linewidth=0.5) # Europe + some of Asia
pl.tight_layout()


#-----------------------CLICK TO ALONG CONTINENTAL DIVIDES-----------------------
print 'Click along watershed. Click middle mouse button to stop.'
coords = pl.ginput(n=0,timeout=0)
pl.show()

# Use MapCoords function to get the longitude/latitude co-ordinates of points:
lonlat = pl.asarray(coords)#lonlat = MapCoords(coords,m)
if shed == 0.:
    boundary = pl.zeros_like(lonlat)
    boundary[:,0], boundary[:,1] = m(lonlat[:,0],lonlat[:,1],inverse=True)
    boundary[0] = isect[0]; boundary[17] = isect[1]; boundary[34] = isect[2]
    boundary[-1] = isect[3]
    lonlat = boundary
    for i in range(lonlat.shape[0]):
        if lonlat[i,0] < 0.: lonlat[i,0] = lonlat[i,0] + 360.
if shed == 1.:
    lonlat[0] = isect[3]; lonlat[-1] = isect[5]; lonlat[9] = isect[4]
if shed == 2.:
    lonlat[0] = isect[5]; lonlat[-1] = isect[0]; lonlat[-2] = isect[-1]
    lonlat[-3] = isect[-2]

# Write clicked points to a file:
f = open(sheddir+'shed_defs/'+label[shed]+'_clicks.txt','w')
f.write('Longitude, Latitude co-ordinates of clicked points discretely defining the '
        + name[shed]+ ' part of the Arctic Ocean drainage divide. \nUsed to release back trajectories. \n \n' 
        + str(len(lonlat)) + '  points defined in this file. \n \n')
pl.savetxt(f,lonlat,fmt='%9.5f')
f.close()


#boundary = pl.zeros_like(lonlat) # set up empty array same size as lonlat
# Get the figure co-ordinates of the clicks:
#boundary[:,0], boundary[:,1] = m(lonlat[:,0],lonlat[:,1],inverse=False)
# PLot the points and line segments:
m.plot(lonlat[:,0],lonlat[:,1],'b',latlon=True) # add line segments
m.plot(lonlat[:,0],lonlat[:,1],'ro',latlon=True) # add markers
m.drawparallels([30,45,60],labels=[1,0,0,1],linewidth=0,ax=None)
m.drawmeridians([-60,-90,-120,-150],labels=[1,0,0,1],linewidth=0,ax=None)
pl.savefig(sheddir+'shed_defs/'+label[shed]+'_clicks.png')

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


f = open(sheddir+'shed_defs/'+label[shed]+'_traj_release.txt','w')
f.write('Longitude, Latitude co-ordinates defining the ' + name[shed] + 
        ' part of the Arctic Ocean drainage divide. \nUsed to release back trajectories. \n \n' 
        + str(len(interp_pts)) + '  points defined in this file. \n \n')
pl.savetxt(f,interp_pts,fmt='%9.5f')
f.close()

# Remove repeated release points:
rp2 = RR2(interp_pts)#; points = pl.zeros_like(rp2)
#for i in range(rp2.shape[0]):
#    points[i] = XYtoLL(rp2[i])

# Convert interpolated points to lon-lat co-ordinates before writing to file:
#interp_lonlat = pl.zeros_like(rp2)
m.plot(rp2[:,0],rp2[:,1],'go',latlon=True) # plot the interpolated points
pl.savefig(sheddir+'shed_defs/'+label[shed]+'_interp_pts.png')
#interp_lonlat[:,0] = rp2[:,0]; interp_lonlat[:,1] = rp2[:,1]


# Write interpolated points to a file:
f = open(sheddir+'shed_defs/'+label[shed]+'_traj_release_new.txt','w')
f.write('Longitude, Latitude co-ordinates defining the ' + name[shed] + 
        ' part of the Arctic Ocean drainage divide. \nUsed to release back trajectories. \n \n' + 
        str(len(rp2)) + '  points defined in this file. \n \n')
pl.savetxt(f,rp2,fmt='%9.5f')
f.close()


#start_time = timeit.default_timer()
#h = RemoveRepeated(interp_pts)
#elapsed = timeit.default_timer() - start_time
#print elapsed