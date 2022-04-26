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
    #endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
    #endpts = pl.asarray(endlist[5:],dtype='float')
    #endpts[:,0] = endpts[:,0] + 360.
    #R = 6.37*10**6 # radius of the Earth
    # convert endpts array to local cartesian co-ordinates:
    loccar = LocalCartesian(endpts)
    #loccar = pl.zeros_like(endpts)
    #endpts = m(endpts[:,0],endpts[:,1])
    #loccar[:,0] = endpts[:,0]; loccar[:,1] = endpts[:,1]
    
    nhat = pl.zeros([loccar.shape[0]-1,2])
    for point in range(loccar.shape[0]-1):
        if loccar[point+1,0] == loccar[point,0]:
            nhat[point] = pl.array([0,-1])
        elif loccar[point+1,1] == loccar[point,1]:
            nhat[point] = pl.array([-1,0])
        else:
            dx = loccar[point+1,0] - loccar[point,0]
            dy = loccar[point+1,1] - loccar[point,1]
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
    loccar[:,0] = R*pl.cos(coords[:,1])*pl.cos(coords[:,0])
    loccar[:,1] = R*pl.cos(coords[:,1])*pl.sin(coords[:,0]) # y = R*lat
    
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
    sheddir = '/home/np838619/Watershed/'
    
    #endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
    #endpts = pl.asarray(endlist[5:],dtype='float')
    
    rlslist = ReadTxtFile(sheddir + loc + '_traj_release.txt')
    rlspts = pl.asarray(rlslist[5:],dtype='float')
    
    seglab = [1] # first release point has to be on the first segment
    #segnos = pl.linspace(1,40,40)
    count = 1
    
    for rls in range(1,rlspts.shape[0]):
        if rlspts[rls,0] == rlspts[rls-1,0]:
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
    a = NearestIndex(lon,relpt[0]) # nearest longitude index
    b = NearestIndex(lat,relpt[1]) # nearest latitude index
    
    if lon[a] < relpt[0]: # nearest lon west of relpt
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
    
    dx = p2[0] - p1[0]; dy = p3[1] - p2[1]
    
    f1 = ((p2[0]-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
    f2 = ((p2[0]-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
    
    F = ((p3[1]-relpt[1])/dy)*f1 + ((relpt[1]-p2[1])/dy)*f2
    
    return F

def MidPts(rlspts):
    """
    """
    #rlspts[:,0] = rlspts[:,0] + 360.
    loccar = LocalCartesian(rlspts)
    #m = Basemap(llcrnrlon=-180,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=80.,
                #lat_ts=20.,resolution='l',projection='mill',suppress_ticks=True)
    #pl.zeros_like(rlspts)
    #loccar[:,0], loccar[:,1] = m(rlspts[:,0],rlspts[:,1])
    
    dl = pl.zeros([rlspts.shape[0]-1])
    for pt in range(rlspts.shape[0]-1):
        dl[pt] = pl.sqrt((loccar[pt+1,0] - loccar[pt,0])**2 +
                            (loccar[pt+1,1] - loccar[pt,1])**2)
    
    midpt = pl.zeros([rlspts.shape[0]])
    
    for pt in range(1,rlspts.shape[0]-1):
        midpt[pt] = 0.5*(dl[pt]+dl[pt-1])
    midpt[0] = dl[0]; midpt[-1] = dl[-1]
    
    return midpt

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

panfs = '/panfs/jasmin/era/era-in/netc/monthly_means/'

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
    filenames[year] = PrintFiles(path,'ggaw')
filenames = pl.sort(filenames,axis=1)

#path = panfs + '2007/'

#ncfile = Dataset(path+'ggaw200707.nc','r')
#lon = ncfile.variables['longitude'][:]
#lat = ncfile.variables['latitude'][:]
#wv_zon = ncfile.variables['TCUQ'][:]
#wv_mer = ncfile.variables['TCVQ'][:]
#ncfile.close()

# empty arrays for zonal & meridional water vapour fluxes:
wv_zon = pl.zeros([filenames.shape[0],filenames.shape[1],1,1,256,512]) # zonal component
wv_mer = pl.zeros_like(wv_zon) # meridional component

#loop over years:
for year in range(len(years)):
    #loop over filenames:
    for name in range(filenames.shape[1]):
        #load ncfile
        ncfile = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                str(year_input[year]) + '/' + str(filenames[year,name]),'r')
        #extract E & TP data
        wv_zon[year,name] = ncfile.variables['TCUQ'][:]
        wv_mer[year,name] = ncfile.variables['TCVQ'][:]
        ncfile.close()

# remove the 1D axes:
wv_zon = pl.squeeze(wv_zon); wv_mer = pl.squeeze(wv_mer)

# climatological monthly means:
#zon_mnths = pl.mean(wv_zon,axis=0); mer_mnths = pl.mean(wv_mer,axis=0)

# annual means:
zon_mean = pl.mean(wv_zon[:31],axis=(0,1)); mer_mean = pl.mean(wv_mer[:31],axis=(0,1))
#zon_mean = wv_zon[:31,:]; mer_mean = wv_mer[:31,:]
W = pl.sqrt(wv_zon**2 + wv_mer**2) # magnitude of water vapour flux

# make quiver plot:

latlon = Dataset('/panfs/jasmin/era/era-in/netc/monthly_means/' + \
                                    year_input[0] + '/' + filenames[0,0],'r')
lon = latlon.variables['longitude'][:]
lat = latlon.variables['latitude'][:]
latlon.close()

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

sheddir = '/home/np838619/Watershed/'
#endpts1 = pl.genfromtxt(sheddir+'NCA_clicks.txt',skip_header=5)
#endpts2 = pl.genfromtxt(sheddir+'SA_clicks.txt',skip_header=5)
#endpts = pl.zeros([len(endpts1)+len(endpts2),2])
#endpts[:len(endpts1)] = endpts1[:]; endpts[len(endpts1):] = endpts2[:]
AD_line = pl.genfromtxt(sheddir + 'AD_wshed.txt')
S16_wshed = pl.zeros_like(AD_line)
S16_wshed[:,0] = pl.flipud(AD_line[:,1]); S16_wshed[:,1] = pl.flipud(AD_line[:,0])
#f = pl.genfromtxt(sheddir+'NCA_traj_release.txt',skip_header=5)
#labs1, rlspts1 = TrajSegLabel('NCA'); labs2, rlspts2 = TrajSegLabel('SA')
#seglab = pl.zeros([len(labs1)+len(labs2)])
#seglab[:len(labs1)] = labs1[:]; seglab[len(labs1):] = labs2[:]
#rlspts = pl.zeros([len(rlspts1)+len(rlspts2),2])
#rlspts[:len(rlspts1)] = rlspts1[:]; rlspts[len(rlspts1):] = rlspts2[:]
#f = S16_wshed; f[:,0] = f[:,0] + 360.
#rlslabs = pl.zeros([len(seglab),3])
#rlslabs[:,0] = rlspts[:,0]; rlslabs[:,1] = rlspts[:,1]; rlslabs[:,2] = seglab[:]

#tst_zon = pl.mean(wv_zon[19:25,:],axis=(0,1)); tst_mer = pl.mean(wv_mer[19:25,:],axis=(0,1))


#rlspts[:,0] = rlspts[:,0]# + 360.


#bf = pl.zeros([182,2])
#for i in range(182):
    

#f = S16_wshed; f[:,0] = f[:,0] + 360.
#g = f[:,0] - 360.
#r,t = m2(rlspts[:,0],rlspts[:,1]); m2.scatter(r,t)


# calculate lengths between points
#rlspts[:,0] = rlspts[:,0] - 360.
#lengths = pl.zeros([len(rlspts)-1])
#for l in range(len(rlspts)-1):
    #lengths[l] = pl.sqrt((r[l+1]-r[l])**2 + (t[l+1]-t[l])**2)
#    lengths[l] = pl.sqrt((lc[l+1,0]-lc[l,0])**2 + (lc[l+1,1]-lc[l,1])**2)

#a = pl.where(lengths==0.); z = a[0]

#p = pl.zeros([len(r),2]); p[:,0] = r; p[:,1] = t
#==============================================================================
# nhat = pl.zeros([p.shape[0]-1,2])
# for point in range(p.shape[0]-1):
#     dx = p[point+1,0] - p[point,0]
#     dy = p[point+1,1] - p[point,1]
#     n = pl.array([-dy,dx])
#     nhat[point] = n/pl.norm(n)
#==============================================================================
nhat = NormalVector(S16_wshed); 
#x=pl.isnan(nhat); c=x[:,0]; c1=pl.where(c==True)

#seglab, rlspts = TrajSegLabel('NCA')

#lengths = list(lengths); bf = list(bf); nhat = list(nhat)
#L2 = []#; flux = [flux_uv[0]]; n2 = []#; LL = [f[0,1]]#; fx =[f[0]]
#repeat = []
#for r in range(1,rlspts.shape[0]):
#    if rlspts[r,0] == rlspts[r-1,0] and rlspts[r,1] == rlspts[r-1,1]:
#        repeat.append(r)
        #l2.append(lengths[l])
        #flux.append(flux_uv[l+1])
        #RP.append(rlslabs[l+1]) 
RP2 = []
#for r in range(rlspts.shape[0]):
#    if r not in repeat:
#        RP2.append(rlspts[r]); L2.append(seglab[r])
#rlspts = pl.asarray(RP2); labs = pl.asarray(L2)

midpt = MidPts(S16_wshed)

flux_uv = pl.zeros_like(S16_wshed)
for i in range(len(flux_uv)):
    a = NearestIndex(lon,S16_wshed[i,0]); b = NearestIndex(lat,S16_wshed[i,1])
    if lon[a] == S16_wshed[i,0]:
        flux_uv[i,0] = zon_mean[b,a]; flux_uv[i,1] = mer_mean[b,a]
    else:
        flux_uv[i,0] = BilinInterp(S16_wshed[i],lon,lat,zon_mean)
        flux_uv[i,1] = BilinInterp(S16_wshed[i],lon,lat,mer_mean)

midflux = pl.zeros_like(S16_wshed)
for f in range(1,midflux.shape[0]-1):
    midflux[f,0] = (flux_uv[f,0]+flux_uv[f-1,0])/2
    midflux[f,1] = (flux_uv[f,1]+flux_uv[f-1,1])/2
midflux[0,:] = flux_uv[0,:]
midflux[-1,:] = flux_uv[-1,:]

FdotN = pl.zeros([midflux.shape[0]-1])
for i in range(FdotN.shape[0]):
    #segno = labs[i] - 1
    FdotN[i] = pl.dot(midflux[i],nhat[i])

fluxes = FdotN[54:91]*midpt[54:91]/(10**9)
print pl.nansum(fluxes[:])
pl.close(); pl.close()

"""
# calculate flux across northern/southern boundaries of Atlantic:
m3 = Basemap(llcrnrlon=-180.0,llcrnrlat=-45.,urcrnrlon=25.,urcrnrlat=75.,\
                lat_ts=20.,resolution='l',projection='merc',suppress_ticks=True)
nth_lat = 70.; nth_lon1 = 360-136.; nth_lon2 = 18.
sth_lat = -35.; sth_lon1 = 360.-70; sth_lon2 = 20.
spacing = 0.7 # approximate spacing in degrees
nlat_ind = NearestIndex(lat,nrth_lat)
nlon_ind1 = NearestIndex(lon,nth_lon1); nlon_ind2 = NearestIndex(lon,nth_lon2)
slat_ind = NearestIndex(lat,sth_lat)
slon_ind1 = NearestIndex(lon,sth_lon1); slon_ind2 = NearestIndex(lon,sth_lon2)

q = lon[nlon_ind1:]; w = lon[:nlon_ind2]
e = pl.zeros([len(q)+len(w)]); e[:len(q)] = q; e[len(q):] = w
#e = lon[slon_ind1:slon_ind2+1]

m3.drawcoastlines()
for i in range(len(e)):
    if e[i] > 180.:
        #e[i] = e[i] - 360.
        a,b = m3(e[i]-360.,nth_lat,inverse=False)
    else:
        a,b = m3(e[i],nth_lat,inverse=False)
    m3.scatter(a,b)

inds = pl.zeros_like(e)
for i in range(len(e)):
    inds[i] = NearestIndex(lon,e[i])

flux_uv = pl.zeros([len(e),2])
for i in range(flux_uv.shape[0]):
    flux_uv[i,0] = tst_zon[nlat_ind,inds[i]]
    flux_uv[i,1] = tst_mer[nlat_ind,inds[i]]

coords = pl.zeros([len(e),2]); coords[:,0] = e[:]; coords[:,1] = nth_lat
loccar = LocalCartesian(coords)
lengths = pl.zeros([len(flux_uv)-1])
for point in range(len(coords)-1):
    lengths[point] = pl.sqrt((loccar[point+1,0]-loccar[point,0])**2 
                            + (loccar[point+1,1]-loccar[point,1])**2)
lmx=pl.where(lengths==lengths.max()); 
lengths[lmx[0][0]] = lengths[lmx[0][0]-1]

midpt = pl.zeros(lengths.shape[0]+1)
for point in range(1,midpt.shape[0]-1):
    midpt[point] = (lengths[point-1]+lengths[point])/2
midpt[0] = lengths[0]; midpt[-1] = lengths[-1]

normal = pl.array([0,1])
FdotN = pl.zeros([flux_uv.shape[0]])
for i in range(FdotN.shape[0]):
    FdotN[i] = pl.dot(flux_uv[i],normal)

print pl.sum(FdotN*midpt)/(10**9)"""