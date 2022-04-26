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
from scipy.stats import pearsonr, linregress

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
    
    rlspts = pl.genfromtxt(sheddir+loc+'_traj_release.txt',skip_header=5)
    
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
    
    return seglab#, rlspts

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
    
    return labs#,rlspts, labs

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
        segno = labs[i] - 1
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
    
    dsig = 2*pl.arcsin(pl.sqrt((pl.sin(dphi/2)**2 + pl.cos(phi1)*pl.cos(phi2)*pl.sin(dlam/2)**2)))
    
    R = 6.37e6
    
    d = R*dsig
    
    return d
    
def ShedFluxes(sheddir,loc,lon,lat,zon,mer):
    """
    """
    endpts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_clicks.txt',skip_header=5)
    labs = TrajSegLabel(loc)
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
    
    #print pl.sum(distances)
    flux_per_len = pl.sum(fluxes)/pl.sum(distances)
    
    return flux_per_len

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
panfs = '/panfs/jasmin/era/era-in/netc/monthly_means/'
sheddir = '/home/np838619/Watershed/'

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

ncfile = Dataset(sheddir+'wvfluxes_7914.nc','r')
lon = ncfile.variables['lon'][:]
lat = ncfile.variables['lat'][:]
#zon_ann = ncfile.variables['tcuq'][:]
#mer_ann = ncfile.variables['tcvq'][:]
ncfile.close()

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
zon_mnths = pl.mean(wv_zon,axis=0); mer_mnths = pl.mean(wv_mer,axis=0)

# annual means:
zon_mean = pl.mean(wv_zon,axis=0); mer_mean = pl.mean(wv_mer,axis=0)
zon_ann = pl.mean(wv_zon,axis=1); mer_ann = pl.mean(wv_mer,axis=1)
zon_mean = wv_zon[28,6]; mer_mean = wv_mer[28,6]
#W = pl.sqrt(wv_zon**2 + wv_mer**2) # magnitude of water vapour flux


# can't remember why these are here, need to check
#lon_0 = lon.mean()
#lat_0 = lat.mean()

#zon_mn, newlons = shiftgrid(180.0, wv_zon, lon, start=False)
#mer_mn, newlons = shiftgrid(180.0, wv_mer, lon, start=False)
#Ws, newlons = shiftgrid(180.0, W, lon, start=False)

#zon_shft = pl.zeros([years.size-1,wv_zon.shape[1],lat.size,lon.size])
#mer_shft = pl.zeros_like(zon_shft)
#zon_shft[:,:7,:,:] = wv_zon[1:,5:,:,:]; zon_shft[:,7:,:,:] = wv_zon[1:,:5,:,:]
#mer_shft[:,:7,:,:] = wv_mer[1:,5:,:,:]; mer_shft[:,7:,:,:] = wv_mer[1:,:5,:,:]
#zon_ann = pl.mean(zon_shft,axis=1); mer_ann = pl.mean(mer_shft,axis=1)

#repeat = []
#for r in range(1,rlspts.shape[0]):
#    if rlspts[r,0] == rlspts[r-1,0] and rlspts[r,1] == rlspts[r-1,1]:
#        repeat.append(r)
#        #l2.append(lengths[l])
#        #flux.append(flux_uv[l+1])
#        #RP.append(rlslabs[l+1]) 
#RP2 = []
#for r in range(rlspts.shape[0]):
#    if r not in repeat:
#        RP2.append(rlspts[r]); L2.append(labs[r])
#rlspts = pl.asarray(RP2); labs = pl.asarray(L2)

f1_yrs = []; f2_yrs = []; f3_yrs = [];
f4_yrs = []; f5_yrs = []; f6_yrs =[]
f7_yrs = []; f8_yrs = []; f9_yrs = []; f10_yrs = []
#for year in range(years.size):
for yr in range(filenames.shape[0]):
    f1_yrs.append(ShedFluxes(sheddir,'Am',lon,lat,zon_ann[yr],mer_ann[yr]))
    f2_yrs.append(ShedFluxes(sheddir,'AfMe',lon,lat,zon_ann[yr],mer_ann[yr]))
    f3_yrs.append(ShedFluxes(sheddir,'EAA',lon,lat,zon_ann[yr],mer_ann[yr]))
    f4_yrs.append(ShedFluxes(sheddir,'ArA',lon,lat,zon_ann[yr],mer_ann[yr]))
    f5_yrs.append(ShedFluxes(sheddir,'ArI',lon,lat,zon_ann[yr],mer_ann[yr]))
    f6_yrs.append(ShedFluxes(sheddir,'ArP',lon,lat,zon_ann[yr],mer_ann[yr]))
    f7_yrs.append(ShedFluxes(sheddir,'SOA',lon,lat,zon_ann[yr],mer_ann[yr]))
    f8_yrs.append(ShedFluxes(sheddir,'SOI',lon,lat,zon_ann[yr],mer_ann[yr]))
    f9_yrs.append(ShedFluxes(sheddir,'SOP',lon,lat,zon_ann[yr],mer_ann[yr]))
    f10_yrs.append(ShedFluxes(sheddir,'NAs',lon,lat,zon_ann[yr],mer_ann[yr]))
f1_yrs = pl.asarray(f1_yrs); f2_yrs = pl.asarray(f2_yrs); f3_yrs = pl.asarray(f3_yrs)
f4_yrs = pl.asarray(f4_yrs); f5_yrs = pl.asarray(f5_yrs); f6_yrs = pl.asarray(f6_yrs)
f7_yrs = pl.asarray(f7_yrs); f8_yrs = pl.asarray(f8_yrs); f9_yrs = pl.asarray(f9_yrs)
f10_yrs = pl.asarray(f10_yrs)

f1_tots = pl.sum(f1_yrs,axis=1); f2_tots = pl.sum(f2_yrs,axis=1)
f3_tots = pl.sum(f3_yrs,axis=1); f4_tots = pl.sum(f4_yrs,axis=1)
f5_tots = pl.sum(f5_yrs,axis=1); f6_tots = pl.sum(f6_yrs,axis=1)
f7_tots = pl.sum(f7_yrs,axis=1); f8_tots = pl.sum(f8_yrs,axis=1)
f9_tots = pl.sum(f9_yrs,axis=1); f10_tots = pl.sum(f10_yrs,axis=1)

#fluxes1 = ShedFluxes(sheddir,'Am',lon,lat,zon_mean,mer_mean)
#fluxes2 = ShedFluxes(sheddir,'AfMe',lon,lat,zon_mean,mer_mean)
#fluxes3 = ShedFluxes(sheddir,'EAA',lon,lat,zon_mean,mer_mean)
#fluxes4 = ShedFluxes(sheddir,'SO',lon,lat,zon_mean,mer_mean)
#fluxes5 = ShedFluxes(sheddir,'Ar',lon,lat,zon_mean,mer_mean)
#fluxes6 = ShedFluxes(sheddir,'NAs',lon,lat,zon_mean,mer_mean)

atl = f1_tots - f2_tots - f4_tots + f7_tots
ind = f2_tots - f3_tots + f8_tots - f5_tots
pac = f3_tots - f1_tots - f6_tots + f9_tots
arc = f4_tots + f5_tots + f6_tots
sou = -1*(f7_tots + f8_tots + f9_tots)

#print 'Americas = ', pl.sum(fluxes1[:]), ' Sv'
#print 'Africa & Middle East = ', pl.sum(fluxes2[:]), ' Sv'
#print 'East Asia & Australia = ', pl.sum(fluxes3[:]), ' Sv'
#print 'Southern Ocean = ', pl.sum(fluxes4[:]), ' Sv'
#print 'Arctic = ', pl.sum(fluxes5[:]), ' Sv'
#print 'North Asia = ', pl.sum(fluxes6[:]), ' Sv'

#nas = f10_tots # north Asia
#wru = pl.sum(f4_yrs[:,137:],axis=1) # west Russia
#ari = f5_tots # Arctic drainage boundary, Indian sector
#chi = pl.sum(f6_yrs[:,:72],axis=1) # China
#divq = -nas+wru+ari+chi

x = pl.linspace(0,11,12); c = pl.linspace(0,35,36)

ax, fig = pl.subplots(5,1,figsize=(9,11))

ax1 = pl.subplot(511)
pl.plot(c,atl,linewidth=2,color='r')
pl.xlim(0,35)
f1_reg = linregress(c,atl)
ax1.plot(c,f1_reg[0]*c+f1_reg[1],lw=2,ls='--',color='r')
ax1.xaxis.set_major_formatter(pl.NullFormatter())
ax1.annotate('(a) Atlantic',(0.01,0.83),xycoords='axes fraction',size=20)
pl.ylim(-0.6,-0.2); ax1.set_yticks(pl.linspace(-0.6,-0.2,5))
#ax1.set_yticklabels([-0.3,' ',-0.2,' ',-0.1,' ',0.0])#,'',-0.2])
pl.ylabel('Sv',fontsize=20,labelpad=3); ax1.grid()
ax1.set_xticks(pl.linspace(0,35,8))
ax1.tick_params(axis='y', which='major', labelsize=14)

ax2 = pl.subplot(512)
pl.plot(c,ind,linewidth=2,color='b')
pl.xlim(0,35)
f2_reg = linregress(c,ind)
ax2.plot(c,f2_reg[0]*c+f2_reg[1],lw=2,ls='--',color='b')
ax2.xaxis.set_major_formatter(pl.NullFormatter())
ax2.annotate('(b) Indian',(0.01,0.83),xycoords='axes fraction',size=20)
pl.ylim(-0.8,-0.4); ax2.set_yticks(pl.linspace(-0.8,-0.4,5))
#ax2.set_yticklabels([-0.3,' ',-0.2,' ',-0.1,' ',0.0])
ax2.grid(); pl.ylabel('Sv',fontsize=20,labelpad=3); #; pl.xlim(0,11)
ax2.set_xticks(pl.linspace(0,35,8))
ax2.tick_params(axis='y', which='major', labelsize=14)

ax3 = pl.subplot(513)
pl.plot(c,pac,linewidth=2,color='g')
pl.xlim(0,35)
f3_reg = linregress(c,pac)
ax3.plot(c,f3_reg[0]*c+f3_reg[1],lw=2,ls='--',color='g')
ax3.xaxis.set_major_formatter(pl.NullFormatter())
ax3.annotate('(c) Pacific',(0.01,0.83),xycoords='axes fraction',size=20)
pl.ylim(-0.2,0.2); ax3.set_yticks(pl.linspace(-0.2,0.2,5))
#ax3.set_yticklabels([0.0,' ',0.1,' ',0.2,'',0.3])
ax3.grid(); pl.ylabel('Sv',fontsize=20,labelpad=3); #; pl.xlim(0,11)
ax3.set_xticks(pl.linspace(0,35,8))
ax3.tick_params(axis='y', which='major', labelsize=14)

ax4 = pl.subplot(514)
pl.plot(c,arc,linewidth=2,color='deeppink',label='$P-E$')
pl.xlim(0,35)
f4_reg = linregress(c,arc)
ax4.plot(c,f4_reg[0]*c+f4_reg[1],lw=2,ls='--',color='deeppink',label='trend')
ax4.xaxis.set_major_formatter(pl.NullFormatter())
ax4.annotate('(d) Arctic',(0.01,0.83),xycoords='axes fraction',size=20)
pl.ylim(0.0,0.4); ax4.set_yticks(pl.linspace(0.0,0.4,5))
#ax4.set_yticklabels([-0.1,' ',0.0,' ',0.1,' ',0.2])
pl.ylabel('Sv',fontsize=20,labelpad=10);
ax4.grid()#; pl.xlim(0,11)
#ax4.set_xticklabels(pl.linspace(1979,2014,8).astype('int'),fontsize=14,rotation=45)
ax4.set_xticks(pl.linspace(0,35,8))
ax4.tick_params(axis='y', which='major', labelsize=14)
ax4.legend(loc=1,fontsize=16,ncol=2,columnspacing=0.4)

ax5 = pl.subplot(515)
pl.plot(c,sou,linewidth=2,color='darkgoldenrod')
pl.xlim(0,35)
f5_reg = linregress(c,sou)
ax5.plot(c,f5_reg[0]*c+f5_reg[1],lw=2,ls='--',color='darkgoldenrod')
#ax5.xaxis.set_major_formatter(pl.NullFormatter())
ax5.annotate('(e) Southern Ocean',(0.01,0.83),xycoords='axes fraction',size=20)
#ax5.set_xticks(x)#;
pl.ylim(0.8,1.2); ax5.set_yticks(pl.linspace(0.8,1.2,5))
#ax5.set_yticklabels([-0.1,' ',0.0,' ',0.1,' ',0.2])
ax5.set_xticks(pl.linspace(0,35,8))
ax5.set_xticklabels(pl.linspace(1979,2014,8).astype('int'),fontsize=18)
ax5.tick_params(axis='x', which='major', pad=10)
ax5.grid(); pl.ylabel('Sv',fontsize=20,labelpad=8); 
ax5.tick_params(axis='y', which='major', labelsize=14)

#ax6 = pl.subplot(336)
#pl.plot(f6_tots,linewidth=2,color='r')
#pl.xlim(0,35)
#f6_reg = linregress(c,f6_tots)
#ax6.plot(c,f6_reg[0]*c+f6_reg[1],lw=2,ls='--',color='r')
#ax6.xaxis.set_major_formatter(pl.NullFormatter())
#ax6.annotate('(f) Atlantic Pacific',(0.01,0.85),xycoords='axes fraction',size=20)
##ax5.set_xticks(x)#;
#pl.ylim(-0.1,0.2); ax6.set_yticks(pl.linspace(-0.1,0.2,7))
#ax6.set_yticklabels([-0.1,'',0.,'',0.1,'',0.2])
#ax6.set_xticks(pl.linspace(0,35,8))
##ax6.set_xticklabels(pl.linspace(1979,2014,8).astype('int'),fontsize=14,rotation=45)
#ax6.tick_params(axis='x', which='major', pad=15)
#ax6.grid()#pl.ylabel('Sv',fontsize=20,labelpad=20);
#ax6.tick_params(axis='y', which='major', labelsize=14)
#
#ax7 = pl.subplot(337)
#pl.plot(f7_tots,linewidth=2,color='darkgoldenrod')
#pl.xlim(0,35)
#f7_reg = linregress(c,f7_tots)
#ax7.plot(c,f7_reg[0]*c+f7_reg[1],lw=2,ls='--',color='darkgoldenrod')
#ax7.annotate('(g) Southern Atlantic',(0.01,0.85),xycoords='axes fraction',size=20)
##ax5.set_xticks(x)#;
#pl.ylim(-0.5,-0.2); ax7.set_yticks(pl.linspace(-0.5,-0.2,7))
#ax7.set_yticklabels([-0.5,'',-0.4,'',-0.3,'',-0.2])
#ax7.set_xticks(pl.linspace(0,35,8))
#ax7.set_xticklabels(pl.linspace(1979,2014,8).astype('int'),fontsize=14,rotation=45)
#ax7.tick_params(axis='x', which='major', pad=2)
#pl.ylabel('Sv',fontsize=20,labelpad=4); ax7.grid()
#ax7.tick_params(axis='y', which='major', labelsize=14)
#
#ax8 = pl.subplot(338)
#pl.plot(f8_tots,linewidth=2,color='darkgoldenrod')
#pl.xlim(0,35)
#f8_reg = linregress(c,f8_tots)
#ax8.plot(c,f8_reg[0]*c+f8_reg[1],lw=2,ls='--',color='darkgoldenrod')
#ax8.annotate('(h) Southern Indian',(0.01,0.85),xycoords='axes fraction',size=20)
##ax5.set_xticks(x)#;
#pl.ylim(-0.5,-0.2); ax8.set_yticks(pl.linspace(-0.5,-0.2,7))
#ax8.set_yticklabels([-0.5,'',-0.4,'',-0.3,'',-0.2])
#ax8.set_xticks(pl.linspace(0,35,8))
#ax8.set_xticklabels(pl.linspace(1979,2014,8).astype('int'),fontsize=14,rotation=45)
#ax8.tick_params(axis='x', which='major', pad=2)
#ax8.grid()#pl.ylabel('Sv',fontsize=20,labelpad=20); 
#ax8.tick_params(axis='y', which='major', labelsize=14)
#
#ax9 = pl.subplot(339)
#pl.plot(f9_tots,linewidth=2,color='darkgoldenrod')
#pl.xlim(0,35)
#f9_reg = linregress(c,f9_tots)
#ax9.plot(c,f9_reg[0]*c+f9_reg[1],lw=2,ls='--',color='darkgoldenrod')
#ax9.annotate('(i) Southern Pacific',(0.01,0.85),xycoords='axes fraction',size=20)
##ax5.set_xticks(x)#;
#pl.ylim(-0.5,-0.2); ax9.set_yticks(pl.linspace(-0.5,-0.2,7))
#ax9.set_yticklabels([-0.5,'',-0.4,'',-0.3,'',-0.2])
#ax9.set_xticks(pl.linspace(0,35,8))
#ax9.set_xticklabels(pl.linspace(1979,2014,8).astype('int'),fontsize=14,rotation=45)
#ax9.tick_params(axis='x', which='major', pad=2)
#ax9.grid()#pl.ylabel('Sv',fontsize=20,labelpad=20);
#ax9.tick_params(axis='y', which='major', labelsize=14)

#pl.suptitle('Ocean Catchments $-$div$Q$',fontsize=24)
#pl.subplots_adjust(left=0.08,right=0.98,top=0.96,bottom=0.06,hspace=0.13)
pl.tight_layout()
#ax6 = pl.subplot(516)
#pl.plot(pl.linspace(1979,2014,36),f6_tots)
#pl.xlim(1979,2014)

f1_std = pl.std(f1_tots); print round(f1_std,3), ' Sv', (f1_std/pl.mean(f1_tots))*100, ' %'
f2_std = pl.std(f2_tots); print round(f2_std,3), ' Sv', (f2_std/pl.mean(f2_tots))*100, ' %'
f3_std = pl.std(f3_tots); print round(f3_std,3), ' Sv', (f3_std/pl.mean(f3_tots))*100, ' %'
f4_std = pl.std(f4_tots); print round(f4_std,3), ' Sv', (f4_std/pl.mean(f4_tots))*100, ' %'
f5_std = pl.std(f5_tots); print round(f5_std,3), ' Sv', (f5_std/pl.mean(f5_tots))*100, ' %'

soa_ave = pl.zeros([36]); soi_ave = pl.zeros([36]); sop_ave = pl.zeros([36])
rls_satl = pl.genfromtxt(sheddir+'shed_defs/'+'SOA'+'_traj_release_new.txt',skip_header=5)
rls_sind = pl.genfromtxt(sheddir+'shed_defs/'+'SOI'+'_traj_release_new.txt',skip_header=5)
rls_spac = pl.genfromtxt(sheddir+'shed_defs/'+'SOP'+'_traj_release_new.txt',skip_header=5)
rls_satl = Add360(rls_satl); rls_sind = Add360(rls_sind); rls_spac = Add360(rls_spac)

for i in range(36):
    soa_ave[i] = FluxPerLength(rls_satl[95:],f7_yrs[i,95:])
    soi_ave[i] = FluxPerLength(rls_sind[:113],f8_yrs[i,:113])
    sop_ave[i] = FluxPerLength(rls_spac[30:180],f8_yrs[i,30:180])

#fig,ax = pl.subplots(1,2,figsize=(13,6))
#
#ax1 = pl.subplot(121)
#ax1.plot(soa_ave*(10**8),label='Atlantic',lw=2.,color='b')
#ax1.plot(soi_ave*(10**8),label='Indian',lw=2.,color='r')
#ax1.plot(sop_ave*(10**8),label='Pacific',lw=2.,color='g')
##pl.xlim(0,11)
#pl.ylabel('10$^8$ kg m$^{-1}$ s$^{-1}$',fontsize=22,labelpad=-1.5)
#ax1.set_yticks(pl.linspace(-5,0,11))
#pl.tick_params(axis='y', which='major', labelsize=14)
##ax1.set_xticks(pl.linspace(0,11,12))
##ax1.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun',
##                    'Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=16)
#ax1.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=16)
#ax1.grid(axis='y')
#ax1.tick_params(axis='x',direction='out', pad=10)
#ax1.annotate('Flux per\nunit length',xy=(2,-0.75),fontsize=22,
#             bbox=dict(facecolor='w', edgecolor='k'))
#
#ax2 = pl.subplot(122)
#ax2.plot(f7_tots,label='Atlantic',lw=2.,color='b')
#ax2.plot(f8_tots,label='Indian',lw=2.,color='r')
#ax2.plot(f9_tots,label='Pacific',lw=2.,color='g')
##pl.xlim(0,11); pl.ylim(-0.5,0)
#pl.ylabel('Sv',fontsize=22,labelpad=-1.5); pl.legend(fontsize=18)
#ax2.set_yticks(pl.linspace(-0.5,0,11))
#pl.tick_params(axis='y', which='major', labelsize=14)
##ax2.set_xticks(pl.linspace(0,11,12))
##ax2.set_xticklabels(['Jan','Feb','Mar','Apr','May','Jun',
##                    'Jul','Aug','Sep','Oct','Nov','Dec'],fontsize=16)
#ax2.set_xticklabels(pl.linspace(1979,2014,8).astype(int),fontsize=16)
#ax2.grid(axis='y')
#ax2.tick_params(axis='x',direction='out', pad=10)
#ax2.annotate('Integrated flux',xy=(2,-0.04),fontsize=22,
#             bbox=dict(facecolor='w', edgecolor='k'))
#
#pl.subplots_adjust(left=0.08,right=0.96)

"""newnc = Dataset(sheddir+'wvfluxes_7914.nc','w')

lat_dim = newnc.createDimension('lat', 256)
lon_dim = newnc.createDimension('lon', 512)
lat_in = newnc.createVariable('lat', pl.float32, ('lat',))
lat_in.units = 'degrees_north'
lat_in.long_name = 'latitude'
lon_in = newnc.createVariable('lon', pl.float32, ('lon',))
lon_in.units = 'degrees_east'
lon_in.long_name = 'longitude'

time_dim = newnc.createDimension('time',36)
time = newnc.createVariable('time', pl.float64, ('time',))
time.units = 'years'
time.long_name = 'time'

tcuq = newnc.createVariable('tcuq',pl.float64,('time','lat','lon'))
tcuq.units = 'kg m**-1 s**-1'
tcuq.standard_name = 'vertical integral of eastward water vapour flux'
lat_in[:] = lat # straight from ERA-Interim nc file
lon_in[:] = lon # straight from ERA-Interim nc file
tcuq[:,:,:] = zon_ann

tcvq = newnc.createVariable('tcvq',pl.float64,('time','lat','lon'))
tcvq.units = 'kg m**-1 s**-1'
tcvq.standard_name = 'vertical integral of northward water vapour flux'
tcvq[:,:,:] = mer_ann

newnc.close()"""