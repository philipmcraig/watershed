# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 10:58:36 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from scipy.interpolate import interp1d
from mpl_toolkits.basemap import Basemap

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
    
    if relpt[0] > lon[-1]:
        dx = (360. + lon[0]) - lon[-1]
    else:
        dx = p2[0] - p1[0]
    dy = p3[1] - p2[1]
    #dx = p2[0] - p1[0];dy = p3[1] - p2[1]
    
    if relpt[0] > lon[-1]:
        f1 = (((360+p2[0])-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
        f2 = (((360+p2[0])-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
    else:
        f1 = ((p2[0]-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
        f2 = ((p2[0]-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
    
    F = ((p3[1]-relpt[1])/dy)*f1 + ((relpt[1]-p2[1])/dy)*f2
    
    return F


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
    loccar[:,0] = R*pl.sin(coords[:,1])*pl.cos(coords[:,0])
    loccar[:,1] = R*pl.sin(coords[:,1])*pl.sin(coords[:,0]) # y = R*lat
    
    return loccar

def NormalVector(sheddir,loc):
    """Function to find the normal vectors to the line segments along a watershed.

	Args:
		sheddir (string): path to directory with end points file
		loc (string): stem indicating which watershed is being used

	Returns:
		n (array): unit normal vectors to line segments
    """
    endpts = pl.genfromtxt(sheddir + 'shed_defs/' + loc + '_clicks.txt',skip_header=5)
    #endpts = pl.asarray(endlist[5:],dtype='float')
    
    # convert endpts array to local cartesian co-ordinates:
    #endpts[:,0] = endpts[:,0] + 360.
    endpts[:,0] = endpts[:,0] + 360.
    endpts = pl.radians(endpts)#loccar = LocalCartesian(endpts)
    R = 6.37*10**6 # radius of the Earth
    
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
    
    rlspts = pl.genfromtxt(sheddir + loc + '_traj_release.txt',skip_header=5)
    rlspts[:,0] = rlspts[:,0] + 360.
    
    seglab = [1] # first release point has to be on the first segment
    #segnos = pl.linspace(1,40,40)
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

def VertInt(variable,surfp,erapres,pressure):
    """
    """
    g = 9.81
    erapres = erapres*100
    pressure = pressure*100
    dp = pl.zeros_like(pressure)
    if surfp == pressure[0]:
        dp[0] = 0.5*(pressure[0] - pressure[1])
        dp[-1] = 0.5*(pressure[-2] - pressure[-1])
        for p in range(1,len(dp)-1):
            dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    elif surfp > pressure[0]:
        dp[0] = surfp - 0.5*(pressure[0]+pressure[1])
        dp[-1] = 0.5*(pressure[-2] - pressure[-1])
        for p in range(1,len(dp)-1):
            dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    elif surfp < pressure[0]:
        b = pl.where(surfp > pressure); b = b[0][0]
        dp[:b] = 0.
        dp[b] = surfp - 0.5*(pressure[b]+pressure[b+1])#- surfp # + pressure[b]#
        dp[-1] = 0.5*(pressure[-2] - pressure[-1])
        for p in range(b+1,len(dp)-1):
            dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    #print dp
    X = pl.zeros([pressure.size])
                
    for p in range(len(dp)):
        a = pl.where(pressure[p]==erapres)
        lev = a[0][0]
        X[p] = variable[lev]*dp[p]#; print variable[lev]
    
    vint = (1/g)*pl.nansum(X,axis=0)
    
    return vint

def VIeta(variable,pressure,surfp,etalevs):
    """
    """
    g = 9.81
    pressure = pressure*100
    deta = etalevs[1] - etalevs[0] # should be constant
    dp = pl.zeros([pressure.shape[0]])
    
    for p in range(1,len(dp)-1):
        dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    dp[0] = 0.5*(surfp-pressure[0])
    dp[-1] =  0.5*(pressure[-2] - pressure[-1])
    
    X = pl.zeros_like(pressure)
    for lev in range(len(etalevs)):
        X[lev] = variable[lev]*dp[lev]
    
    vint = (1/g)*pl.sum(X,axis=0)
    
    return vint

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

def Add360(rlspts):
    """
    """
    for i in range(rlspts.shape[0]):
        if rlspts[i,0] < 0.:
            rlspts[i,0] = rlspts[i,0] + 360.
    
    return rlspts

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')


sheddir = '/home/np838619/Watershed/'
clusdir = '/net/glusterfs_essc/scenario/users/np838619/tjnew/TESTS/'
panfs = '/panfs/jasmin/era/era-in/netc/'
loc = 'Am'

ncfile = Dataset(sheddir+'ggap200707211200.nc','r')
#q = ncfile.variables['Q'][:]
#u = ncfile.variables['U'][:]
#v = ncfile.variables['V'][:]
pres = ncfile.variables['p'][:]
lon = ncfile.variables['longitude'][:]
lat = ncfile.variables['latitude'][:]
ncfile.close()
#q = pl.squeeze(q); u = pl.squeeze(u); v = pl.squeeze(v)
#
#nc2 = Dataset(clusdir+'ggas201408310000.nc','r')
#sp = nc2.variables['SP'][:]
#lon = nc2.variables['longitude'][:]
#lat = nc2.variables['latitude'][:]
#nc2.close()
#sp = pl.squeeze(sp)

years = pl.linspace(1979,2009,31); year_input = [str(int(i)) for i in years]
months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
filenames1 = PrintFiles(panfs+'ggap/2014/dec2014','ggap')
filenames2 = PrintFiles(panfs+'ggas/2014/dec2014','ggas')

q = pl.zeros([len(filenames1),37,256,512])
u = pl.zeros_like(q); v = pl.zeros_like(q)
surfp = pl.zeros([len(filenames2),256,512])


#"""for yr in range(len(years)):
#q_yr = pl.zeros([len(months),len(pres),len(lat),len(lon)])
#u_yr = pl.zeros_like(q_yr);  v_yr = pl.zeros_like(q_yr)
#sp_yr = pl.zeros([len(months),len(lat),len(lon)])
#for mt in range(len(months)):
#path1 = panfs + 'ggap' +'/' + '2007' + '/' + 'jul2007' + '/'
#filenames1 = PrintFiles(path1,'ggap')
#q_mt = pl.zeros([len(filenames1),len(pres),len(lat),len(lon)])
#u_mt = pl.zeros_like(q_mt); v_mt = pl.zeros_like(q_mt)
#qu_mt = pl.zeros_like(q_mt); qv_mt = pl.zeros_like(q_mt)
#path2 = panfs + 'ggas' +'/' + '2007' + '/' + 'jul2007' + '/'
#filenames2 = PrintFiles(path2,'ggas')
#sp_mt = pl.zeros([len(filenames2),len(lat),len(lon)])
for name in range(len(filenames1)):
    ncfile = Dataset(panfs+'ggap/2014/dec2014/'+filenames1[name],'r')
    #q_read = ncfile.variables['Q'][:]; q_read = pl.squeeze(q_read)
    q[name] = ncfile.variables['Q'][:]
    #u_read = ncfile.variables['U'][:]; u_read = pl.squeeze(u_read)
    u[name] = ncfile.variables['U'][:]
    #v_read = ncfile.variables['V'][:]; v_read = pl.squeeze(v_read)
    v[name] = ncfile.variables['V'][:]
    ncfile.close()
    other = Dataset(panfs+'ggas/2014/dec2014/'+filenames2[name],'r')
    surfp[name] = other.variables['SP'][:]#; sp_read = pl.squeeze(sp_read)
    other.close()
#    sp_mt[name,:,:] = sp_read; q_mt[name,:,:] = q_read
#    qu_mt[name,:,:] = q_read*u_read; qv_mt[name,:,:] = q_read*v_read
#    q_yr[mt] = pl.mean(q_mt,axis=0)
#    u_yr[mt] = pl.mean(u_mt,axis=0)
#    v_yr[mt] = pl.mean(v_mt,axis=0)
#    sp_yr[mt] = pl.mean(sp_mt,axis=0)
#    q[yr] = pl.mean(q_yr,axis=0); u[yr] = pl.mean(u_yr,axis=0); v[yr] = pl.mean(v_yr,axis=0)
#    surfp[yr] = pl.mean(sp_yr,axis=0)


q = pl.squeeze(q); u = pl.squeeze(u); v = pl.squeeze(v); surfp = pl.squeeze(surfp)

#q_mn = pl.mean(q_mt,axis=0); u_mn = pl.mean(u_mt,axis=0); v_mn = pl.mean(v_mt,axis=0) 
#sp_mn = pl.mean(sp_mt,axis=0)#"""

#q_mn = pl.squeeze(q_mn); u_mn = pl.squeeze(u_mn); v_mn = pl.squeeze(v_mn); sp_mn = pl.squeeze(sp_mn)

trajlevs = pl.linspace(950,150,17)

rlspts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_traj_release.txt',skip_header=5)
#rlspts[:,0] = rlspts[:,0] + 360.
rlspts = Add360(rlspts)

#rlspts = pl.zeros([len(rlspts1)+len(rlspts2),2])
#rlspts[:len(rlspts1)] = rlspts1[:]; rlspts[len(rlspts1):] = rlspts2[:]
#rlspts[:,0] = rlspts[:,0]
labs = TrajSegLabel(loc); L2 = []
#RP2 = [rlspts[0]]
repeat = []
for r in range(1,rlspts.shape[0]):
    if rlspts[r,0] == rlspts[r-1,0] and rlspts[r,1] == rlspts[r-1,1]:
        repeat.append(r)
        #L2.append(labs[r])

RP2 = []
for r in range(rlspts.shape[0]):
    if r not in repeat:
        RP2.append(rlspts[r]); L2.append(labs[r])
rlspts = pl.asarray(RP2); labs = pl.asarray(L2)

#qu = q_mn*u_mn; qv = q_mn*v_mn
#qn = q_mn
#for step in range(len(filenames1)):
#    for p in range(len(pres)):
#        for i in range(len(lat)):
#            for j in range(len(lon)):
#                if pres[p]*100 > sp_mn[i,j]:
#                    qn[p,i,j] = pl.float64('nan')
#                    qu[p,i,j] = pl.float64('nan')
#                    qv[p,i,j] = pl.float64('nan')
                
#QU = VertInt(qu,pres,trajlevs); QV = VertInt(qv,pres,trajlevs)

coeffs = pl.genfromtxt('/home/np838619/Trajectory/coeffs.txt')
a = coeffs[:,1]; b = coeffs[:,2]

pout = pl.genfromtxt('/home/np838619/Trajectory/pres_out.txt')
pout = pl.reshape(pout,(17,182))

qu = q*u; qv = q*v; qus = pl.zeros([rlspts.shape[0],pres.size]); qvs = pl.zeros_like(qus)
QU = pl.zeros([len(filenames1),lat.size,lon.size]); QV = pl.zeros_like(QU); Q = pl.zeros_like(QU)
for name in range(len(filenames1)):
    for i in range(lat.size):
        for j in range(lon.size):
            Q[name,i,j] = VertInt(q[name,:,i,j],surfp[-1,i,j],pres[:25],pres[:25])
            QU[name,i,j] = VertInt(qu[name,:,i,j],surfp[-1,i,j],pres[:25],pres[:25])
            #QU[i,j] = VIeta(qu[:,i,j],)
            QV[name,i,j] = VertInt(qv[name,:,i,j],surfp[-1,i,j],pres[:25],pres[:25])

QU_mn = pl.mean(QU,axis=0); QV_mn = pl.mean(QV,axis=0); Q_mn = pl.mean(Q,axis=0)

QU_shed = pl.zeros([rlspts.shape[0]]); QV_shed = pl.zeros([rlspts.shape[0]])
Q_shed = pl.zeros([rlspts.shape[0]])
for pt in range(rlspts.shape[0]):
    QU_shed[pt] = BilinInterp(rlspts[pt],lon,lat,QU_mn)
    QV_shed[pt] = BilinInterp(rlspts[pt],lon,lat,QV_mn)
    Q_shed[pt] = BilinInterp(rlspts[pt],lon,lat,Q_mn)
    #for lev in range(pres.size):
     #   qus[pt,lev] = BilinInterp(rlspts[pt],lon,lat,qu[lev])
     #   qvs[pt,lev] = BilinInterp(rlspts[pt],lon,lat,qv[lev])

QUV_shed = pl.zeros([rlspts.shape[0],2])
QUV_shed[:,0] = QU_shed; QUV_shed[:,1] = QV_shed

nhat = NormalVector(sheddir,loc)#; nhat2 = NormalVector(sheddir,'SA')
#nhat = pl.zeros([len(nhat1)+len(nhat2),2])
#nhat[:len(nhat1)] = nhat1[:]; nhat[len(nhat1):] = nhat2[:]
#; labs2 = TrajSegLabel('SA') + labs1.max()
#labs = pl.zeros([len(labs1)+len(labs2)])
#labs[:len(labs1)] = labs1[:]; labs[len(labs1):] = labs2[:]

midpt = MidPts(rlspts)

midflux = pl.zeros_like(rlspts)
for f in range(1,midflux.shape[0]-1):
    midflux[f,0] = (QUV_shed[f,0]+QUV_shed[f-1,0])/2
    midflux[f,1] = (QUV_shed[f,1]+QUV_shed[f-1,1])/2
midflux[0,:] = QUV_shed[0,:]
midflux[-1,:] = QUV_shed[-1,:]

NSEG = labs.max()

fluxes = pl.zeros([midflux.shape[0]])#[midflux.shape[0]] 116:209
magflx = pl.zeros_like(fluxes)

for pt in range(fluxes.shape[0]):#[midflux.shape[0]]
    segno = labs[pt] - 1
    fluxes[pt] = pl.dot(midflux[pt],nhat[segno])*midpt[pt]
    magflx[pt] = pl.dot(pl.absolute(midflux[pt]),pl.absolute(nhat[segno]))*midpt[pt]

print pl.sum(fluxes)/(10**9)
pl.plot(fluxes/(10**9))



"""newnc = Dataset('sp_annmean.nc',mode='w',format='NETCDF4')
lat_dim = newnc.createDimension('lat', 256)
lon_dim = newnc.createDimension('lon', 512)
#pres_dim = newnc.createDimension('pres', 37)
lat_in = newnc.createVariable('lat', pl.float32, ('lat',))
lat_in.units = 'degrees_north'
lat_in.long_name = 'latitude'
lon_in = newnc.createVariable('lon', pl.float32, ('lon',))
lon_in.units = 'degrees_east'
lon_in.long_name = 'longitude'
#prss = newnc.createVariable('pres', pl.float64, ('pres',))
#prss.units = 'Hectoascals'
#prss.long_name = 'pressure'

lat_in[:] = lat # straight from ERA-Interim nc file
lon_in[:] = lon # straight from ERA-Interim nc file
#prss[:] = pres

s_pres = newnc.createVariable('surfp',pl.float64,('lat','lon'))
s_pres.units = 'Pa'
s_pres.standard_name = 'annual mean surface pressure'
s_pres[:,:] = sp_mn

newnc.close()"""