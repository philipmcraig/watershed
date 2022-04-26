# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 14:39:31 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs

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

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
exec(open('/home/np838619/Trajectory/trajfuncs.py').read())

pl.close('all')

panfs = '/panfs/jasmin/era/era-in/netc/monthly_means/'
sheddir = '/home/np838619/Watershed/'

loc = 'ArI'
rlspts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_traj_release_new.txt',skip_header=5)
endpts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_clicks.txt',skip_header=5)
labs = TrajSegLabel(loc)
l2 = pl.zeros([labs[0].size,3]); l2[:,0] = labs[0]; l2[:,1:] = labs[1]
labs = l2.copy()
nhat = NormalVector(sheddir+'shed_defs/',loc)
labs = RR2(labs,loc)

# list of years
years = pl.linspace(1979,2014,36)
# need next part to get rid of decimal point
year_input = [str(int(i)) for i in years]

# empty array for filenames
filenames1 = pl.zeros([len(years),12],dtype='S13')
filenames2 = pl.zeros_like(filenames1)

for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year] + '/'
    #filenames[year] = PrintFiles(path,type)
    filenames1[year] = PrintFiles(path,'ggap')
    filenames2[year] = PrintFiles(path,'ggas')
    filenames1[year] = pl.sort(filenames1[year])
    filenames2[year] = pl.sort(filenames2[year])

#path = panfs + '2007/'

ncfile = Dataset(sheddir+'ggap201107211200.nc','r')
lon = ncfile.variables['longitude'][:]
lat = ncfile.variables['latitude'][:]
pres = ncfile.variables['p'][:]
#zon_ann = ncfile.variables['tcuq'][:]
#mer_ann = ncfile.variables['tcvq'][:]
ncfile.close()

# empty arrays for zonal & meridional water vapour fluxes:
#wv_zon = pl.zeros([filenames1.shape[0],filenames1.shape[1],1,1,256,512]) # zonal component
#wv_mer = pl.zeros_like(wv_zon) # meridional component

v = pl.zeros([filenames1.shape[0],filenames1.shape[1],1,len(pres),256,512])
q = pl.zeros_like(v)
surfp = pl.zeros([filenames2.shape[0],filenames2.shape[1],1,256,512])

#loop over years:q
for year in range(len(years)):
#    loop over filenames:
    for name in range(filenames1.shape[1]):
        #load ncfile
        nc1 = Dataset(panfs + year_input[year] + '/' + str(filenames1[year,name]),'r')
        #extract E & TP data
        v[year,name] = nc1.variables['V'][:]
        q[year,name] = nc1.variables['Q'][:]
        nc1.close()
        nc2 = Dataset(panfs+year_input[year]+'/'+str(filenames2[year,name]),'r')
        surfp[year,name] = nc2.variables['SP'][:]
        nc2.close()

# remove the 1D axes:
#wv_zon = pl.squeeze(wv_zon); wv_mer = pl.squeeze(wv_mer)
v = pl.squeeze(v); q = pl.squeeze(q); surfp = pl.squeeze(surfp)
#zon_flt = pl.reshape(wv_zon,(36*12,256,512)); mer_flt = pl.reshape(wv_mer,(36*12,256,512))
#zon_ann = pl.mean(wv_zon,axis=(0,1)); mer_ann = pl.mean(wv_mer,axis=(0,1))


#U_mn = pl.mean(U,axis=(0,1)); V_mn = pl.mean(V,axis=(0,1)); Q_mn = pl.mean(Q,axis=(0,1))
#UQ_mn = pl.mean(U*Q)

qv = q*v
#QU = pl.zeros([filenames1.shape[0],filenames1.shape[1],lat.size,lon.size])

v_mn = pl.mean(v,axis=(0,1)); v_pt = v - v_mn
del v, v_mn
q_mn = pl.mean(q,axis=(0,1)); q_pt = q - q_mn
del q, q_mn
sp_mn = pl.mean(surfp,axis=(0,1)); del surfp
#u_pt = u - u_mn; q_pt = q - q_mn
#uq_mn = pl.mean(qu,axis=(0,1))
#uq_pt = qu - uq_mn

#QU_prod = pl.zeros([lat.size,lon.size])
QV_pert = pl.zeros([lat.size,lon.size])


print 'Starting loops now!'
for i in range(lat.size):
    for j in range(lon.size):
        #QU_prod[i,j] = VertInt(q_mn[:,i,j]*u_mn[:,i,j],sp_mn[i,j],pres,pres)
        QV_pert[i,j] = VertInt(pl.mean(q_pt[:,:,:,i,j]*v_pt[:,:,:,i,j],axis=(0,1)),sp_mn[i,j],pres,pres)
        #for year in range(filenames2.shape[0]):
        #    for month in range(filenames2.shape[1]):
        #        QU[year,month,i,j] = VertInt(qu[year,month,:,i,j],surfp[year,month,i,j],pres[:],pres[:])

# Save to netCDF

newnc = Dataset(sheddir+'reynolds/reyn_av_quv_7914.nc','a')

QV_in = newnc.createVariable('Vertically-integrated meridional moisture flux (eddies)',
                                             pl.float64,('lat','lon'))
QV_in.units = 'kg/m/s'
QV_in.standard_name = 'meridional moisture flux (eddies)'
QV_in[:,:] = QV_pert

newnc.close()

#QU_mn = pl.mean(QU,axis=(0,1))


#rlspts = Add360(rlspts)
#midpt = MidPts(rlspts)
#midflux = MidFlux(rlspts,lon,lat,pl.mean(QU_mn,axis=0),pl.mean(QV_mn,axis=0))
#FdotN = NormalFlux(midflux,labs,nhat)
#fluxes = FdotN[:]*midpt[:]/(10**9); print pl.sum(fluxes)
#
## mean flow:
#mfx_mf = MidFlux(rlspts,lon,lat,QU_prod,QV_prod)
#FdN_mf = NormalFlux(mfx_mf,labs,nhat)
#flx_mf = FdN_mf[:]*midpt[:]/(10**9); print pl.sum(flx_mf)
#
## eddies:
#mfx_ed = MidFlux(rlspts,lon,lat,QU_pert,QV_pert)
#FdN_ed = NormalFlux(mfx_ed,labs,nhat)
#flx_ed = FdN_ed[:]*midpt[:]/(10**9); print pl.sum(flx_ed)