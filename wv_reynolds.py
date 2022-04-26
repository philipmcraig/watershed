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

def SeasonalCyc(variable_array):
    """Function to calculate the climatological seasonal means (DJF,MAM,JJA,SON)
    of a variable.
    
    Args:
        variable_array (array): variable of which the climatological seasonal 
                                means are required
    
    Returns:
        variable_seasons (array): climatological seasonal means of the input 
                                  variable
    """
    # Calculate the mean of each trio of months for each year, missing out the 
    # first year as there is no data for December the year before the data starts
    MAM = variable_array[2:5] # March, April, May
    JJA = variable_array[5:8] # June, July, August
    SON = variable_array[8:11] # September, October, November
    DJF = pl.zeros_like(MAM) # December, January, February
    DJF[0] = variable_array[-1]; DJF[1:] = variable_array[:2]
    
    # Calculate the climatological mean of each season:  
    MAM_mn = pl.mean(MAM,axis=0) # only need axis=0 for profiles!!!!
    JJA_mn = pl.mean(JJA,axis=0)
    SON_mn = pl.mean(SON,axis=0)
    DJF_mn = pl.mean(DJF,axis=0)
    
    # Stick all seasons in one array before outputting:
    variable_seasons = pl.array([DJF_mn,MAM_mn,JJA_mn,SON_mn])
    
    return variable_seasons

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
exec(open('/home/np838619/Trajectory/trajfuncs.py').read())

pl.close('all')

panfs = '/panfs/jasmin/era/era-in/netc/monthly_means/'
sheddir = '/home/np838619/Watershed/'

# list of years
years = pl.linspace(1979,2014,36)
# need next part to get rid of decimal point
year_input = [str(int(i)) for i in years]

# empty array for filenames
filenames1 = pl.zeros([len(years),12],dtype='S13')
filenames2 = pl.zeros_like(filenames1)
filenames3 = pl.zeros_like(filenames1)
filenames4 = pl.zeros_like(filenames1)

for year in range(len(years)):
    #path = path to mm folder + year
    path = '/panfs/jasmin/era/era-in/netc/monthly_means/' + year_input[year] + '/'
    #filenames[year] = PrintFiles(path,type)
    filenames1[year] = PrintFiles(path,'ggap')
    filenames2[year] = PrintFiles(path,'ggas')
    filenames3[year] = PrintFiles(path,'ggaw')
    filenames1[year] = pl.sort(filenames1[year])
    filenames2[year] = pl.sort(filenames2[year])
    filenames3[year] = pl.sort(filenames3[year])

#path = panfs + '2007/'

ncfile = Dataset(sheddir+'ggap201107211200.nc','r')
lon = ncfile.variables['longitude'][:]
lat = ncfile.variables['latitude'][:]
pres = ncfile.variables['p'][:]
#zon_ann = ncfile.variables['tcuq'][:]
#mer_ann = ncfile.variables['tcvq'][:]
ncfile.close()

#ncfile = Dataset(sheddir+'wvfluxes_7914.nc','r')
#zon_ann = ncfile.variables['tcuq'][:]
#mer_ann = ncfile.variables['tcvq'][:]
#ncfile.close()

tcuq = pl.zeros([filenames3.shape[0],filenames3.shape[1],1,256,512])
tcvq = pl.zeros_like(tcuq)

u = pl.zeros([filenames1.shape[0],filenames1.shape[1],1,len(pres),256,512])
v = pl.zeros_like(u); q = pl.zeros_like(u)
surfp = pl.zeros([filenames2.shape[0],filenames2.shape[1],1,256,512])

#loop over years:q
for year in range(len(years)):
#    loop over filenames:
    for name in range(filenames1.shape[1]):
        #load ncfile
        nc1 = Dataset(panfs + year_input[year] + '/' + str(filenames1[year,name]),'r')
        #extract E & TP data
        u[year,name] = nc1.variables['U'][:]
        v[year,name] = nc1.variables['V'][:]
        q[year,name] = nc1.variables['Q'][:]
        nc1.close()
        nc2 = Dataset(panfs+year_input[year]+'/'+str(filenames2[year,name]),'r')
        surfp[year,name] = nc2.variables['SP'][:]
        nc2.close()
        nc3 = Dataset(panfs+year_input[year]+'/'+str(filenames3[year,name]),'r')
        tcuq[year,name] = nc3.variables['TCUQ'][:]
        tcvq[year,name] = nc3.variables['TCVQ'][:]
        nc3.close()

# remove the 1D axes:
#wv_zon = pl.squeeze(wv_zon); wv_mer = pl.squeeze(wv_mer)
u = pl.squeeze(u); v = pl.squeeze(v); q = pl.squeeze(q); surfp = pl.squeeze(surfp)
tcuq = pl.squeeze(tcuq); tcvq = pl.squeeze(tcvq)

tcuq_mn = pl.mean(tcuq,axis=0); tcvq_mn = pl.mean(tcvq,axis=0)
tcuq_DJF, tcuq_MAM, tcuq_JJA, tcuq_SON = SeasonalCyc(tcuq_mn)
tcvq_DJF, tcvq_MAM, tcvq_JJA, tcvq_SON = SeasonalCyc(tcvq_mn)

#U_mn = pl.mean(U,axis=(0,1)); V_mn = pl.mean(V,axis=(0,1)); Q_mn = pl.mean(Q,axis=(0,1))
#UQ_mn = pl.mean(U*Q)

#qv = q*v; qv = q*v
QU = pl.zeros([filenames1.shape[0],filenames1.shape[1],lat.size,lon.size])
QV = pl.zeros_like(QU)#; Q = pl.zeros_like(QU)

u_mn = pl.mean(u,axis=0); v_mn = pl.mean(v,axis=0)
q_mn = pl.mean(q,axis=0); sp_mn = pl.mean(surfp,axis=0)

uq_mn = pl.mean(u_mn*q_mn,axis=0); vq_mn = pl.mean(v_mn*q_mn,axis=0)

u_DJF, u_MAM, u_JJA, u_SON = SeasonalCyc(u_mn)
v_DJF, v_MAM, v_JJA, v_SON = SeasonalCyc(v_mn)
q_DJF, q_MAM, q_JJA, q_SON = SeasonalCyc(q_mn)
sp_DJF, sp_MAM, sp_JJA, sp_SON = SeasonalCyc(sp_mn)

QU_prod1 = pl.zeros([lat.size,lon.size]); QV_prod1 = pl.zeros_like(QU_prod)
QU_prod2 = pl.zeros_like(QU_prod); QV_prod2 = pl.zeros_like(QU_prod)
QU_prod3 = pl.zeros_like(QU_prod); QV_prod3 = pl.zeros_like(QU_prod)
QU_prod4 = pl.zeros_like(QU_prod); QV_prod4 = pl.zeros_like(QU_prod)

print 'Starting loops now!'
for i in range(lat.size):
    for j in range(lon.size):
        QU_prod1[i,j] = VertInt(u_DJF[:,i,j]*q_DJF[:,i,j],sp_DJF[i,j],pres,pres)
        QV_prod1[i,j] = VertInt(v_DJF[:,i,j]*q_DJF[:,i,j],sp_DJF[i,j],pres,pres)
        QU_prod2[i,j] = VertInt(u_MAM[:,i,j]*q_MAM[:,i,j],sp_MAM[i,j],pres,pres)
        QV_prod2[i,j] = VertInt(v_MAM[:,i,j]*q_MAM[:,i,j],sp_MAM[i,j],pres,pres)
        QU_prod3[i,j] = VertInt(u_JJA[:,i,j]*q_JJA[:,i,j],sp_JJA[i,j],pres,pres)
        QV_prod3[i,j] = VertInt(v_JJA[:,i,j]*q_JJA[:,i,j],sp_JJA[i,j],pres,pres)
        QU_prod4[i,j] = VertInt(u_SON[:,i,j]*q_SON[:,i,j],sp_SON[i,j],pres,pres)
        QV_prod4[i,j] = VertInt(v_SON[:,i,j]*q_SON[:,i,j],sp_SON[i,j],pres,pres)

#Q_mn = pl.mean(Q,axis=(0,1))
#QU_mn = pl.mean(zon_ann,axis=0); QV_mn = pl.mean(mer_ann,axis=0)
#QU_pert = QU_mn - QU_prod; QV_pert = QV_mn - QV_prod

QU_pert1 = tcuq_DJF - QU_prod1; QV_pert1 = tcvq_DJF - QV_prod1
QU_pert2 = tcuq_MAM - QU_prod2; QV_pert2 = tcvq_MAM - QV_prod2
QU_pert3 = tcuq_JJA - QU_prod3; QV_pert3 = tcvq_JJA - QV_prod3
QU_pert4 = tcuq_SON - QU_prod4; QV_pert4 = tcvq_SON - QV_prod4

#newnc = Dataset(sheddir+'reynolds/reyn_av_quv_SON.nc',mode='w',format='NETCDF4')
##
#lat_dim = newnc.createDimension('lat',lat.size)
#lon_dim = newnc.createDimension('lon',lon.size)
#lat_in = newnc.createVariable('lat',pl.float64,('lat',))
#lat_in.units = 'degrees_north'
#lat_in.long_name = 'latitude'
#lon_in = newnc.createVariable('lon',pl.float64,('lon',))
#lon_in.units = 'degrees_east'
#lon_in.long_name = 'longitude'
#lat_in = lat
#lon_in = lon
#
#QU_in = newnc.createVariable('Vertically-integrated zonal moisture flux (full field)',
#                                             pl.float64,('lat','lon'))
#QU_in.units = 'kg/m/s'
#QU_in.standard_name = 'zonal moisture flux'
#QU_in[:,:] = tcuq_SON
#
#QU_mf = newnc.createVariable('Vertically-integrated zonal moisture flux (mean flow)',
#                                             pl.float64,('lat','lon'))
#QU_mf.units = 'kg/m/s'
#QU_mf.standard_name = 'zonal moisture flux (mean flow)'
#QU_mf[:,:] = QU_prod4
#
#
#QU_ed = newnc.createVariable('Vertically-integrated zonal moisture flux (eddies)',
#                                             pl.float64,('lat','lon'))
#QU_ed.units = 'kg/m/s'
#QU_ed.standard_name = 'zonal moisture flux (eddies)'
#QU_ed[:,:] = QU_pert4
#
#QV_in = newnc.createVariable('Vertically-integrated meridional moisture flux (full field)',
#                                             pl.float64,('lat','lon'))
#QV_in.units = 'kg/m/s'
#QV_in.standard_name = 'meridional moisture flux'
#QV_in[:,:] = tcvq_SON
#
#QV_mf = newnc.createVariable('Vertically-integrated meridional moisture flux (mean flow)',
#                                             pl.float64,('lat','lon'))
#QV_mf.units = 'kg/m/s'
#QV_mf.standard_name = 'meridional moisture flux (mean flow)'
#QV_mf[:,:] = QV_prod4
#
#
#QV_ed = newnc.createVariable('Vertically-integrated meridional moisture flux (eddies)',
#                                             pl.float64,('lat','lon'))
#QV_ed.units = 'kg/m/s'
#QV_ed.standard_name = 'meridional moisture flux (eddies)'
#QV_ed[:,:] = QV_pert4
#newnc.close()