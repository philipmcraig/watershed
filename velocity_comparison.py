# -*- coding: utf-8 -*-
"""
Created on Wed May 18 15:30:58 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap, shiftgrid

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


exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
pl.close('all')
# need to read in one particular file: 21 July 2007 1200
trajdir = '/home/np838619/Trajectory/'
eradir = '/panfs/jasmin/era/era-in/netc/'; sheddir = '/home/np838619/Watershed/'
file_plev = 'ggap'; file_surf = 'ggas'; file_inv = 'ggis'
year = '2007'; mnth_name = 'jul'
#filedir = eradir+filetype+'/'+year+'/'+mnth_name+year+'/'
 

mnth_no = '07'; day = '21'; time = '1200'
fname1 = file_plev + year + mnth_no + day + time
fname2 = file_surf + year + mnth_no + day + time
fname3 = file_inv + '198901010000.nc'

nc1 = Dataset(sheddir+fname1+'.nc','r')
lon = nc1.variables['longitude'][:]
lat = nc1.variables['latitude'][:]
pres = nc1.variables['p'][:]
u = nc1.variables['U'][:]; v = nc1.variables['V'][:]
q = nc1.variables['Q'][:]
nc1.close()

nc2 = Dataset(sheddir+fname2+'.nc','r')
surfp = nc2.variables['SP'][:]
nc2.close()

nc3 = Dataset(sheddir+fname3,'r')
geop = nc3.variables['Z'][:]
nc3.close()

nc4 = Dataset('/home/np838619/Downloads/etopo05.nc','r')
topo = nc4.variables['ROSE'][:]
toplon = nc4.variables['ETOPO05_X'][:]
toplat = nc4.variables['ETOPO05_Y'][:]
nc4.close()
toplat = pl.flipud(toplat); topo = pl.flipud(topo)

# remove 1D axes:
pres = pl.squeeze(pres); u = pl.squeeze(u); v = pl.squeeze(v); q = pl.squeeze(q)
surfp = pl.squeeze(surfp); geop = pl.squeeze(geop)

# calculate orography field by dividing geopotential by gravity:
orog = geop/9.80665

trajlevs = pl.linspace(950,200,16)# pressure levels used in traj model

#read in velocity data from traj model:

trajvels = pl.genfromtxt(trajdir + 'uv_3day07.txt')#; trajvels = -trajvels
trajq = pl.genfromtxt(trajdir + 'q_3day07.txt')
trajp = pl.genfromtxt(trajdir + 'p_3day07.txt')

# read in traj release points:
relpts = pl.genfromtxt(sheddir + 'NCA_traj_release.txt',skip_header=5)
rl2 = pl.genfromtxt(sheddir + 'NCA_traj_release_new.txt',skip_header=5)
relpts[:,0] = relpts[:,0] + 360.; rl2[:,0] = rl2[:,0] + 360.
trajvels = pl.reshape(trajvels,(len(trajlevs),len(relpts),2))
trajq = pl.reshape(trajq,(len(trajlevs),len(relpts)))
trajp = pl.reshape(trajp,(len(trajlevs),len(relpts)))

# set up array for interpolated velocities:
vel_int = pl.zeros([len(trajlevs),len(relpts),2])
# interpolated specific humidities:
q_int = pl.zeros([len(trajlevs),len(relpts)])
surfp_int = pl.zeros([len(relpts)]); orog_int = pl.zeros([len(rl2)])
topo_int = pl.zeros_like(orog_int)
for lev in range(len(trajlevs)):
    L = pl.where(pres==trajlevs[lev])#; print lev
    for point in range(len(relpts)):
        vel_int[lev,point,0] = BilinInterp(relpts[point],lon,lat,u[L[0][0]])
        vel_int[lev,point,1] = BilinInterp(relpts[point],lon,lat,v[L[0][0]])
        q_int[lev,point] = BilinInterp(relpts[point],lon,lat,q[L[0][0]])

for point in range(len(relpts)):
    surfp_int[point] = BilinInterp(relpts[point],lon,lat,surfp)
for point in range(len(rl2)):
    orog_int[point] = BilinInterp(rl2[point],lon,lat,orog)
    topo_int[point] = BilinInterp(rl2[point],toplon,toplat,topo)
#pl.figure(1)
#for i in range(182):
#    pl.plot(trajq[:,i],trajp[:,i],marker='x',ls='None')
#pl.plot(pl.mean(trajq,axis=1),trajp[:,0],color='k',linewidth=2)
#pl.plot(trajq.max(axis=1),trajp[:,0],marker='+',ls='None')
#pl.ylim(1100,0)

#for point in range(len(relpts)):
#    for lev in range(len(trajlevs)):
#        if trajlevs[lev] > surfp_int[point]/100:
#            #vel_int[lev,point,:] = pl.float64('nan')
#            #trajvels[lev,point,:] = pl.float64('nan')
#            q_int[lev,point] = pl.float64('nan')
#            #trajq[lev,point] = pl.float64('nan')

# Calculate errors:
u_diff = trajvels[:,:,0] - vel_int[:,:,0]; mean_err_u = pl.nanmean(u_diff,axis=1)
v_diff = trajvels[:,:,1] - vel_int[:,:,1]; mean_err_v = pl.nanmean(v_diff,axis=1)
q_diff = trajq - q_int; mean_err_q = pl.mean(q_diff,axis=1)

# Standard deviations:
u_std = pl.std(u_diff,axis=1); v_std = pl.std(v_diff,axis=1)
q_std = pl.std(q_diff,axis=1)

fig, ax = pl.subplots(1,3,figsize=(9,6)); lw = 2.
ax1 = pl.subplot(1,3,1)
pl.plot(pl.mean(trajvels[:,:,0],axis=1),pl.mean(trajp,axis=1),label='$u_{traj}$',color='b',lw=lw)
pl.plot(pl.mean(vel_int[:,:,0],axis=1),trajlevs,label='$u_{ERA}$',ls='--',color='b',lw=lw)
pl.plot(pl.mean(trajvels[:,:,1],axis=1),pl.mean(trajp,axis=1),label='$v_{traj}$',color='g',lw=lw)
pl.plot(pl.mean(vel_int[:,:,1],axis=1),trajlevs,label='$v_{ERA}$',ls='--',color='g',lw=lw)
pl.ylabel('Pressure (hPa)',fontsize=20); pl.ylim(1000,200)
pl.xlabel('m/s',fontsize=20); pl.legend(loc=0)
pl.title('Mean vertical profiles'); pl.axvline(x=0,color='k',ls='--')

ax2 = pl.subplot(1,3,2)
pl.plot(mean_err_u,trajlevs,label='$u$ error',lw=lw)
#pl.errorbar(mean_err_u,trajlevs,xerr=u_std,label='$u$ error',capsize=5); pl.ylim(1000,200)
pl.plot(mean_err_v,trajlevs,label='$v$ error',lw=lw); pl.ylim(1000,200)
#pl.errorbar(mean_err_v,trajlevs,xerr=v_std,label='$v$ error',capsize=5); pl.ylim(1000,200)
pl.legend(loc=2); pl.axvline(x=0,ls='--',color='k')
pl.title('Mean difference'); pl.xlabel('m/s',fontsize=20)
ax2.yaxis.set_major_formatter(pl.NullFormatter())
#pl.savefig(trajdir+'velocity_diffs.png')

ax3 = pl.subplot(1,3,3)#pl.figure(2)
pl.plot(u_std,trajlevs,label='$u_\sigma$',lw=lw)
pl.plot(v_std,trajlevs,label='$v_\sigma$',lw=lw)
pl.ylim(1000,200); pl.legend(loc=0)
pl.title('Standard deviation');pl.xlabel('m/s',fontsize=20)
ax3.yaxis.set_major_formatter(pl.NullFormatter())
#pl.savefig(trajdir+'velocity_std.png')

#pl.close(); pl.close()

fig,ax = pl.subplots(1,2,figsize=(9,6))#pl.figure(4)
pl.subplot(1,2,1)
pl.plot(pl.mean(q_int,axis=1)*1000,trajlevs,label='$q$ analysis')
pl.plot(pl.mean(trajq,axis=1)*1000,trajlevs,label='$q$ trajectory')
pl.ylim(1000,200); pl.legend(); pl.title('Mean vertical profiles',fontsize=16)
pl.xlabel('g/kg',fontsize=20); pl.ylabel('hPa',fontsize=20)
pl.subplot(1,2,2)
pl.plot(mean_err_q*1000,trajlevs,color='r',linewidth=2,label='$\Delta q$')
pl.plot(q_std*1000,trajlevs,color='purple',linewidth=2,label='$\Delta q_\sigma$')
pl.ylim(1000,200); pl.title('Trajectory minus analysis',fontsize=16)
pl.axvline(x=0,ls='--',color='k'); pl.xlabel('g/kg',fontsize=20); pl.legend(loc=0)


X = pl.linspace(70,10,7)
Y = pl.zeros_like(X)
for i in range(len(X)):
    Y[i] = NearestIndex(relpts[:,1],X[i])

tj_bnds_mns = pl.zeros([len(X),len(trajlevs),2])
vi_bnds_mns = pl.zeros_like(tj_bnds_mns)
uv_sd = pl.zeros_like(tj_bnds_mns)
tj_bnds_mns[0,:,0] = pl.mean(trajvels[:,:Y[0]+1,0],axis=1)
tj_bnds_mns[0,:,1] = pl.mean(trajvels[:,:Y[0]+1,1],axis=1)
vi_bnds_mns[0,:,0] = pl.mean(vel_int[:,:Y[0]+1,0],axis=1)
vi_bnds_mns[0,:,1] = pl.mean(vel_int[:,:Y[0]+1,1],axis=1)
uv_sd[0,:,0] = pl.std(u_diff[:,:Y[0]+1],axis=1)
uv_sd[0,:,1] = pl.std(v_diff[:,:Y[0]+1],axis=1)
for bnd in range(0,len(X)-1):
    tj_bnds_mns[bnd+1,:,0] = pl.mean(trajvels[:,Y[bnd]:Y[bnd+1]+1,0],axis=1)
    tj_bnds_mns[bnd+1,:,1] = pl.mean(trajvels[:,Y[bnd]:Y[bnd+1]+1,1],axis=1)
    vi_bnds_mns[bnd+1,:,0] = pl.mean(vel_int[:,Y[bnd]:Y[bnd+1]+1,0],axis=1)
    vi_bnds_mns[bnd+1,:,1] = pl.mean(vel_int[:,Y[bnd]:Y[bnd+1]+1,1],axis=1)
    uv_sd[bnd+1,:,0] = pl.std(u_diff[:,Y[bnd]:Y[bnd+1]+1],axis=1)
    uv_sd[bnd+1,:,1] = pl.std(v_diff[:,Y[bnd]:Y[bnd+1]+1],axis=1)
#tj_bnds_mns[-1,:,0] = pl.mean(trajvels[:,Y[-2]:Y[-1]+1,0],axis=1)
#tj_bnds_mns[-1,:,1] = pl.mean(trajvels[:,Y[-2]:Y[-1]+1,1],axis=1)

diffs = tj_bnds_mns - vi_bnds_mns

ax, fig = pl.subplots(2,4,figsize=(20,14))
bnds = ['Alaska','70N-60N','60N-50N','50N-40N','40N-30N','30N-20N','20N-10N','']
for i in range(len(X)):
    pl.subplot(2,4,i+1)
    #pl.plot(tj_bnds_mns[i,:,0],trajlevs,label='traj')
    #pl.plot(vi_bnds_mns[i,:,0],trajlevs,label='ERA')
    pl.plot(diffs[i,:,0],trajlevs,label='$u_{diff}$',color='b')
    pl.plot(diffs[i,:,1],trajlevs,label='$v_{diff}$',color='g')
    pl.plot(uv_sd[i,:,0],trajlevs,label='$u_\sigma$',color='b',ls='--')
    pl.plot(uv_sd[i,:,1],trajlevs,label='$v_\sigma$',color='g',ls='--')
    pl.ylim(1000,200)#; pl.legend(loc=0)
    pl.xlim(-5,5); pl.title(bnds[i],fontsize=18)
    pl.axvline(x=0,color='k',ls='--')
    if i in (0,4):
        pl.ylabel('Pressure (hPa)',fontsize=18)
    if i in (4,5,6,7):
        pl.xlabel('$u$ (m/s)',fontsize=18)
pl.legend(loc=(1.5,0.5),fontsize=20)
pl.suptitle('Trajectory minus analysis',fontsize=20)

ax, fig = pl.subplots(3,2,figsize=(24,15))
m = Basemap(projection='cyl',resolution='l',llcrnrlat=0.,urcrnrlat=80.,
                llcrnrlon=-180.,urcrnrlon=-50.,lat_ts=20)

u, lons = shiftgrid(180.0, u, lon, start=False) # shifting the grid works for some reason
lons, lats = pl.meshgrid(lons,lat)
X, Y = m(lons,lats)
levs = pl.arange(-30,35,5)#pl.linspace(-30,30,11)
levinds = [11,15,19]; levnames = [str(int(pres[11])),str(int(pres[15])),str(int(pres[19]))]
for i in range(len(levinds)):
    pl.subplot(3,2,2*i+1)
    cs=pl.contourf(X,Y,u[levinds[i]],cmap='seismic',norm=pl.Normalize(-30,30),extend='both',levels=levs)
    pl.plot(relpts[:,0]-360,relpts[:,1],color='k',linewidth=1.5)
    m.drawcoastlines(); pl.title('$u$'+levnames[i])
    m.drawparallels([10,20,30,40,50,60,70],labels=[1,1,1,1],linewidth=0.5,ax=None)
    
    pl.subplot(3,2,2*i+2)
    cs=pl.contourf(X,Y,v[levinds[i]],cmap='seismic',norm=pl.Normalize(-30,30),extend='both',levels=levs)
    pl.plot(relpts[:,0]-360,relpts[:,1],color='k',linewidth=1.5)
    m.drawcoastlines(); pl.title('$v$'+levnames[i])
    m.drawparallels([10,20,30,40,50,60,70],labels=[1,1,1,1],linewidth=0.5,ax=None)
#pl.subplot(3,1,2)
#pl.contourf(X,Y,u[15],cmap='seismic',norm=pl.Normalize(-30,30),extend='both',levels=levs)
#m.drawcoastlines(); pl.title('500 hPa')
#m.drawparallels([10,20,30,40,50,60,70],labels=[1,1,1,1],linewidth=0.5,ax=None)
#pl.subplot(3,1,3)
#pl.contourf(X,Y,u[19],cmap='seismic',norm=pl.Normalize(-30,30),extend='both',levels=levs)
#m.drawcoastlines(); pl.title('300 hPa')
#m.drawparallels([10,20,30,40,50,60,70],labels=[1,1,1,1],linewidth=0.5,ax=None)

f=pl.gcf();
colax = f.add_axes([0.41,0.1,0.21,0.03])
#pcmbounds = pl.linspace(-30,30,11)
clb = pl.colorbar(cs,extend='max',cax=colax,orientation='horizontal',ticks=[-30,-20,-10,0,10,20,30])
clb.set_label('m/s',fontsize=20)
clb.ax.tick_params(labelsize=15)
pl.subplots_adjust(top=0.95,bottom=0.15,wspace=-0.5)