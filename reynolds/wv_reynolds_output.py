# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 11:20:05 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs

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
    
    return seglab

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

def ShedFluxes(sheddir,loc,lon,lat,zon,mer):
    """
    """
    endpts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_clicks.txt',skip_header=5)
    labs = TrajSegLabel(loc)
    l2 = pl.zeros([labs[0].size,3]); l2[:,0] = labs[0]; l2[:,1:] = labs[1]
    labs = l2.copy()
    nhat = NormalVector(sheddir+'shed_defs/',loc)
    labs = RR2(labs,loc)
    rlspts = pl.genfromtxt(sheddir+'shed_defs/'+loc+'_traj_release_new.txt',skip_header=5)
    rlspts = Add360(rlspts)
    midpt = MidPts(rlspts)
    midflux = MidFlux(rlspts,lon,lat,zon,mer)
    FdotN = NormalFlux(midflux,labs,nhat)
    fluxes = FdotN[:]#*midpt[:]/(10**9)
    
    return fluxes

pl.close('all')

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
exec(open('/home/np838619/Trajectory/trajfuncs.py').read())

sheddir = '/home/np838619/Watershed/'
reyndir = sheddir + 'reynolds/'

ncfile = Dataset(reyndir+'reyn_av_quv_7914_new.nc','r')
#lat = ncfile.variables['lat'][:]
#lon = ncfile.variables['lon'][:]
qu = ncfile.variables['Vertically-integrated zonal moisture flux (full field)'][:]
qu_mf = ncfile.variables['Vertically-integrated zonal moisture flux (mean flow)'][:]
qu_ed = ncfile.variables['Vertically-integrated zonal moisture flux (eddies)'][:]
qv = ncfile.variables['Vertically-integrated meridional moisture flux (full field)'][:]
qv_mf = ncfile.variables['Vertically-integrated meridional moisture flux (mean flow)'][:]
qv_ed = ncfile.variables['Vertically-integrated meridional moisture flux (eddies)'][:]
ncfile.close()

ncfile = Dataset(sheddir+'wvfluxes_7914.nc','r')
lon = ncfile.variables['lon'][:]
lat = ncfile.variables['lat'][:]
zon_ann = ncfile.variables['tcuq'][:]
mer_ann = ncfile.variables['tcvq'][:]
ncfile.close()
zon_ann = pl.mean(zon_ann,axis=0); mer_ann = pl.mean(mer_ann,axis=0)

lons,lats = pl.meshgrid(lon,lat)
proj = ccrs.Robinson(central_longitude=0)

###############################################################################
norm = pl.Normalize(-200,200)
levels = [-200,-150,-100,-50,-20,-10,10,20,50,100,150,200]
fig,ax = pl.subplots(3,2,figsize=(13,10.5))

ax1 = pl.subplot(321,projection=proj); ax1.coastlines(); ax1.gridlines()
cs1 = ax1.contourf(lons, lats, qu,transform=ccrs.PlateCarree(),levels=levels,
                  norm=norm,extend='both',cmap='seismic')
pl.title('(a) $\overline{qu}$',fontsize=20)

ax2 = pl.subplot(323,projection=proj); ax2.coastlines(); ax2.gridlines()
cs2 = ax2.contourf(lons,lats,qu_mf,transform=ccrs.PlateCarree(),levels=levels,
                   norm=norm,extend='both',cmap='seismic')
pl.title('(c) $\overline{q}\,\overline{u}$',fontsize=20)

ax3 = pl.subplot(325,projection=proj); ax3.coastlines(); ax3.gridlines()
cs3 = ax3.contourf(lons,lats,qu_ed,transform=ccrs.PlateCarree(),levels=levels,
                   norm=norm,extend='both',cmap='seismic')
pl.title('(e) $\overline{q\\backprime u\\backprime}$',fontsize=20)

#f = pl.gcf()
##colax = f.add_axes([0.29,0.03,0.45,0.05])
#colax = f.add_axes([0.83,0.22,0.05,0.6])
#clb = pl.colorbar(cs1, cax=colax,orientation='vertical')
#clb.set_ticks(levels)
#ticklabs = pl.asarray(levels)
#clb.set_ticklabels(ticklabs.astype('int'))
#clb.ax.tick_params(labelsize=14)
#clb.update_ticks()#; cb.ax.set_aspect(0.09)
#clb.set_label('kg/m/s',fontsize=20,labelpad=-2)
#
#pl.subplots_adjust(top=0.95,bottom=0.05,left=0.00)
#pl.savefig(reyndir+'reyn_qu_7914_new.png')
###############################################################################

###############################################################################
#fig,ax = pl.subplots(3,1,figsize=(6.5,9))

ax4 = pl.subplot(322,projection=proj); ax4.coastlines(); ax4.gridlines()
cs4 = ax4.contourf(lons, lats, qv,transform=ccrs.PlateCarree(),levels=levels,
                  norm=norm,extend='both',cmap='seismic')
pl.title('(b) $\overline{qv}$',fontsize=20)

ax5 = pl.subplot(324,projection=proj); ax5.coastlines(); ax5.gridlines()
cs5 = ax5.contourf(lons,lats,qv_mf,transform=ccrs.PlateCarree(),levels=levels,
                   norm=norm,extend='both',cmap='seismic')
pl.title('(d) $\overline{q}\,\overline{v}$',fontsize=20)

ax6 = pl.subplot(326,projection=proj); ax6.coastlines(); ax6.gridlines()
cs6 = ax6.contourf(lons,lats,qv_ed,transform=ccrs.PlateCarree(),levels=levels,
                   norm=norm,extend='both',cmap='seismic')
pl.title('(f) $\overline{q\\backprime v\\backprime}$',fontsize=20)

f = pl.gcf()
#colax = f.add_axes([0.29,0.03,0.45,0.05])
colax = f.add_axes([0.2,0.05,0.6,0.02])
clb = pl.colorbar(cs1, cax=colax,orientation='horizontal')
clb.set_ticks(levels)
ticklabs = pl.asarray(levels)
clb.set_ticklabels(ticklabs.astype('int'))
clb.ax.tick_params(labelsize=14)
clb.update_ticks()#; cb.ax.set_aspect(0.09)
clb.set_label('kg/m/s',fontsize=20,labelpad=-2)

pl.subplots_adjust(top=0.98,bottom=0.07,left=0.09,right=0.93,hspace=0.17,wspace=0.05)
pl.savefig(reyndir+'reyn_qu_qv_7914_new.png')
##############################################################################

###############################################################################
qu_zm = pl.mean(qu,axis=1); qu_za = qu - qu_zm[:,None]
qv_zm = pl.mean(qv,axis=1); qv_za = qv - qv_zm[:,None]


fig,ax = pl.subplots(2,1)

ax1 = pl.subplot(211,projection=proj); ax1.coastlines(); ax1.gridlines()
cs1 = ax1.contourf(lons,lats,qu_za,transform=ccrs.PlateCarree(),levels=levels,
                   norm=norm,extend='both',cmap='seismic')
pl.title('(a) $\overline{qu}^{*}$',fontsize=20)

ax2 = pl.subplot(212,projection=proj); ax2.coastlines(); ax2.gridlines()
cs2 = ax2.contourf(lons,lats,qv_za,transform=ccrs.PlateCarree(),levels=levels,
                   norm=norm,extend='both',cmap='seismic')
pl.title('(b) $\overline{qv}^{*}$',fontsize=20)

f = pl.gcf()
#colax = f.add_axes([0.29,0.03,0.45,0.05])
colax = f.add_axes([0.83,0.22,0.05,0.6])
clb = pl.colorbar(cs1, cax=colax,orientation='vertical')
clb.set_ticks(levels)
ticklabs = pl.asarray(levels)
clb.set_ticklabels(ticklabs.astype('int'))
clb.ax.tick_params(labelsize=14)
clb.update_ticks()#; cb.ax.set_aspect(0.09)
clb.set_label('kg/m/s',fontsize=20,labelpad=-2)

#pl.savefig(reyndir+'zonalmean_devs7914.png')
###############################################################################

amr = ShedFluxes(sheddir,'Am',lon,lat,qu,qv)
amr2 = ShedFluxes(sheddir,'Am',lon,lat,qu_mf,qv_mf)
amr3 = ShedFluxes(sheddir,'Am',lon,lat,qu_ed,qv_ed)

afr = ShedFluxes(sheddir,'AfMe',lon,lat,qu,qv)
afr2 = ShedFluxes(sheddir,'AfMe',lon,lat,qu_mf,qv_mf)
afr3 = ShedFluxes(sheddir,'AfMe',lon,lat,qu_ed,qv_ed)

eaa = ShedFluxes(sheddir,'EAA',lon,lat,qu,qv)
eaa2 = ShedFluxes(sheddir,'EAA',lon,lat,qu_mf,qv_mf)
eaa3 = ShedFluxes(sheddir,'EAA',lon,lat,qu_ed,qv_ed)

ara = ShedFluxes(sheddir,'ArA',lon,lat,qu,qv)
ara2 = ShedFluxes(sheddir,'ArA',lon,lat,qu_mf,qv_mf)
ara3 = ShedFluxes(sheddir,'ArA',lon,lat,qu_ed,qv_ed)

ari = ShedFluxes(sheddir,'ArI',lon,lat,qu,qv)
ari2 = ShedFluxes(sheddir,'ArI',lon,lat,qu_mf,qv_mf)
ari3 = ShedFluxes(sheddir,'ArI',lon,lat,qu_ed,qv_ed)

arp = ShedFluxes(sheddir,'ArP',lon,lat,qu,qv)
arp2 = ShedFluxes(sheddir,'ArP',lon,lat,qu_mf,qv_mf)
arp3 = ShedFluxes(sheddir,'ArP',lon,lat,qu_ed,qv_ed)

soa = ShedFluxes(sheddir,'SOA',lon,lat,qu,qv)
soa2 = ShedFluxes(sheddir,'SOA',lon,lat,qu_mf,qv_mf)
soa3 = ShedFluxes(sheddir,'SOA',lon,lat,qu_ed,qv_ed)

soi = ShedFluxes(sheddir,'SOI',lon,lat,qu,qv)
soi2 = ShedFluxes(sheddir,'SOI',lon,lat,qu_mf,qv_mf)
soi3 = ShedFluxes(sheddir,'SOI',lon,lat,qu_ed,qv_ed)

sop = ShedFluxes(sheddir,'SOP',lon,lat,qu,qv)
sop2 = ShedFluxes(sheddir,'SOP',lon,lat,qu_mf,qv_mf)
sop3 = ShedFluxes(sheddir,'SOP',lon,lat,qu_ed,qv_ed)

F = pl.array([[pl.sum(amr),pl.sum(amr2),pl.sum(amr3)],
              [pl.sum(afr),pl.sum(afr2),pl.sum(afr3)],
                [pl.sum(eaa),pl.sum(eaa2),pl.sum(eaa3)],
                [pl.sum(ara),pl.sum(ara2),pl.sum(ara3)],
                [pl.sum(ari),pl.sum(ari2),pl.sum(ari3)],
                [pl.sum(arp),pl.sum(arp2),pl.sum(arp3)],
                [pl.sum(soa),pl.sum(soa2),pl.sum(soa3)],
                [pl.sum(soi),pl.sum(soi2),pl.sum(soi3)],
                [pl.sum(sop),pl.sum(sop2),pl.sum(sop3)]])

#f = open(reyndir+'reyn_fluxes_new.csv','w')
#pl.savetxt(f,F,fmt='%9.5f',delimiter=',')
#f.close()

amr_pts = pl.genfromtxt(sheddir+'shed_defs/Am_traj_release_new.txt',skip_header=5)
afr_pts = pl.genfromtxt(sheddir+'shed_defs/AfMe_traj_release_new.txt',skip_header=5)
eaa_pts = pl.genfromtxt(sheddir+'shed_defs/EAA_traj_release_new.txt',skip_header=5)
ara_pts = pl.genfromtxt(sheddir+'shed_defs/ArA_traj_release_new.txt',skip_header=5)
ari_pts = pl.genfromtxt(sheddir+'shed_defs/ArI_traj_release_new.txt',skip_header=5)
arp_pts = pl.genfromtxt(sheddir+'shed_defs/ArP_traj_release_new.txt',skip_header=5)
soa_pts = pl.genfromtxt(sheddir+'shed_defs/SOA_traj_release_new.txt',skip_header=5)
soi_pts = pl.genfromtxt(sheddir+'shed_defs/SOI_traj_release_new.txt',skip_header=5)
sop_pts = pl.genfromtxt(sheddir+'shed_defs/SOP_traj_release_new.txt',skip_header=5)

labs = ['$\overline{\mathbf{Q}}\cdot\mathbf{\hat{n}}$',
        '$\overline{q}\;\overline{\mathbf{u}}\cdot\mathbf{\hat{n}}$',
        '$\overline{q\\backprime \mathbf{u}\\backprime}\cdot\mathbf{\hat{n}}$']
fig, ax = pl.subplots(3,3,figsize=(12.5,12))

ax1 = pl.subplot(331)
ax1.plot(amr,color='b',lw=2,zorder=3,label=labs[0])
ax1.plot(amr2,color='darkmagenta',lw=2,zorder=2,label=labs[1])
ax1.plot(amr3,color='lightblue',lw=2,zorder=1,label=labs[2])
pl.ylim(-300,200); pl.xlim(0,155); ax1.grid(axis='y')
pl.yticks(pl.arange(-300,210,50),fontsize=13); pl.ylabel('kg/m/s',fontsize=16)
pl.text(3,-290,'(a) Americas',fontsize=20)
ax1.legend(loc=1,ncol=2,columnspacing=0.5,handletextpad=0.5)
a = pl.arange(0,155,20); b = pl.around(amr_pts[a,1],decimals=0)
ax1.set_xticks(a); ax1.set_xticklabels(b.astype('int'),fontsize=14)
ax1.tick_params(axis='y',pad=5); ax1.tick_params(axis='x',pad=7)

ax2 = pl.subplot(332)
ax2.plot(afr,color='k',lw=2,zorder=3,label=labs[0])
ax2.plot(afr2,color='dimgray',lw=2,zorder=2,label=labs[1])
ax2.plot(afr3,color='silver',lw=2.,zorder=1,label=labs[2])
pl.ylim(-300,200); pl.xlim(0,214); ax2.grid(axis='y')
pl.yticks(pl.arange(-300,210,50))
ax2.yaxis.set_major_formatter(pl.NullFormatter())
pl.text(3,-290,'(b) Africa',fontsize=20)
ax2.legend(loc=4,ncol=1,columnspacing=0.5,handletextpad=0.5)
a = pl.arange(0,214,30); b = pl.around(afr_pts[a,1],decimals=0)
ax2.set_xticks(a); ax2.set_xticklabels(b.astype('int'),fontsize=14)
pl.xlabel('latitude',fontsize=16); ax2.tick_params(axis='x',pad=7)

ax3 = pl.subplot(333)
ax3.plot(eaa,color='g',lw=2,zorder=3,label=labs[0])
ax3.plot(eaa2,color='olive',lw=2,zorder=2,label=labs[1])
ax3.plot(eaa3,color='lime',lw=2.,zorder=1,label=labs[2])
pl.ylim(-300,200); pl.xlim(0,165); ax3.grid(axis='y')
pl.yticks(pl.arange(-300,210,50))
ax3.yaxis.set_major_formatter(pl.NullFormatter())
pl.text(3,-290,'(c) South-East Asia',fontsize=20)
ax3.legend(loc=(0.03,0.14),ncol=2,columnspacing=0.5,handletextpad=0.5)
a = pl.arange(0,165,20); b = pl.around(eaa_pts[a,1],decimals=0)
ax3.set_xticks(a); ax3.set_xticklabels(b.astype('int'),fontsize=14)
ax3.tick_params(axis='x',pad=7)

ax4 = pl.subplot(334)
ax4.plot(ara,color='r',lw=2,zorder=3,label=labs[0])
ax4.plot(ara2,color='maroon',lw=2,zorder=2,label=labs[1])
ax4.plot(ara3,color='lightsalmon',lw=2.,zorder=1)
pl.ylim(-100,80); pl.xlim(0,197); ax4.grid(axis='y'); pl.yticks(fontsize=13)
pl.ylabel('kg/m/s',fontsize=16)
pl.text(50,-95,'(d) Arctic Atlantic',fontsize=20)
ax4.legend(loc=2,ncol=1,columnspacing=0.5,handletextpad=0.5)
a = pl.arange(0,197,30); b = pl.around(ara_pts[a,0],decimals=0)
ax4.set_xticks(a); ax4.set_xticklabels(b.astype('int'),fontsize=14)
ax4.tick_params(axis='y',pad=5); ax4.tick_params(axis='x',pad=7)

ax5 = pl.subplot(335)
ax5.plot(ari,color='r',lw=2,zorder=3)
ax5.plot(ari2,color='maroon',lw=2,zorder=2)
ax5.plot(ari3,color='lightsalmon',lw=2.,zorder=1,label=labs[2])
pl.ylim(-100,80); pl.xlim(0,102); ax5.grid(axis='y')
ax5.yaxis.set_major_formatter(pl.NullFormatter())
pl.text(30,-95,'(e) Arctic Indian',fontsize=20)
ax5.legend(loc=1,ncol=2,columnspacing=0.5,handletextpad=0.5)
a = pl.arange(0,102,20); b = pl.around(ari_pts[a,0],decimals=0)
ax5.set_xticks(a); ax5.set_xticklabels(b.astype('int'),fontsize=14)
pl.xlabel('longitude',fontsize=16); ax5.tick_params(axis='x',pad=7)

ax6 = pl.subplot(336)
ax6.plot(arp,color='r',lw=2,zorder=3)
ax6.plot(arp2,color='maroon',lw=2,zorder=2)
ax6.plot(arp3,color='lightsalmon',lw=2.,zorder=1)
pl.ylim(-100,80); pl.xlim(0,220); ax6.grid(axis='y')
ax6.yaxis.set_major_formatter(pl.NullFormatter())
pl.text(70,-95,'(f) Arctic Pacific',fontsize=20)
a = pl.arange(0,220,30); b = pl.around(arp_pts[a,0],decimals=0)
ax6.set_xticks(a); ax6.set_xticklabels(b.astype('int'),fontsize=14)
ax6.tick_params(axis='x',pad=7)

ax7 = pl.subplot(337)
ax7.plot(soa,color='darkgoldenrod',lw=2,zorder=3,label=labs[0])
ax7.plot(soa2,color='darkorange',lw=2,zorder=2)
ax7.plot(soa3,color='gold',lw=2.,zorder=1)
pl.ylim(-200,100); pl.xlim(0,185); ax7.grid(axis='y'); pl.yticks(fontsize=13)
pl.ylabel('kg/m/s',fontsize=16)
pl.text(25,-190,'(g) Southern Atlantic',fontsize=20)
ax7.legend(loc=1,ncol=1,columnspacing=0.5,handletextpad=0.5)
a = pl.arange(0,185,30); b = pl.around(soa_pts[a,0],decimals=0)
ax7.set_xticks(a); ax7.set_xticklabels(b.astype('int'),fontsize=14)
ax7.tick_params(axis='y',pad=5); ax7.tick_params(axis='x',pad=7)

ax8 = pl.subplot(338)
ax8.plot(soi,color='darkgoldenrod',lw=2,zorder=3)
ax8.plot(soi2,color='darkorange',lw=2,zorder=2,label=labs[1])
ax8.plot(soi3,color='gold',lw=2.,zorder=1,label=labs[2])
pl.ylim(-200,100); pl.xlim(0,184); ax8.grid(axis='y')
ax8.yaxis.set_major_formatter(pl.NullFormatter())
pl.text(30,-190,'(h) Southern Indian',fontsize=20)
ax8.legend(loc=2,ncol=2,columnspacing=0.5,handletextpad=0.5)
a = pl.arange(0,184,30); b = pl.around(soi_pts[a,0],decimals=0)
ax8.set_xticks(a); ax8.set_xticklabels(b.astype('int'),fontsize=14)
pl.xlabel('longitude',fontsize=16); ax8.tick_params(axis='x',pad=7)

ax9 = pl.subplot(339)
ax9.plot(sop,color='darkgoldenrod',lw=2,zorder=3)
ax9.plot(sop2,color='darkorange',lw=2,zorder=2)
ax9.plot(sop3,color='gold',lw=2.,zorder=1)
pl.ylim(-200,100); pl.xlim(0,212); ax9.grid(axis='y')
ax9.yaxis.set_major_formatter(pl.NullFormatter())
pl.text(40,-190,'(i) Southern Pacific',fontsize=20)
sop_pts[63:,0] = sop_pts[63:,0] - 360
a = pl.arange(0,212,30); b = pl.around(sop_pts[a,0],decimals=0)
ax9.set_xticks(a); ax9.set_xticklabels(b.astype('int'),fontsize=14)
ax9.tick_params(axis='x',pad=7)

pl.tight_layout()
#pl.savefig(reyndir+'reyn_profs_7914_new.png')