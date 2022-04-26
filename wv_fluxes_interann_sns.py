# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:15:47 2017

@author: np838619
"""

import pylab as pl
from netCDF4 import Dataset

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
    fluxes = FdotN[:]*midpt[:]/(10**9)
    
    return fluxes

def Seasons(flux):
    """
    """
    flx_sns = pl.zeros([flux.shape[0],4])
    flx_sns[:,1] = pl.mean(flux[:,2:5],axis=1)
    flx_sns[:,2] = pl.mean(flux[:,5:8],axis=1)
    flx_sns[:,3] = pl.mean(flux[:,8:11],axis=1)
    flx_sns[1:,0] = (flux[:-1,-1]+flux[1:,1]+flux[1:,2])/3
    flx_sns[0,0] = pl.float32('nan')
    
    return flx_sns

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
zon_ann = pl.mean(wv_zon,axis=1); mer_ann = pl.mean(wv_mer,axis=1)

sheds = ['Am','AfMe','EAA','ArA','ArI','ArP','SOA','SOI','SOP']
F1 = pl.zeros([36,12]); F2 = pl.zeros([36,12]); F3 = pl.zeros([36,12])
F4 = pl.zeros([36,12]); F5 = pl.zeros([36,12]); F6 = pl.zeros([36,12])
F7 = pl.zeros([36,12]); F8 = pl.zeros([36,12]); F9 = pl.zeros([36,12])
for yr in range(filenames.shape[0]):
    for mt in range(12):
        F1[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[0],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # amr
        F2[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[1],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # afr
        F3[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[2],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # eaa
        F4[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[3],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # ara
        F5[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[4],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # ari
        F6[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[5],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # arp
        F7[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[6],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # soa
        F8[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[7],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # soi
        F9[yr,mt]=pl.sum(ShedFluxes(sheddir,sheds[8],lon,lat,wv_zon[yr,mt],wv_mer[yr,mt])) # sop

amr_sns = Seasons(F1); afr_sns = Seasons(F2); eaa_sns = Seasons(F3)
ara_sns = Seasons(F4); ari_sns = Seasons(F5); arp_sns = Seasons(F6)
soa_sns = Seasons(F7); soi_sns = Seasons(F8); sop_sns = Seasons(F9)

#F1 = pl.asarray(F1); F2 = pl.asarray(F2); F3 = pl.asarray(F3)
#F4 = pl.asarray(F4); F5 = pl.asarray(F5); F6 = pl.asarray(F6)
#F7 = pl.asarray(F7); F8 = pl.asarray(F8); F9 = pl.asarray(F9)
#
#F1s = pl.sum(F1,axis=1); F2s = pl.sum(F2,axis=1); F3s = pl.sum(F3,axis=1)
#F4s = pl.sum(F4,axis=1); F5s = pl.sum(F5,axis=1); F6s = pl.sum(F6,axis=1)
#F7s = pl.sum(F7,axis=1); F8s = pl.sum(F8,axis=1); F9s = pl.sum(F9,axis=1)
#
#F = pl.array([F1s,F2s,F3s,F4s,F5s,F6s,F7s,F8s,F9s])
#
#atl = F1s - F2s - F4s + F7s
#ind = F2s - F3s - F5s + F8s
#pac = F3s - F1s - F6s + F9s
#arc = F4s + F5s + F6s
#sou = -1*(F7s + F8s + F9s)
#D = pl.array([atl,ind,pac,arc,sou])

#f = open(sheddir+'variability/'+'shedfluxes_7914_var.csv','w')
#f.write('atl ind pac arc sou\n')
#pl.savetxt(f,F.T,fmt='%9.5f',delimiter=' ')
#f.close()

years = pl.linspace(1979,2014,36)
C = ['b','purple','r','k']
lw=2

fig,ax = pl.subplots(3,3,figsize=(14,10))

ax1 = pl.subplot(331)
for i in range(4):
    ax1.plot(years,amr_sns[:,i],color=C[i],linewidth=lw)
ax1.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(-0.6,0)
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylabel('Sv',fontsize=16); pl.title('(a) Americas',fontsize=16,loc='left')
pl.yticks(fontsize=14); ax1.grid(axis='both')

ax2 = pl.subplot(332)
for i in range(4):
    ax2.plot(years,afr_sns[:,i],color=C[i],linewidth=lw)
ax2.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(-0.4,0.2)
ax2.xaxis.set_major_formatter(pl.NullFormatter())
ax2.legend(['DJF','MAM','JJA','SON'],ncol=2,loc=2)
pl.title('(b) Africa',fontsize=16,loc='left')
pl.yticks(fontsize=14); ax2.grid(axis='both')

ax3 = pl.subplot(333)
for i in range(4):
    ax3.plot(years,eaa_sns[:,i],color=C[i],linewidth=lw)
ax3.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(-0.4,0.8)
ax3.xaxis.set_major_formatter(pl.NullFormatter())
pl.title('(c) South-East Asia',fontsize=16,loc='left')
pl.yticks(fontsize=14); ax3.grid(axis='both')

ax4 = pl.subplot(334)
for i in range(4):
    ax4.plot(years,ara_sns[:,i],color=C[i],linewidth=lw)
ax4.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(0,0.3)
ax4.xaxis.set_major_formatter(pl.NullFormatter())
pl.ylabel('Sv',fontsize=16)
pl.title('(d) Arctic Atlantic',fontsize=16,loc='left')
pl.yticks(fontsize=14); ax4.grid(axis='both')

ax5 = pl.subplot(335)
for i in range(4):
    ax5.plot(years,ari_sns[:,i],color=C[i],linewidth=lw)
ax5.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(-0.1,0.2)
ax5.xaxis.set_major_formatter(pl.NullFormatter())
pl.title('(e) Arctic Indian',fontsize=16,loc='left')
pl.yticks(fontsize=14); ax5.grid(axis='both')

ax6 = pl.subplot(336)
for i in range(4):
    ax6.plot(years,arp_sns[:,i],color=C[i],linewidth=lw)
ax6.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(-0.1,0.2)
ax6.xaxis.set_major_formatter(pl.NullFormatter())
pl.title('(f) Arctic Pacific',fontsize=16,loc='left')
pl.yticks(fontsize=14); ax6.grid(axis='both')

ax7 = pl.subplot(337)
for i in range(4):
    ax7.plot(years,soa_sns[:,i],color=C[i],linewidth=lw)
ax7.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(-0.4,-0.1)
pl.xticks(rotation=45,fontsize=16)
pl.ylabel('Sv',fontsize=16)
pl.title('(g) Southern Atlantic',fontsize=14,loc='left')
pl.yticks(fontsize=14); ax7.grid(axis='both')

ax8 = pl.subplot(338)
for i in range(4):
    ax8.plot(years,soi_sns[:,i],color=C[i],linewidth=lw)
ax8.set_xticks(years[::5]); pl.xlim(1979,2014)
pl.xticks(rotation=45,fontsize=16)
pl.title('(h) Southern Indian',fontsize=14,loc='left')
pl.yticks(fontsize=14); ax8.grid(axis='both')

ax9 = pl.subplot(339)
for i in range(4):
    ax9.plot(years,sop_sns[:,i],color=C[i],linewidth=lw)
ax9.set_xticks(years[::5]); pl.xlim(1979,2014); pl.ylim(-0.5,-0.2)
pl.xticks(rotation=45,fontsize=16)
pl.title('(i) Southern Pacific',fontsize=14,loc='left')
pl.yticks(fontsize=14); ax9.grid(axis='both')

pl.tight_layout()
#pl.savefig(sheddir+'shedfluxes_ssns_iv.png')