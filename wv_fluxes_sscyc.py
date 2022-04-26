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
    fluxes = FdotN[:]#*midpt[:]/(10**9)
    
    return fluxes

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
zon_mts = pl.mean(wv_zon,axis=0); mer_mts = pl.mean(wv_mer,axis=0)

sheds = ['Am','AfMe','EAA','ArA','ArI','ArP','SOA','SOI','SOP']
F1 = []; F2 = []; F3 = []; F4 = []; F5 = []; F6 = []; F7 = []; F8 = []; F9 = []
for mon in range(filenames.shape[1]):
    F1.append(ShedFluxes(sheddir,sheds[0],lon,lat,zon_mts[mon],mer_mts[mon]))
    F2.append(ShedFluxes(sheddir,sheds[1],lon,lat,zon_mts[mon],mer_mts[mon]))
    F3.append(ShedFluxes(sheddir,sheds[2],lon,lat,zon_mts[mon],mer_mts[mon]))
    F4.append(ShedFluxes(sheddir,sheds[3],lon,lat,zon_mts[mon],mer_mts[mon]))
    F5.append(ShedFluxes(sheddir,sheds[4],lon,lat,zon_mts[mon],mer_mts[mon]))
    F6.append(ShedFluxes(sheddir,sheds[5],lon,lat,zon_mts[mon],mer_mts[mon]))
    F7.append(ShedFluxes(sheddir,sheds[6],lon,lat,zon_mts[mon],mer_mts[mon]))
    F8.append(ShedFluxes(sheddir,sheds[7],lon,lat,zon_mts[mon],mer_mts[mon]))
    F9.append(ShedFluxes(sheddir,sheds[8],lon,lat,zon_mts[mon],mer_mts[mon]))

F1 = pl.asarray(F1); F2 = pl.asarray(F2); F3 = pl.asarray(F3)
F4 = pl.asarray(F4); F5 = pl.asarray(F5); F6 = pl.asarray(F6)
F7 = pl.asarray(F7); F8 = pl.asarray(F8); F9 = pl.asarray(F9)

F1s = pl.sum(F1,axis=1); F2s = pl.sum(F2,axis=1); F3s = pl.sum(F3,axis=1)
F4s = pl.sum(F4,axis=1); F5s = pl.sum(F5,axis=1); F6s = pl.sum(F6,axis=1)
F7s = pl.sum(F7,axis=1); F8s = pl.sum(F8,axis=1); F9s = pl.sum(F9,axis=1)

F = pl.array([F1s,F2s,F3s,F4s,F5s,F6s,F7s,F8s,F9s])

atl = F1s - F2s - F4s + F7s
ind = F2s - F3s + F8s + F5s
pac = F3s - F1s - F6s + F9s
arc = F4s + F5s + F6s
sou = -1*(F7s + F8s + F9s)

D = pl.array([atl,ind,pac,arc,sou])

#f = open(sheddir+'variability/divQ_annmean_sscyc.csv','w')
#f.write('atl ind pac arc sou\n')
#pl.savetxt(f,D.T,fmt='%9.5f',delimiter=' ')
#f.close()

#F1_sns = pl.zeros([4]); F2_sns = pl.zeros([4]); F3_sns = pl.zeros([4])
#F4_sns = pl.zeros([4]); F5_sns = pl.zeros([4]); F6_sns = pl.zeros([4])
#F7_sns = pl.zeros([4]); F8_sns = pl.zeros([4]); F9_sns = pl.zeros([4])

F1_sns = SeasonalCyc(F1);  F2_sns = SeasonalCyc(F2); F3_sns = SeasonalCyc(F3)
F4_sns = SeasonalCyc(F4);  F5_sns = SeasonalCyc(F5); F6_sns = SeasonalCyc(F6)
F7_sns = SeasonalCyc(F7);  F8_sns = SeasonalCyc(F8); F9_sns = SeasonalCyc(F9)

#f = open(sheddir+'variability/seasonal_profs_sop.csv','w')
#pl.savetxt(f,F9_sns.T,fmt='%9.5f')
#f.close()

#S = pl.array([F1_sns,F2_sns,F3_sns,F4_sns,F5_sns,F6_sns,F7_sns,F8_sns,F9_sns])

rlspts1 = pl.genfromtxt(sheddir+'shed_defs/Am_traj_release_new.txt',skip_header=5)
rlspts2 = pl.genfromtxt(sheddir+'shed_defs/AfMe_traj_release_new.txt',skip_header=5)
rlspts3 = pl.genfromtxt(sheddir+'shed_defs/EAA_traj_release_new.txt',skip_header=5)
rlspts4 = pl.genfromtxt(sheddir+'shed_defs/ArA_traj_release_new.txt',skip_header=5)
rlspts5 = pl.genfromtxt(sheddir+'shed_defs/ArI_traj_release_new.txt',skip_header=5)
rlspts6 = pl.genfromtxt(sheddir+'shed_defs/ArP_traj_release_new.txt',skip_header=5)
rlspts7 = pl.genfromtxt(sheddir+'shed_defs/SOA_traj_release_new.txt',skip_header=5)
rlspts8 = pl.genfromtxt(sheddir+'shed_defs/SOI_traj_release_new.txt',skip_header=5)
rlspts9 = pl.genfromtxt(sheddir+'shed_defs/SOP_traj_release_new.txt',skip_header=5)

names = ['DJF','MAM','JJA','SON']
C = ['b','purple','r','dimgrey']#,dtype='|S7')
#pos = pl.arange(4); width=0.75
Q = []

fig,ax = pl.subplots(3,3,figsize=(16,11))

ax1 = pl.subplot(331)
for i in range(4):
    M=ax1.plot(F1_sns[i],color=C[i],lw=2)
    Q.append(M[0])
pl.ylim(-350,100); pl.xlim(0,155); pl.grid(axis='y')
pl.ylabel('kg/m/s',fontsize=22,labelpad=1)#; ax1.text(0.1,-0.35,'Americas',fontsize=17)
pl.title('(a) Americas',fontsize=16)
pl.yticks(fontsize=14)
#pl.subplots_adjust(left=0.14,right=0.91)
a = pl.arange(0,155,20); b = pl.around(rlspts1[a,1],decimals=0)
ax1.set_xticks(a); ax1.set_xticklabels(b.astype('int'),fontsize=14)
ax1.legend(Q,names,loc=3,ncol=2,fontsize=14,columnspacing=0.7)
ax1.tick_params(axis='y',pad=5); ax1.tick_params(axis='x',pad=7)

ax2 = pl.subplot(332)
for i in range(4):
    ax2.plot(F2_sns[i],color=C[i],lw=2)
pl.ylim(-250,100.); pl.xlim(0,214); pl.grid(axis='y')
#pl.ylabel('Sv',fontsize=22,labelpad=1)#; pl.yticks([0,0.05,0.1,0.15,0.2],fontsize=16)
pl.title('(b) Africa',fontsize=16)
a = pl.arange(0,214,30); b = pl.around(rlspts2[a,1],decimals=0)
ax2.set_xticks(a); ax2.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
#pl.subplots_adjust(left=0.14,right=0.91)
ax2.tick_params(axis='y',pad=5); ax2.tick_params(axis='x',pad=7)
pl.xlabel('latitude',fontsize=18)

ax3 = pl.subplot(333)
for i in range(4):
    ax3.plot(F3_sns[i],color=C[i],lw=2)
pl.ylim(-300,400); pl.xlim(0,165); pl.grid(axis='y')
#pl.ylabel('Sv',fontsize=22,labelpad=1)#; pl.yticks([0,0.05,0.1,0.15,0.2],fontsize=16)
pl.title('(c) South-East Asia',fontsize=16)
a = pl.arange(0,165,20); b = pl.around(rlspts3[a,1],decimals=0)
ax3.set_xticks(a); ax3.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
ax3.tick_params(axis='y',pad=5); ax3.tick_params(axis='x',pad=7)

ax4 = pl.subplot(334)
for i in range(4):
    ax4.plot(F4_sns[i],color=C[i],lw=2)
pl.ylim(-100,100); pl.xlim(0,197); pl.grid(axis='y')
pl.ylabel('kg/m/s',fontsize=22,labelpad=1)#; pl.yticks([0,0.05,0.1,0.15,0.2],fontsize=16)
pl.title('(d) Arctic Atlantic',fontsize=16)
a = pl.arange(0,197,30); b = pl.around(rlspts4[a,0],decimals=0)
ax4.set_xticks(a); ax4.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
ax4.tick_params(axis='y',pad=5); ax4.tick_params(axis='x',pad=7)

ax5 = pl.subplot(335)
for i in range(4):
    ax5.plot(F5_sns[i],color=C[i],lw=2)
pl.ylim(-100,100); pl.xlim(0,103); pl.grid(axis='y')
#pl.ylabel('Sv',fontsize=22,labelpad=1)#; pl.yticks([0,0.05,0.1,0.15,0.2],fontsize=16)
pl.title('(e) Arctic Indian',fontsize=16)
a = pl.arange(0,102,20); b = pl.around(rlspts5[a,0],decimals=0)
ax5.set_xticks(a); ax5.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
ax5.tick_params(axis='y',pad=5); ax5.tick_params(axis='x',pad=7)
pl.xlabel('longitude',fontsize=18)

ax6 = pl.subplot(336)
for i in range(4):
    ax6.plot(F6_sns[i],color=C[i],lw=2)
pl.ylim(-100,100); pl.xlim(0,220); pl.grid(axis='y')
#pl.ylabel('Sv',fontsize=22,labelpad=1)#; pl.yticks([0,0.05,0.1,0.15,0.2],fontsize=16)
pl.title('(f) Arctic Pacific',fontsize=16)
a = pl.arange(0,220,30); b = pl.around(rlspts6[a,0],decimals=0)
ax6.set_xticks(a); ax6.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
ax6.tick_params(axis='y',pad=5); ax6.tick_params(axis='x',pad=7)

ax7 = pl.subplot(337)
for i in range(4):
    ax7.plot(F7_sns[i],color=C[i],lw=2)
pl.ylim(-250,150); pl.xlim(0,185); pl.grid(axis='y')
pl.ylabel('kg/m/s',fontsize=22,labelpad=1)
pl.title('(g) Southern Atlantic',fontsize=16)
a = pl.arange(0,185,30); b = pl.around(rlspts7[a,0],decimals=0)
ax7.set_xticks(a); ax7.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
ax7.tick_params(axis='y',pad=5); ax7.tick_params(axis='x',pad=7)

ax8 = pl.subplot(338)
for i in range(4):
    ax8.plot(F8_sns[i],color=C[i],lw=2)
pl.ylim(-150,100); pl.xlim(0,184); pl.grid(axis='y')
#pl.ylabel('Sv',fontsize=22,labelpad=1)
pl.title('(h) Southern Indian',fontsize=16)
a = pl.arange(0,184,30); b = pl.around(rlspts8[a,0],decimals=0)
ax8.set_xticks(a); ax8.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
ax8.tick_params(axis='y',pad=5); ax8.tick_params(axis='x',pad=7)
pl.xlabel('longitude',fontsize=18)

ax9 = pl.subplot(339)
for i in range(4):
    ax9.plot(F9_sns[i],color=C[i],lw=2)
pl.ylim(-100,100.); pl.xlim(0,212); pl.grid(axis='y')
#pl.ylabel('Sv',fontsize=22,labelpad=1)
pl.title('(i) Southern Pacific',fontsize=16)
a = pl.arange(0,210,30); b = pl.around(rlspts9[a,0],decimals=0)
ax9.set_xticks(a); ax9.set_xticklabels(b.astype('int'),fontsize=14)
pl.yticks(fontsize=14)
ax9.tick_params(axis='y',pad=5); ax9.tick_params(axis='x',pad=7)

pl.tight_layout()

#pl.savefig(sheddir+'shedfluxes_bars_sscyc.png')