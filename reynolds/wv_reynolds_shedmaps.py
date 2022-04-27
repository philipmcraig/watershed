# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 11:57:43 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset
import os

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
    flux_uv = pl.zeros_like(rlspts)
    for i in range(len(flux_uv)):
        #print i
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
        #print i
        segno = labs[i,0] - 1
        FdotN[i] = pl.dot(midflux[i],nhat[segno])
        #if segno not in labs:
        #    pass
        #else:
            #print i, segno
            #FdotN[i] = pl.dot(midflux[i],nhat[segno])
    
    return FdotN

#def TrajSegLabel(loc):
#    """Function to assign each trajectory release point a label referring to which
#	segment of the watershed it is from.
#
#	Args:
#		loc (string): stem of continental watershed e.g. NCA = North/Central America
#	
#	Returns:
#             seglab (array): segment labels (integers) of each release point
#             rlspts (array): trajectory release points
#    """
#    sheddir = '/home/np838619/Watershed/shed_defs/'
#    
#    #endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
#    #endpts = pl.asarray(endlist[5:],dtype='float')
#    
#    rlslist = ReadTxtFile(sheddir + loc + '_traj_release.txt')
#    rlspts = pl.asarray(rlslist[5:],dtype='float')
#    
#    seglab = [1] # first release point has to be on the first segment
#    #segnos = pl.linspace(1,40,40)
#    count = 1
#    
#    for rls in range(1,rlspts.shape[0]):
#        if rlspts[rls,0] == rlspts[rls-1,0] and rlspts[rls,1] == rlspts[rls-1,1]:
#            count = count + 1
#            seglab.append(count)
#        else:
#            count = count
#            seglab.append(count)
#    seglab = pl.asarray(seglab)
#    
#    return seglab, rlspts 

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

def SplitCalcs(splitlocs,sheddir,loc,lon,lat,zon,mer):
    """
    """
    F = ShedFluxes(sheddir,loc,lon,lat,zon,mer)
    
    if len(splitlocs) == 1.:
        Fs1 = F[:splitlocs[0]+1]; Fs1_tot = pl.sum(Fs1)
        Fs2 = F[splitlocs[0]+1:]; Fs2_tot = pl.sum(Fs2)
        splits = pl.array([Fs1_tot,Fs2_tot])
    elif len(splitlocs) == 2.:
        Fs1 = F[:splitlocs[0]+1]; Fs1_tot = pl.sum(Fs1)
        Fs2 = F[splitlocs[0]+1:splitlocs[1]+1]; Fs2_tot = pl.sum(Fs2)
        Fs3 = F[splitlocs[1]+1:]; Fs3_tot = pl.sum(Fs3)
        splits = pl.array([Fs1_tot,Fs2_tot,Fs3_tot])

    return splits

def ShedPlot(axx,f1,f2,f3,f4,f5,f6,f7,f8,f9,PmE):
    """
    """
    m = Basemap(projection='robin',lon_0=-180.,resolution='l')
    m.drawcoastlines(linewidth=0.4,color='lightgray',zorder=1)
    #m.drawcountries(color='lightgray')
    
    sheddir = '/home/users/np838619/Watershed/shed_defs/'
    filex = '_traj_release_new.txt'
    
    Am = pl.genfromtxt(sheddir + 'Am' + filex,skip_header=5)
    AfMe = pl.genfromtxt(sheddir + 'AfMe' + filex,skip_header=5)
    EAA = pl.genfromtxt(sheddir + 'EAA' + filex,skip_header=5)
    Ar = pl.genfromtxt(sheddir + 'Ar' + filex,skip_header=5)
    SO = pl.genfromtxt(sheddir + 'SO' + filex,skip_header=5)
    NAs = pl.genfromtxt(sheddir + 'NAs' + filex,skip_header=5)
    
    lw = 3
    m.plot(Am[:,0],Am[:,1],latlon=True,color='b',linewidth=lw)
    m.plot(AfMe[:,0],AfMe[:,1],latlon=True,color='k',linewidth=lw)
    m.plot(EAA[:,0],EAA[:,1],latlon=True,color='g',linewidth=lw)
    m.plot(NAs[:,0],NAs[:,1],latlon=True,color='maroon',ls='--',linewidth=2)
    
    F = pl.zeros_like(SO); F[:,0], F[:,1] = m(SO[:,0],SO[:,1])
    m.plot(F[:246,0], F[:246,1],color='darkgoldenrod',linewidth=lw)
    m.plot(F[246:555,0], F[246:555,1],color='darkgoldenrod',linewidth=lw)
    m.plot(F[555:,0],F[555:,1],color='darkgoldenrod',linewidth=lw)
    
    m.plot(Ar[:157,0],Ar[:157,1],latlon=True,linewidth=lw,color='r',zorder=50)
    m.plot(Ar[157:,0],Ar[157:,1],latlon=True,linewidth=lw,color='r',zorder=50)
    
    isect1 = pl.where(Am[0]==Ar); isect1=(Ar[isect1[0][0],0],Ar[isect1[0][0],1])
    isect2 = pl.where(AfMe[0]==Ar); isect2=(Ar[isect2[0][0],0],Ar[isect2[0][0],1])
    isect3 = pl.where(EAA[0]==Ar); isect3=(Ar[isect3[0][0],0],Ar[isect3[0][0],1])
    m.plot(isect1[0],isect1[1],'r.',markersize=13,mew=3,latlon=True)
    m.plot(isect2[0],isect2[1],'r.',markersize=13,mew=3,latlon=True)
    m.plot(isect3[0],isect3[1],'r.',markersize=13,mew=3,latlon=True)
    
    isect4 = pl.where(Am[-1]==SO); isect4=(SO[isect4[0][0],0],SO[isect4[0][0],1])
    isect5 = pl.where(AfMe[-1]==SO); isect5=(SO[isect5[0][0],0],SO[isect5[0][0],1])
    isect6 = pl.where(EAA[-1]==SO); isect6=(SO[isect6[0][0],0],SO[isect6[0][0],1])
    m.plot(isect4[0],isect4[1],color='darkgoldenrod',marker='.',markersize=13,mew=3,latlon=True)
    m.plot(isect5[0],isect5[1],color='darkgoldenrod',marker='.',markersize=13,mew=3,latlon=True)
    m.plot(isect6[0],isect6[1],color='darkgoldenrod',marker='.',markersize=13,mew=3,latlon=True)
    
    isect7 = pl.where(NAs[0]==Ar); isect7=(Ar[isect7[0][0],0],Ar[isect7[0][0],1])
    isect8 = pl.where(NAs[-1]==Ar); isect8=(Ar[isect8[0][0],0],Ar[isect8[0][0],1])
    m.plot(isect7[0],isect7[1],'rx',markersize=13,mew=3,latlon=True)
    m.plot(isect8[0],isect8[1],'rx',markersize=13,mew=3,latlon=True)
    
    m.plot(Am[67,0]+360,Am[67,1],'bx',markersize=13,mew=3,latlon=True) # amr split
    m.plot(AfMe[51,0],AfMe[51,1],'kx',markersize=13,mew=3,latlon=True) # afr split 1
    m.plot(AfMe[182,0],AfMe[182,1],'kx',markersize=13,mew=3,latlon=True) # afr split 2
    m.plot(EAA[51,0],EAA[51,1],'gx',markersize=13,mew=3,latlon=True) # eaa split
    m.plot(Ar[121,0],Ar[121,1],'rx',markersize=13,mew=3,latlon=True) # ara aplit
    m.plot(Ar[291,0],Ar[291,1],'rx',markersize=13,mew=3,latlon=True) # ari split 1
    m.plot(Ar[325,0],Ar[325,1],'rx',markersize=13,mew=3,latlon=True) # arp split 2
    m.plot(Ar[431,0],Ar[431,1],'rx',markersize=13,mew=3,latlon=True) # arp split 1
    m.plot(Ar[517,0],Ar[517,1],'rx',markersize=13,mew=3,latlon=True) # arp split 2
    m.plot(SO[489,0],SO[489,1],'x',mec='darkgoldenrod',markersize=13,mew=3,latlon=True) # soa split
    m.plot(SO[489,0],SO[489,1],'x',mec='darkgoldenrod',markersize=13,mew=3,latlon=True) # soi split
    m.plot(SO[112,0],SO[112,1],'x',mec='darkgoldenrod',markersize=13,mew=3,latlon=True) # soi split
    m.plot(SO[213,0],SO[213,1],'x',mec='darkgoldenrod',markersize=13,mew=3,latlon=True) # sop split 1
    m.plot(SO[363,0],SO[363,1],'x',mec='darkgoldenrod',markersize=13,mew=3,latlon=True) # sop split 2

    lon = pl.array([-85,55,80,-97,66,188,-22,80,-130,77,42,101])
    lat = pl.array([18.4,2.5,6.1,60,23,60,-30,-30,-30,45,48,38])
    x,y = m(lon, lat)
    
    amr_lon = [245,290]; amr_lat = [34,0]; x1,y1 = m(amr_lon, amr_lat)
    afr_lon = [26,55,16]; afr_lat = [27,2.5,-29]; x2,y2 = m(afr_lon, afr_lat)
    eaa_lon = [85,123]; eaa_lat = [16.5,-2]; x3,y3 = m(eaa_lon, eaa_lat)
    ara_lon = [266,9,36]; ara_lat = [70,72,49]; x4,y4 = m(ara_lon, ara_lat)
    ari_lon = [44,65,88]; ari_lat = [33,34,26]; x5,y5 = m(ari_lon, ari_lat)
    arp_lon = [100,130,219]; arp_lat = [40,60,57]; x6,y6 = m(arp_lon,arp_lat)
    soa_lon = [308,338]; soa_lat = [-12,-30]; x7,y7 = m(soa_lon,soa_lat)
    soi_lon = [76,122]; soi_lat = [-30,-21]; x8,y8 = m(soi_lon,soi_lat)
    sop_lon = [145,220,286]; sop_lat = [-32,-32,-26]; x9,y9 = m(sop_lon,sop_lat)
    
    pl.arrow(x1[0],y1[0],1500000,0,fc="none", ec="b", linewidth = 4, head_width=600000,
         head_length=600000,width=9000,head_starts_at_zero=False)
    axx.annotate('{:.2f}'.format(pl.absolute(f1[0])),xy=(0.61,0.685),
                 xycoords='axes fraction',color='b',size=20)
    if f1[1] > 0:
        P = 1; x1a = x1[1] - 1500000; y1a = y1[1]
    else:
        P = -1; x1a = x1[1]; y1a = y1[1]
    pl.arrow(x1a,y1a,P*1500000,0,fc="none", ec="b", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f1[1])),xy=(0.68,0.485),
                xycoords='axes fraction',color='b',size=20)

    if f2[0] > 0:
        P = 1; x2a = x2[0]; y2a = y2[0]
    else:
        P = -1; x2a = x2[0] + 1500000; y2a = y2[0] + 600000
    pl.arrow(x2a,y2a,P*990000,-P*250000,fc="none", ec="k", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f2[0])),xy=(0.02,0.66),
                xycoords='axes fraction',color='k',size=20)
    if f2[1] < 0:
        P = 1; x2b= x2[1]; y2b = y2[1]
    else:
        P = -1; x2b = x2[1] - 3000000; y2b = y2[1]
    pl.arrow(x2b,y2b,-P*2500000,0,fc="none", ec="k", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f2[1])),xy=(0.11,0.455),
                xycoords='axes fraction',color='k',size=20)
    pl.arrow(x2[2],y2[2],2000000,0,fc="none", ec="k", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f2[2])),xy=(0.095,0.35),
                xycoords='axes fraction',color='k',size=20)

    pl.arrow(x3[0],y3[0],1700000,0,fc="none", ec="g", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f3[0])),xy=(0.21,0.54),
                xycoords='axes fraction',color='g',size=20)
    pl.arrow(x3[1],y3[1],-1000000,-1000000,fc="none", ec="g", linewidth = 4,
             head_width=600000,head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f3[1])),xy=(0.35,0.47),
                xycoords='axes fraction',color='g',size=20)

    if f4[0] < 0:
        P = 1; x4a = x4[0]; y4a = y4[0]
    else:
        P = -1; x4a = x4[0] + 500000; y4a = y4[0] - 1000000
    pl.arrow(x4a,y4a,P*450000,-P*850000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f4[0])),xy=(0.70,0.83),
                xycoords='axes fraction',color='r',size=20)
    pl.arrow(x4[1],y4[1],1600000,0,fc="none", ec="r", linewidth = 4, head_width=650000,
             head_length=650000,width=10000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f4[1])),xy=(0.24,0.92),
                xycoords='axes fraction',color='r',size=20)
    pl.arrow(x4[2],y4[2],1100000,-150000,fc="none", ec="r", linewidth = 4, head_width=650000,
             head_length=650000,width=10000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f4[2])),xy=(0.08,0.77),
                xycoords='axes fraction',color='r',size=20)

    pl.arrow(x5[0],y5[0],1100000,400000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
#    pl.annotate('{:.2f}'.format(pl.absolute(f5[0])),xy=(0.144,0.666),
#                xycoords='axes fraction',color='r',size=18,zorder=10)
    if f5[1] < 0:
        P = -1; x5a = x5[1]; y5a = y5[1]
    else:
        P = 1; x5a = x5[1]; y5a = y5[1] - 900000
    pl.arrow(x5a,y5a,0,P*700000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
#    pl.annotate('{:.2f}'.format(pl.absolute(f5[1])),xy=(0.19,0.6),
#                xycoords='axes fraction',color='r',size=18,zorder=10)
    pl.arrow(x5[2],y5[2],0,700000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
#    pl.annotate('{:.2f}'.format(pl.absolute(f5[2])),xy=(0.25,0.735),
#                xycoords='axes fraction',color='r',size=18,zorder=10)

    if f6[0] < 0:
        P = 1; x6a = x6[0]; y6a = y6[0]
    else:
        P = -1; x6a = x6[0] + 900000; y6a = y6[0] - 450000
    pl.arrow(x6a,y6a,P*900000,-P*450000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f6[0])),xy=(0.32,0.65),
                xycoords='axes fraction',color='r',size=20,zorder=10)
    if f6[1] < 0:
        P = 1; x6b = x6[1]; y6b = y6[1]
    else:
        P = -1; x6b = x6[1] + 900000; y6b = y6[1] - 510000
    pl.arrow(x6b,y6b,P*900000,-P*510000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f6[1])),xy=(0.43,0.83),
                xycoords='axes fraction',color='r',size=20,zorder=10)
    pl.arrow(x6[2],y6[2],900000,510000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f6[2])),xy=(0.55,0.79),
                xycoords='axes fraction',color='r',size=20,zorder=10)
    
    if  f7[0] < 0:
        P = 1; x7a = x7[0]; y7a = y7[0]#   
    else:
        P = -1; x7a = x7[0] - 450000; y7a = y7[0] - 980000
    pl.arrow(x7a,y7a,-P*450000,-P*750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f7[0])),xy=(0.865,0.42),
                xycoords='axes fraction',color='darkgoldenrod',size=18)
    pl.arrow(x7[1],y7[1],0,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f7[1])),xy=(0.90,0.33),
                xycoords='axes fraction',color='darkgoldenrod',size=20)

    pl.arrow(x8[0],y8[0],0,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f8[0])),xy=(0.19,0.33),
                xycoords='axes fraction',color='darkgoldenrod',size=20)
    pl.arrow(x8[1],y8[1],450000,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f8[1])),xy=(0.28,0.34),
                xycoords='axes fraction',color='darkgoldenrod',size=18)

    if f9[0] > 0:
        P = 1; x9a = x9[0]; y9a = y9[0]
    else:
        P = -1; x9a = x9[0] + 650000; y9a = y9[0]  + 950000
    pl.arrow(x9a,y9a,P*650000,P*650000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f9[0])),xy=(0.44,0.33),
                xycoords='axes fraction',color='darkgoldenrod',size=18)
    pl.arrow(x9[1],y9[1],0,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f9[1])),xy=(0.58,0.33),
                xycoords='axes fraction',color='darkgoldenrod',size=20)
    if f9[2] < 0:
        P = 1; x9b = x9[2]; y9b = y9[2]
    else:
        P = -1; x9b = x9[2] + 1000000; y9b = y9[2]
    pl.arrow(x9b,y9b,P*650000,0,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(pl.absolute(f9[2])),xy=(0.725,0.33),
                xycoords='axes fraction',color='darkgoldenrod',size=18)
    
    a,b = m(180,15)
    pl.text(a,b,'{:.2f}'.format(PmE[2]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    a,b = m(310,15)
    pl.text(a,b,'{:.2f}'.format(PmE[0]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    a,b = m(65,-13)
    pl.text(a,b,'{:.2f}'.format(PmE[1]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    a,b = m(100,70)
    pl.text(a,b,'{:.2f}'.format(PmE[3]),bbox={'facecolor':'white', 'alpha':1.0, 'pad':5},fontsize=23)
    a,b = m(180,-55)
    pl.text(a,b,'{:.2f}'.format(PmE[4]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    
    a,b = m(50,-75)
    pl.text(a,b,'Arctic Indian sector fluxes: '+'{:.2f}'.format(pl.absolute(f5[0]))+', '+
                '{:.2f}'.format(pl.absolute(f5[1]))+', '+'{:.2f}'.format(pl.absolute(f5[2])),
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=14,color='r')

pl.close('all')

exec(open('/home/users/np838619/Trajectory/trajfuncs.py').read())

sheddir = '/home/users/np838619/Watershed/'
reyndir = sheddir + 'reynolds/'

ncfile = Dataset(reyndir+'reyn_av_quv_7914_new.nc','r')
#eralat = ncfile.variables['lat'][:]
#eralon = ncfile.variables['lon'][:]
qu = ncfile.variables['Vertically-integrated zonal moisture flux (full field)'][:]
qu_mf = ncfile.variables['Vertically-integrated zonal moisture flux (mean flow)'][:]
qu_ed = ncfile.variables['Vertically-integrated zonal moisture flux (eddies)'][:]
qv = ncfile.variables['Vertically-integrated meridional moisture flux (full field)'][:]
qv_mf = ncfile.variables['Vertically-integrated meridional moisture flux (mean flow)'][:]
qv_ed = ncfile.variables['Vertically-integrated meridional moisture flux (eddies)'][:]
ncfile.close()

ncfile = Dataset(sheddir+'wvfluxes_7914.nc','r')
eralon = ncfile.variables['lon'][:]
eralat = ncfile.variables['lat'][:]
ncfile.close()

splits = [[67],
          [51,182],
          [51],
          [59,137],
          [33,67],
          [71,157],
          [95],
          [112],
          [30,180]]
locs = ['Am','AfMe','EAA','ArA','ArI','ArP','SOA','SOI','SOP']

flx_mf = []; flx_ed = []

for i in range(len(splits)):
    A = SplitCalcs(splits[i],sheddir,locs[i],eralon,eralat,qu_mf,qv_mf)
    B = SplitCalcs(splits[i],sheddir,locs[i],eralon,eralat,qu_ed,qv_ed)
    flx_mf.append(A); flx_ed.append(B)

atl_mf = -pl.sum(flx_mf[3]) - pl.sum(flx_mf[1]) + pl.sum(flx_mf[6]) + pl.sum(flx_mf[0])
atl_ed = -pl.sum(flx_ed[3]) - pl.sum(flx_ed[1]) + pl.sum(flx_ed[6]) + pl.sum(flx_ed[0])

ind_mf = -pl.sum(flx_mf[4]) - pl.sum(flx_mf[2]) + pl.sum(flx_mf[7]) + pl.sum(flx_mf[1])
ind_ed = -pl.sum(flx_ed[4]) - pl.sum(flx_ed[2]) + pl.sum(flx_ed[7]) + pl.sum(flx_ed[1])

pac_mf = -pl.sum(flx_mf[5]) - pl.sum(flx_mf[0]) + pl.sum(flx_mf[-1]) + pl.sum(flx_mf[2])
pac_ed = -pl.sum(flx_ed[5]) - pl.sum(flx_ed[0]) + pl.sum(flx_ed[-1]) + pl.sum(flx_ed[2])

arc_mf = pl.sum(flx_mf[3]) + pl.sum(flx_mf[4]) + pl.sum(flx_mf[5])
arc_ed = pl.sum(flx_ed[3]) + pl.sum(flx_ed[4]) + pl.sum(flx_ed[5])

sou_mf = -pl.sum(flx_mf[6]) - pl.sum(flx_mf[7]) - pl.sum(flx_mf[8])
sou_ed = -pl.sum(flx_ed[6]) - pl.sum(flx_ed[7]) - pl.sum(flx_ed[8])

PmE = pl.array([[atl_mf,atl_ed],
                [ind_mf,ind_ed],
                [pac_mf,pac_ed],
                [arc_mf,arc_ed],
                [sou_mf,sou_ed]])

fig, ax = pl.subplots(2,1,figsize=(12,12))

ax1 = pl.subplot(211)
ShedPlot(ax1,flx_mf[0],flx_mf[1],flx_mf[2],flx_mf[3],flx_mf[4],flx_mf[5],
         flx_mf[6],flx_mf[7],flx_mf[8],PmE[:,0])
pl.title('(a) $\overline{q}\,\overline{\mathbf{v}}\cdot\mathbf{\hat{n}}$',fontsize=22)

ax2 = pl.subplot(212)
ShedPlot(ax2,flx_ed[0],flx_ed[1],flx_ed[2],flx_ed[3],flx_ed[4],flx_ed[5],
         flx_ed[6],flx_ed[7],flx_ed[8],PmE[:,1])
pl.title('(b) $\overline{q\\backprime \mathbf{v}\\backprime}\cdot\mathbf{\hat{n}}$',
         fontsize=22)
 
pl.tight_layout()
#pl.savefig(reyndir+'split_fluxes_reyn.png')