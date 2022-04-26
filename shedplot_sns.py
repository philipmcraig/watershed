    # -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 18:14:32 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset

def SSCyc(loc):
    """
    """
    sheddir = '/home/users/np838619/Watershed/fluxsplits/'
    flux = pl.genfromtxt(sheddir+'splits_'+loc+'_sscyc.csv',delimiter=',')
    
    MAM = pl.mean(flux[2:5],axis=0)
    JJA = pl.mean(flux[5:8],axis=0)
    SON = pl.mean(flux[8:11],axis=0)
    DJF = (flux[0]+flux[-2]+flux[-1])/3
    
    seasons = pl.array([DJF,MAM,JJA,SON])
    sns_abs = pl.absolute(seasons)
    
    return seasons, sns_abs

def ShedPlot(axx,f1,f2,f3,f4,f5,f6,f7,f8,f9,i,PmE,S):
    """
    """
    m = Basemap(projection='robin',lon_0=-180.,resolution='l')
    m.drawcoastlines(linewidth=0.4,color='lightgray',zorder=1)
    m.drawcountries(color='lightgray',zorder=1)
    
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
    afr_lon = [26,50,16]; afr_lat = [27,2.5,-29]; x2,y2 = m(afr_lon, afr_lat)
    eaa_lon = [87,123]; eaa_lat = [16.5,-2]; x3,y3 = m(eaa_lon, eaa_lat)
    ara_lon = [266,9,36]; ara_lat = [70,72,49]; x4,y4 = m(ara_lon, ara_lat)
    ari_lon = [44,65,88]; ari_lat = [33,34,26]; x5,y5 = m(ari_lon, ari_lat)
    arp_lon = [100,130,219]; arp_lat = [40,60,57]; x6,y6 = m(arp_lon,arp_lat)
    soa_lon = [308,338]; soa_lat = [-12,-30]; x7,y7 = m(soa_lon,soa_lat)
    soi_lon = [76,122]; soi_lat = [-30,-21]; x8,y8 = m(soi_lon,soi_lat)
    sop_lon = [145,220,286]; sop_lat = [-32,-32,-26]; x9,y9 = m(sop_lon,sop_lat)
    
    pl.arrow(x1[0],y1[0],1500000,0,fc="none", ec="b", linewidth = 4, head_width=600000,
         head_length=600000,width=9000,head_starts_at_zero=False)
    axx.annotate('{:.2f}'.format(f1[0]),xy=(0.61,0.685),xycoords='axes fraction',color='b',
                size=15)
    pl.arrow(x1[1],y1[1],-1500000,0,fc="none", ec="b", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f1[1]),xy=(0.68,0.485),xycoords='axes fraction',color='b',
                size=15)
    
    pl.arrow(x2[0],y2[0],990000,-250000,fc="none", ec="k", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f2[0]),xy=(0.03,0.66),xycoords='axes fraction',color='k',
                size=15)
    pl.arrow(x2[1],y2[1],-2000000,0,fc="none", ec="k", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f2[1]),xy=(0.01,0.46),xycoords='axes fraction',color='k',
                size=15)
    pl.arrow(x2[2],y2[2],2000000,0,fc="none", ec="k", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f2[2]),xy=(0.095,0.35),xycoords='axes fraction',color='k',
                size=15)

    if i == 0:
        P = -1; x3a = x3[0] + 2000000
    else:
        P = 1; x3a = x3[0]#; y3a = y3[0]
    pl.arrow(x3a,y3[0],P*1700000,0,fc="none", ec="g", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f3[0]),xy=(0.33,0.58),xycoords='axes fraction',color='g',
                size=15)
    if i == 2:
        P = 1; x3b = x3[1] - 1000000; y3b = y3[1] - 1000000
    else:
        P = -1; x3b = x3[1]; y3b = y3[1]
    pl.arrow(x3b,y3b,P*1000000,P*1000000,fc="none", ec="g", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f3[1]),xy=(0.36,0.47),xycoords='axes fraction',color='g',
                size=15)

    pl.arrow(x4[0],y4[0],450000,-850000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f4[0]),xy=(0.71,0.81),xycoords='axes fraction',color='r',
                size=15)
    pl.arrow(x4[1],y4[1],1600000,0,fc="none", ec="r", linewidth = 4, head_width=650000,
             head_length=650000,width=10000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f4[1]),xy=(0.24,0.92),xycoords='axes fraction',color='r',
                size=15)
    pl.arrow(x4[2],y4[2],1100000,-150000,fc="none", ec="r", linewidth = 4, head_width=650000,
             head_length=650000,width=10000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f4[2]),xy=(0.08,0.77),xycoords='axes fraction',color='r',
                size=15)

    pl.arrow(x5[0],y5[0],1100000,250000,fc="none", ec="r", linewidth = 2, head_width=500000,
             head_length=700000,width=9000,head_starts_at_zero=False)
#    pl.annotate('{:.2f}'.format(f5[0]),xy=(0.144,0.666),xycoords='axes fraction',color='r',
#                size=18,zorder=10)
    pl.arrow(x5[1],y5[1],500000,-700000,fc="none", ec="r", linewidth = 2, head_width=500000,
             head_length=700000,width=9000,head_starts_at_zero=False)
#    pl.annotate('{:.2f}'.format(f5[1]),xy=(0.19,0.6),xycoords='axes fraction',color='r',
#                size=18,zorder=10)
    pl.arrow(x5[2],y5[2],0,700000,fc="none", ec="r", linewidth = 2, head_width=400000,
             head_length=600000,width=9000,head_starts_at_zero=False)
#    pl.annotate('{:.2f}'.format(f5[2]),xy=(0.25,0.735),xycoords='axes fraction',color='r',
#                size=18,zorder=10)

    pl.arrow(x6[0],y6[0],900000,-450000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f6[0]),xy=(0.32,0.65),xycoords='axes fraction',color='r',
                size=15,zorder=10)
    pl.arrow(x6[1],y6[1],900000,-510000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f6[1]),xy=(0.43,0.77),xycoords='axes fraction',color='r',
                size=15,zorder=10)
    pl.arrow(x6[2],y6[2],900000,510000,fc="none", ec="r", linewidth = 4, head_width=600000,
             head_length=700000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f6[2]),xy=(0.55,0.79),xycoords='axes fraction',color='r',
                size=15,zorder=10)
    
    if  i == 2:
        P = 1; x7a = x7[0] - 450000; y7a = y7[0] - 980000
    else:
        P = -1; x7a = x7[0]; y7a = y7[0]
    pl.arrow(x7a,y7a,P*450000,P*750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f7[0]),xy=(0.865,0.42),xycoords='axes fraction',
                color='darkgoldenrod',size=15)
    pl.arrow(x7[1],y7[1],0,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f7[1]),xy=(0.90,0.33),xycoords='axes fraction',
                color='darkgoldenrod',size=15)

    pl.arrow(x8[0],y8[0],0,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f8[0]),xy=(0.19,0.33),xycoords='axes fraction',
                color='darkgoldenrod',size=15)
    pl.arrow(x8[1],y8[1],450000,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f8[1]),xy=(0.28,0.34),xycoords='axes fraction',
                color='darkgoldenrod',size=15)

    if i == 1:
        P = -1; x9a = x9[0] + 650000; y9a = y9[0] + 950000
    else:
        P = 1; x9a = x9[0]; y9a = y9[0]
    pl.arrow(x9a,y9a,P*650000,P*650000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f9[0]),xy=(0.445,0.35),xycoords='axes fraction',
                color='darkgoldenrod',size=15)
    pl.arrow(x9[1],y9[1],0,-750000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f9[1]),xy=(0.58,0.33),xycoords='axes fraction',
                color='darkgoldenrod',size=15)
    pl.arrow(x9[2],y9[2],650000,0,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=600000,
             head_length=600000,width=9000,head_starts_at_zero=False)
    pl.annotate('{:.2f}'.format(f9[2]),xy=(0.725,0.33),xycoords='axes fraction',
                color='darkgoldenrod',size=15)
    
    a,b = m(180,15)
    pl.text(a,b,'{:.2f}'.format(PmE[2]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    a,b = m(310,15)
    pl.text(a,b,'{:.2f}'.format(PmE[0]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    a,b = m(55,0)
    pl.text(a,b,'{:.2f}'.format(PmE[1]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    a,b = m(90,70)
    pl.text(a,b,'{:.2f}'.format(PmE[3]),bbox={'facecolor':'white', 'alpha':1.0, 'pad':5},fontsize=23)
    a,b = m(140,-55)
    pl.text(a,b,'{:.2f}'.format(PmE[4]),bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
    
    a,b = m(180,-2)
    pl.text(a,b,'('+str(S[2])+')',fontsize=17)
    a,b = m(310,-2)
    pl.text(a,b,'('+str(S[0])+')',fontsize=17)
    a,b = m(55,-17)
    pl.text(a,b,'('+str(S[1])+')',fontsize=17)
    a,b = m(150,75)
    pl.text(a,b,'('+'{:.2f}'.format(S[3])+')',fontsize=17)
    a,b = m(190,-55)
    pl.text(a,b,'('+str(S[4])+')',fontsize=17)
    
    a,b = m(50,-75)
    pl.text(a,b,'Arctic Indian sector fluxes: '+'{:.2f}'.format(f5[0])+', '+
                '{:.2f}'.format(f5[1])+', '+'{:.2f}'.format(f5[2]),
            bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=14,color='r')


pl.close('all')

amr_sns, amr_abs = SSCyc('Am'); afr_sns, afr_abs = SSCyc('AfMe')
eaa_sns, eaa_abs = SSCyc('EAA'); ara_sns, ara_abs = SSCyc('ArA')
ari_sns, ari_abs = SSCyc('ArI'); arp_sns, arp_abs = SSCyc('ArP')
soa_sns, soa_abs = SSCyc('SOA'); soi_sns, soi_abs = SSCyc('SOI')
sop_sns, sop_abs = SSCyc('SOP')

atl_sns = amr_sns[:,0] + amr_sns[:,1] - ara_sns[:,0] - ara_sns[:,1] - \
        ara_sns[:,2] - afr_sns[:,0] - afr_sns[:,1] - afr_sns[:,2] + \
        soa_sns[:,0] + soa_sns[:,1]
ind_sns = -ari_sns[:,0] - ari_sns[:,1] - ari_sns[:,2] - eaa_sns[:,0] - \
        eaa_sns[:,1] + soi_sns[:,0] + soi_sns[:,1] + afr_sns[:,0] + \
        afr_sns[:,1] + afr_sns[:,2]
pac_sns = -arp_sns[:,0] - arp_sns[:,1] - arp_sns[:,2] - amr_sns[:,0] - \
        amr_sns[:,1] + sop_sns[:,0] + sop_sns[:,1] + sop_sns[:,2] + \
        eaa_sns[:,0] + eaa_sns[:,1]
arc_sns = ara_sns[:,0] + ara_sns[:,1] + ara_sns[:,2] + ari_sns[:,0] + \
        ari_sns[:,1] + ari_sns[:,2] + arp_sns[:,0] + arp_sns[:,1] + \
        arp_sns[:,2]
sou_sns = -soa_sns[:,0] - soa_sns[:,1] - soi_sns[:,0] - soi_sns[:,1] - \
        sop_sns[:,0] - sop_sns[:,1] - sop_sns[:,2]

PmE = pl.array([atl_sns,ind_sns,pac_sns,arc_sns,sou_sns])
PmE = PmE.T

storage = pl.array([[0.01,-0.07,-0.03,0.08], # Atlantic
                    [-0.01,0.01,0.01,-0.01], # Indian
                    [-0.01,-0.06,-0.02,0.08], # Pacific
                    [-0.00,-0.04,-0.01,0.05], # Arctic
                    [-0.02,0.06,0.02,-0.05] # Southern
                    ])
storage = storage.T

fig,ax = pl.subplots(2,2,figsize=(23.5,12))

ax1 = pl.subplot(221); i = 0
ShedPlot(ax1,amr_abs[i],afr_abs[i],eaa_abs[i],ara_abs[i],ari_abs[i],arp_abs[i],
         soa_abs[i],soi_abs[i],sop_abs[i],i,PmE[i],storage[i])
pl.title('(a) DJF',fontsize=22)

ax2 = pl.subplot(222); i = 1
ShedPlot(ax2,amr_abs[i],afr_abs[i],eaa_abs[i],ara_abs[i],ari_abs[i],arp_abs[i],
         soa_abs[i],soi_abs[i],sop_abs[i],i,PmE[i],storage[i])
pl.title('(b) MAM',fontsize=22)

ax3 = pl.subplot(223); i = 2
ShedPlot(ax3,amr_abs[i],afr_abs[i],eaa_abs[i],ara_abs[i],ari_abs[i],arp_abs[i],
         soa_abs[i],soi_abs[i],sop_abs[i],i,PmE[i],storage[i])
pl.title('(c) JJA',fontsize=22)

ax4 = pl.subplot(224); i = 3
ShedPlot(ax4,amr_abs[i],afr_abs[i],eaa_abs[i],ara_abs[i],ari_abs[i],arp_abs[i],
         soa_abs[i],soi_abs[i],sop_abs[i],i,PmE[i],storage[i])
pl.title('(d) SON',fontsize=22)

#fluxes = pl.array([0.23,0.17,0.16,0.98,0.17,0.01,0.05,0.07,0.29,0.27,0.38,0.09,0.12,0.07])


#a,b = m(180,15)
#pl.text(a,b,'0.00',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10},fontsize=30)
#a,b = m(320,15)
#pl.text(a,b,'-0.47',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10},fontsize=30)
#a,b = m(70,-15)
#pl.text(a,b,'-0.64',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10},fontsize=30)
#a,b = m(100,70)
#pl.text(a,b,'0.16',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=30)
#a,b = m(180,-55)
#pl.text(a,b,'0.95',bbox={'facecolor':'white', 'alpha':0.5, 'pad':10},fontsize=30)

pl.subplots_adjust(top=0.95,bottom=0.05,left=0.05,right=0.95,hspace=0.12,wspace=0.00)
#pl.tight_layout()