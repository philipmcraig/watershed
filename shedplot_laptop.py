    # -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 18:14:32 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset


pl.close('all')
pl.figure(figsize=(23.5,13))
m = Basemap(projection='robin',lon_0=-180.,resolution='l')
            #llcrnrlat=-80.,urcrnrlon=360,urcrnrlat=80.)
m.drawcoastlines(linewidth=0.5,color='lightgray',zorder=1)
m.drawcountries(color='lightgray',zorder=1)#; pl.tight_layout()

sheddir = '/home/users/qx911590/np838619/Watershed/shed_defs/'
filex = '_traj_release_new.txt'
file_inv = 'ggis'; fname3 = file_inv + '198901010000.nc'

nc3 = Dataset('/home/users/qx911590/np838619/Watershed/'+fname3,'r')
geop = nc3.variables['Z'][:]
nc3.close()
geop = pl.squeeze(geop)
orog = geop/9.81

Am = pl.genfromtxt(sheddir + 'Am' + filex,skip_header=5)
AfMe = pl.genfromtxt(sheddir + 'AfMe' + filex,skip_header=5)
EAA = pl.genfromtxt(sheddir + 'EAA' + filex,skip_header=5)
#LP = pl.genfromtxt(sheddir + 'LP' + filex,skip_header=5)
#Oz = pl.genfromtxt(sheddir + 'Oz' + filex,skip_header=5)
Ar = pl.genfromtxt(sheddir + 'Ar' + filex,skip_header=5)
SO = pl.genfromtxt(sheddir + 'SO' + filex,skip_header=5)
NAs = pl.genfromtxt(sheddir + 'NAs' + filex,skip_header=5)

#NAs[-6:,0] = NAs[-6:,0] + 360.
lw = 3
m.plot(Am[:,0],Am[:,1],latlon=True,color='b',linewidth=lw)
#m.plot(SAm[:,0],SAm[:,1],latlon=True,color='r',linewidth=lw)
m.plot(AfMe[:,0],AfMe[:,1],latlon=True,color='k',linewidth=lw)
#m.plot(Eu[:,0],Eu[:,1],latlon=True,color='darkgreen',linewidth=lw)
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
m.plot(isect1[0],isect1[1],'r.',markersize=15,mew=5,latlon=True)
m.plot(isect2[0],isect2[1],'r.',markersize=15,mew=5,latlon=True)
m.plot(isect3[0],isect3[1],'r.',markersize=15,mew=5,latlon=True)

isect4 = pl.where(Am[-1]==SO); isect4=(SO[isect4[0][0],0],SO[isect4[0][0],1])
isect5 = pl.where(AfMe[-1]==SO); isect5=(SO[isect5[0][0],0],SO[isect5[0][0],1])
isect6 = pl.where(EAA[-1]==SO); isect6=(SO[isect6[0][0],0],SO[isect6[0][0],1])
m.plot(isect4[0],isect4[1],color='darkgoldenrod',marker='.',markersize=15,mew=5,latlon=True)
m.plot(isect5[0],isect5[1],color='darkgoldenrod',marker='.',markersize=15,mew=5,latlon=True)
m.plot(isect6[0],isect6[1],color='darkgoldenrod',marker='.',markersize=15,mew=5,latlon=True)

isect7 = pl.where(NAs[0]==Ar); isect7=(Ar[isect7[0][0],0],Ar[isect7[0][0],1])
isect8 = pl.where(NAs[-1]==Ar); isect8=(Ar[isect8[0][0],0],Ar[isect8[0][0],1])
m.plot(isect7[0],isect7[1],'rx',markersize=15,mew=5,latlon=True)
m.plot(isect8[0],isect8[1],'rx',markersize=15,mew=5,latlon=True)

m.plot(Am[67,0]+360,Am[67,1],'bx',markersize=15,mew=5,latlon=True) # amr split
m.plot(AfMe[51,0],AfMe[51,1],'kx',markersize=15,mew=5,latlon=True) # afr split 1
m.plot(AfMe[182,0],AfMe[182,1],'kx',markersize=15,mew=5,latlon=True) # afr split 2
m.plot(EAA[51,0],EAA[51,1],'gx',markersize=15,mew=5,latlon=True) # eaa split
m.plot(Ar[121,0],Ar[121,1],'rx',markersize=15,mew=5,latlon=True) # ara aplit
m.plot(Ar[291,0],Ar[291,1],'rx',markersize=15,mew=5,latlon=True) # ari split 1
m.plot(Ar[325,0],Ar[325,1],'rx',markersize=15,mew=5,latlon=True) # ari split 2
m.plot(Ar[431,0],Ar[431,1],'rx',markersize=15,mew=5,latlon=True) # arp split 1
m.plot(Ar[517,0],Ar[517,1],'rx',markersize=15,mew=5,latlon=True) # arp split 2
m.plot(SO[489,0],SO[489,1],'x',mec='darkgoldenrod',markersize=15,mew=5,latlon=True) # soa split
m.plot(SO[489,0],SO[489,1],'x',mec='darkgoldenrod',markersize=15,mew=5,latlon=True) # soi split
m.plot(SO[112,0],SO[112,1],'x',mec='darkgoldenrod',markersize=15,mew=5,latlon=True) # soi split
m.plot(SO[213,0],SO[213,1],'x',mec='darkgoldenrod',markersize=15,mew=5,latlon=True) # sop split 1
m.plot(SO[363,0],SO[363,1],'x',mec='darkgoldenrod',markersize=15,mew=5,latlon=True) # sop split 2


fluxes = pl.array([0.23,0.17,0.16,0.98,0.17,0.01,0.05,0.07,0.29,0.27,0.38,0.09,0.12,0.07])
amr_flx = [0.17,0.40]; afr_flx = [0.10,0.34,0.08]; eaa_flx = [0.19,0.03]
ara_flx = [0.09,0.08,0.12]; ari_flx = [0.09,0.06,0.02]; arp_flx = [0.07,0.04,0.11]
soa_flx = [0.02,0.27]; soi_flx = [0.20,0.07]; sop_flx = [0.03,0.37,0.05]

lon = pl.array([-85,55,80,-97,66,188,-22,80,-130,77,42,101])
lat = pl.array([18.4,2.5,6.1,60,23,60,-30,-30,-30,45,48,38])
x,y = m(lon, lat)

amr_lon = [245,290]; amr_lat = [34,0]; x1,y1 = m(amr_lon, amr_lat)
afr_lon = [26,55,16]; afr_lat = [27,2.5,-29]; x2,y2 = m(afr_lon, afr_lat)
eaa_lon = [85,123]; eaa_lat = [16,-2]; x3,y3 = m(eaa_lon, eaa_lat)
ara_lon = [266,9,36]; ara_lat = [70,72,49]; x4,y4 = m(ara_lon, ara_lat)
ari_lon = [44,65,88]; ari_lat = [33,34,26]; x5,y5 = m(ari_lon, ari_lat)
arp_lon = [100,130,219]; arp_lat = [40,60,57]; x6,y6 = m(arp_lon,arp_lat)
soa_lon = [308,327]; soa_lat = [-12,-30]; x7,y7 = m(soa_lon,soa_lat)
soi_lon = [76,122]; soi_lat = [-30,-21]; x8,y8 = m(soi_lon,soi_lat)
sop_lon = [145,220,286]; sop_lat = [-32,-32,-26]; x9,y9 = m(sop_lon,sop_lat)
#
pl.arrow(x1[0],y1[0],1500000,0,fc="none",ec="b",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(amr_flx[0]),xy=(0.62,0.69),xycoords='figure fraction',
            color='b',size=18)
pl.arrow(x1[1],y1[1],-1500000,0,fc="none",ec="b",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate('{:.2f}'.format(amr_flx[1]),xy=(0.685,0.485),
            xycoords='figure fraction',color='b',size=18)
            
pl.arrow(x2[0],y2[0],990000,-250000,fc="none",ec="k",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate('{:.2f}'.format(afr_flx[0]),xy=(0.07,0.66),
            xycoords='figure fraction',color='k',size=18)
pl.arrow(x2[1],y2[1],-2500000,0,fc="none", ec="k",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(afr_flx[1]),xy=(0.15,0.48),xycoords='figure fraction',
            color='k',size=18)
pl.arrow(x2[2],y2[2],2000000,0,fc="none",ec="k",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(afr_flx[2]),xy=(0.05,0.32),xycoords='figure fraction',
            color='k',size=18)

pl.arrow(x3[0],y3[0],1700000,0,fc="none",ec="g",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(eaa_flx[0]),xy=(0.25,0.56),xycoords='figure fraction',
            color='g',size=18)
pl.arrow(x3[1],y3[1],-1000000,-1000000,fc="none",ec="g",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(eaa_flx[1]),xy=(0.27,0.40),xycoords='figure fraction',
            color='g',size=18)
          
pl.arrow(x4[0],y4[0],450000,-850000,fc="none",ec="r",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(ara_flx[0]),xy=(0.69,0.81),xycoords='figure fraction',
            color='r',size=18)
pl.arrow(x4[1],y4[1],1600000,0,fc="none",ec="r",lw=3,head_width=550000,
         head_length=650000,width=10000,head_starts_at_zero=False)
pl.annotate(str(ara_flx[1]),xy=(0.255,0.90),xycoords='figure fraction',
            color='r',size=18)
pl.arrow(x4[2],y4[2],1100000,-150000,fc="none",ec="r",lw=3,head_width=550000,
         head_length=650000,width=10000,head_starts_at_zero=False)
pl.annotate('{:.2f}'.format(ara_flx[2]),xy=(0.135,0.80),
            xycoords='figure fraction',color='r',size=18)

pl.arrow(x5[0],y5[0],1100000,250000,fc="none",ec="r",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(ari_flx[0]),xy=(0.14,0.666),xycoords='figure fraction',
            color='r',size=18,zorder=10)
pl.arrow(x5[1],y5[1],500000,-700000,fc="none",ec="r",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(ari_flx[1]),xy=(0.21,0.6),xycoords='figure fraction',color='r',
            size=18,zorder=10)
pl.arrow(x5[2],y5[2],0,700000,fc="none",ec="r",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(ari_flx[2]),xy=(0.25,0.74),xycoords='figure fraction',
            color='r',size=18,zorder=10)

pl.arrow(x6[0],y6[0],900000,-450000,fc="none",ec="r",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(arp_flx[0]),xy=(0.32,0.66),xycoords='figure fraction',
            color='r',size=18,zorder=10)
pl.arrow(x6[1],y6[1],900000,-510000,fc="none",ec="r",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(arp_flx[1]),xy=(0.43,0.775),xycoords='figure fraction',
            color='r',size=18,zorder=10)
pl.arrow(x6[2],y6[2],900000,510000,fc="none",ec="r",lw=3,head_width=500000,
         head_length=700000,width=9000,head_starts_at_zero=False)
pl.annotate(str(arp_flx[2]),xy=(0.55,0.805),xycoords='figure fraction',
            color='r',size=18,zorder=10)
  
pl.arrow(x7[0],y7[0],-450000,-750000,fc="none",ec="darkgoldenrod",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(soa_flx[0]),xy=(0.81,0.44),xycoords='figure fraction',
            color='darkgoldenrod',size=20)
pl.arrow(x7[1],y7[1],0,-750000,fc="none",ec="darkgoldenrod",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(soa_flx[1]),xy=(0.85,0.33),xycoords='figure fraction',
            color='darkgoldenrod',size=20)

pl.arrow(x8[0],y8[0],0,-750000,fc="none",ec="darkgoldenrod",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate('{:.2f}'.format(soi_flx[0]),xy=(0.215,0.33),
            xycoords='figure fraction',color='darkgoldenrod',size=20)
pl.arrow(x8[1],y8[1],450000,-750000,fc="none",ec="darkgoldenrod",lw=3,
         head_width=600000,head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(soi_flx[1]),xy=(0.355,0.265),xycoords='figure fraction',
            color='darkgoldenrod',size=20)

pl.arrow(x9[0],y9[0],650000,650000,fc="none",ec="darkgoldenrod",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(sop_flx[0]),xy=(0.43,0.38),xycoords='figure fraction',
            color='darkgoldenrod',size=20)
pl.arrow(x9[1],y9[1],0,-750000,fc="none",ec="darkgoldenrod",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(sop_flx[1]),xy=(0.585,0.32),xycoords='figure fraction',
            color='darkgoldenrod',size=20)
pl.arrow(x9[2],y9[2],650000,0,fc="none",ec="darkgoldenrod",lw=3,head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False)
pl.annotate(str(sop_flx[2]),xy=(0.725,0.335),xycoords='figure fraction',
            color='darkgoldenrod',size=20)

#
#pl.arrow(x[6],y[6],0,-850000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=900000,
#         head_length=900000,width=10000,head_starts_at_zero=False)
#pl.annotate(str(fluxes[8]),xy=(0.83,0.25),xycoords='figure fraction',color='darkgoldenrod',
#            size=30)
#
#pl.arrow(x[7],y[7],0,-850000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=900000,
#         head_length=900000,width=10000,head_starts_at_zero=False,zorder=10)
#pl.annotate(str(fluxes[9]),xy=(0.19,0.25),xycoords='figure fraction',color='darkgoldenrod',
#            size=30)
#
#pl.arrow(x[8],y[8],0,-850000,fc="none", ec="darkgoldenrod", linewidth = 4, head_width=900000,
#         head_length=900000,width=10000,head_starts_at_zero=False,zorder=10)
#pl.annotate(str(fluxes[10]),xy=(0.56,0.25),xycoords='figure fraction',color='darkgoldenrod',
#            size=30)
#
pl.arrow(x[-3],y[-3],350000,650000,fc="none", ec="maroon", linewidth = 3, head_width=500000,
         head_length=600000,width=9000,head_starts_at_zero=False,zorder=10)
pl.annotate('{:.2f}'.format(fluxes[-3]),xy=(0.285,0.82),xycoords='figure fraction',color='maroon',
            size=20)
#            
#pl.arrow(x[-2],y[-2],600000,-150000,fc="none", ec="maroon", linewidth = 4, head_width=800000,
#         head_length=900000,width=10000,head_starts_at_zero=False)
#pl.annotate(str(fluxes[-2]),xy=(0.12,0.75),xycoords='figure fraction',color='maroon',
#            size=30)
#
#pl.arrow(x[-1],y[-1],600000,-150000,fc="none", ec="maroon", linewidth = 4, head_width=800000,
#         head_length=900000,width=10000,head_starts_at_zero=False)
#pl.annotate(str(fluxes[-1]),xy=(0.35,0.65),xycoords='figure fraction',color='maroon',
#            size=30)

a,b = m(180,15)
pl.text(a,b,'0.00',bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
a,b = m(320,15)
pl.text(a,b,'-0.47',bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
a,b = m(61,-12)
pl.text(a,b,'-0.64',bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)
a,b = m(100,70)
pl.text(a,b,'0.16',bbox={'facecolor':'white', 'alpha':1.0, 'pad':5},fontsize=23)
a,b = m(180,-55)
pl.text(a,b,'0.95',bbox={'facecolor':'white', 'alpha':0.5, 'pad':5},fontsize=23)


a,b = m(33,41)
pl.text(a,b,'A',bbox={'facecolor':'white','boxstyle':'circle','alpha':0.5},fontsize=22)
a,b = m(15,-44)
pl.text(a,b,'B',bbox={'facecolor':'white','boxstyle':'circle','alpha':0.5},fontsize=22)
a,b = m(98,29)
pl.text(a,b,'C',bbox={'facecolor':'white','boxstyle':'circle','alpha':0.5},fontsize=22)
a,b = m(133,-18)
pl.text(a,b,'D',bbox={'facecolor':'white','boxstyle':'circle','alpha':0.5},fontsize=22)
a,b = m(226,47)
pl.text(a,b,'E',bbox={'facecolor':'white','boxstyle':'circle','alpha':0.5},fontsize=22)
a,b = m(277,-18)
pl.text(a,b,'F',bbox={'facecolor':'white','boxstyle':'circle','alpha':0.5},fontsize=22)


#pl.annotate('AM1',xy=(0.62,0.685),xycoords='figure fraction',color='b',
#            size=50)
#pl.annotate('AM2',xy=(0.71,0.485),xycoords='figure fraction',color='b',
#            size=50)
#pl.annotate('AF1',xy=(0.055,0.66),xycoords='figure fraction',color='k',
#            size=50)
#pl.annotate('AF2',xy=(0.11,0.48),xycoords='figure fraction',color='k',
#            size=50)
#pl.annotate('AF3',xy=(0.11,0.31),xycoords='figure fraction',color='k',
#            size=50)
#pl.annotate('EA1',xy=(0.29,0.58),xycoords='figure fraction',color='g',
#            size=50)
#pl.annotate('EA2',xy=(0.33,0.46),xycoords='figure fraction',color='g',
#            size=50)
#pl.annotate('AA1',xy=(0.69,0.84),xycoords='figure fraction',color='r',
#            size=50,zorder=5)
#pl.annotate('AA2',xy=(0.23,0.90),xycoords='figure fraction',color='r',
#            size=50)
#pl.annotate('AA3',xy=(0.10,0.77),xycoords='figure fraction',color='r',
#            size=50)
#pl.annotate('AI1',xy=(0.13,0.666),xycoords='figure fraction',color='r',
#            size=35,zorder=10)
#pl.annotate('AI2',xy=(0.20,0.64),xycoords='figure fraction',color='r',
#            size=35,zorder=10)
#pl.annotate('AI3',xy=(0.2375,0.71),xycoords='figure fraction',color='r',
#            size=35,zorder=10)
#pl.annotate('AP1',xy=(0.31,0.69),xycoords='figure fraction',color='r',
#            size=50,zorder=10)
#pl.annotate('AP2',xy=(0.42,0.80),xycoords='figure fraction',color='r',
#            size=50,zorder=10)
#pl.annotate('AP3',xy=(0.545,0.82),xycoords='figure fraction',color='r',
#            size=50,zorder=10)
#pl.annotate('SA1',xy=(0.81,0.42),xycoords='figure fraction',color='darkgoldenrod',
#            size=50)
#pl.annotate('SA2',xy=(0.89,0.31),xycoords='figure fraction',color='darkgoldenrod',
#            size=50)
#pl.annotate('SI1',xy=(0.23,0.31),xycoords='figure fraction',color='darkgoldenrod',
#            size=50)
#pl.annotate('SI2',xy=(0.325,0.36),xycoords='figure fraction',color='darkgoldenrod',
#            size=50)
#pl.annotate('SP1',xy=(0.425,0.345),xycoords='figure fraction',color='darkgoldenrod',
#            size=50)
#pl.annotate('SP2',xy=(0.58,0.31),xycoords='figure fraction',color='darkgoldenrod',
#            size=50)
#pl.annotate('SP3',xy=(0.74,0.34),xycoords='figure fraction',color='darkgoldenrod',
#            size=50)

pl.tight_layout()
pl.savefig('/home/users/qx911590/Figure1.pdf',dpi=350)