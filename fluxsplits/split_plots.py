# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 15:22:57 2017

@author: np838619
"""

import pylab as pl

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

pl.close('all')
sheddir = '/home/np838619/Watershed/fluxsplits/'

#amr = pl.genfromtxt(sheddir+'shed_defs/'+'Am_traj_release_new.txt',skip_header=5)
#afr = pl.genfromtxt(sheddir+'shed_defs/'+'AfMe_traj_release_new.txt',skip_header=5)
#eaa = pl.genfromtxt(sheddir+'shed_defs/'+'EAA_traj_release_new.txt',skip_header=5)
#ara = pl.genfromtxt(sheddir+'shed_defs/'+'ArA_traj_release_new.txt',skip_header=5)
#ari = pl.genfromtxt(sheddir+'shed_defs/'+'ArI_traj_release_new.txt',skip_header=5)
#arp = pl.genfromtxt(sheddir+'shed_defs/'+'ArP_traj_release_new.txt',skip_header=5)
#soa = pl.genfromtxt(sheddir+'shed_defs/'+'SOA_traj_release_new.txt',skip_header=5)
#soi = pl.genfromtxt(sheddir+'shed_defs/'+'SOI_traj_release_new.txt',skip_header=5)
#sop = pl.genfromtxt(sheddir+'shed_defs/'+'SOP_traj_release_new.txt',skip_header=5)

amr_sp = pl.genfromtxt(sheddir+'splits_Am_interann.csv',delimiter=',')
afr_sp = pl.genfromtxt(sheddir+'splits_AfMe_interann.csv',delimiter=',')
eaa_sp = pl.genfromtxt(sheddir+'splits_EAA_interann.csv',delimiter=',')
ara_sp = pl.genfromtxt(sheddir+'splits_ArA_interann.csv',delimiter=',')
ari_sp = pl.genfromtxt(sheddir+'splits_ArI_interann.csv',delimiter=',')
arp_sp = pl.genfromtxt(sheddir+'splits_ArP_interann.csv',delimiter=',')
soa_sp = pl.genfromtxt(sheddir+'splits_SOA_interann.csv',delimiter=',')
soi_sp = pl.genfromtxt(sheddir+'splits_SOI_interann.csv',delimiter=',')
sop_sp = pl.genfromtxt(sheddir+'splits_SOP_interann.csv',delimiter=',')

years = pl.linspace(1979,2014,36)

###############################################################################
ax,fig = pl.subplots(3,3,figsize=(12,10))

ax1 = pl.subplot(331)
ax1.plot(years,amr_sp[:,0],color='b',ls='-',lw=2,label='AM1')
ax1.plot(years,amr_sp[:,1],color='b',ls='--',lw=2,label='AM2')
pl.xlim(1979,2014,35); ax1.set_xticks(pl.linspace(1979,2014,8))
pl.ylim(-0.5,0.3); pl.ylabel('Sv',fontsize=22)
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.yticks(fontsize=14)#; pl.xticks([-0.5,' ',-0.3,])
ax1.grid()#ax1.tick_params(axis='x', pad=10); 
ax1.legend(loc=6,fontsize=18,ncol=1,columnspacing=0.5)
pl.title('(a) Americas',loc='left',fontsize=18)

ax2 = pl.subplot(332)
ax2.plot(years,afr_sp[:,0],color='k',ls='-',lw=2,label='AF1')
ax2.plot(years,afr_sp[:,1],color='k',ls='--',lw=2,label='AF2')
ax2.plot(years,afr_sp[:,2],color='k',ls=':',ms=5,lw=2,label='AF3')
pl.xlim(1979,2014,35); ax2.set_xticks(pl.linspace(1979,2014,8))
pl.ylim(-0.5,0.3); ax2.yaxis.set_major_formatter(pl.NullFormatter())
ax2.xaxis.set_major_formatter(pl.NullFormatter())
#pl.xticks(fontsize=15,rotation=45)
ax2.grid()#; ax2.tick_params(axis='x', pad=10)
ax2.legend(loc=6,fontsize=16,ncol=2,columnspacing=0.5)
pl.title('(b) Africa',loc='left',fontsize=18)

ax3 = pl.subplot(333)
ax3.plot(years,eaa_sp[:,0],color='g',ls='-',lw=2.,label='EA1')
ax3.plot(years,eaa_sp[:,1],color='g',ls='--',lw=2.,label='EA2')
ax3.xaxis.set_major_formatter(pl.NullFormatter())
pl.xlim(1979,2014,35); ax3.set_xticks(pl.linspace(1979,2014,8))
pl.ylim(-0.5,0.3); ax3.yaxis.set_major_formatter(pl.NullFormatter())
#pl.xticks(fontsize=15,rotation=45)
ax3.grid()#; ax3.tick_params(axis='x', pad=10)
ax3.legend(loc=0,fontsize=18)
pl.title('(c) South-East Asia',loc='left',fontsize=18)

#pl.tight_layout()
#pl.savefig(sheddir+'amr_afr_eaa_splits.png')
###############################################################################

###############################################################################
#ax,fig = pl.subplots(1,3,figsize=(17,5))

ax4 = pl.subplot(334)
ax4.plot(years,ara_sp[:,0],color='r',ls='-',lw=2,label='AA1')
ax4.plot(years,ara_sp[:,1],color='r',ls='--',lw=2,label='AA2')
ax4.plot(years,ara_sp[:,2],color='r',ls=':',lw=2,label='AA3')
pl.xlim(1979,2014,35); ax4.set_xticks(pl.linspace(1979,2014,8))
pl.ylim(-0.2,0.2); pl.ylabel('Sv',fontsize=22)
pl.yticks(fontsize=14)#; pl.xticks(fontsize=15,rotation=45)
ax4.xaxis.set_major_formatter(pl.NullFormatter())
ax4.grid()#ax1.tick_params(axis='x', pad=10); 
ax4.legend(loc=6,ncol=2,fontsize=15,columnspacing=0.5)
pl.title('(d) Atlantic sector',loc='left',fontsize=18)

ax5 = pl.subplot(335)
ax5.plot(years,ari_sp[:,0],color='r',ls='-',lw=2,label='AI1')
ax5.plot(years,ari_sp[:,1],color='r',ls='--',lw=2,label='AI2')
ax5.plot(years,ari_sp[:,2],color='r',ls=':',lw=2,label='AI3')
pl.xlim(1979,2014,35); ax5.set_xticks(pl.linspace(1979,2014,8))
pl.ylim(-0.2,0.2)
pl.yticks(fontsize=16)#; pl.xticks(fontsize=15,rotation=45)
ax5.xaxis.set_major_formatter(pl.NullFormatter())
ax5.grid()#ax2.tick_params(axis='x', pad=10); 
ax5.yaxis.set_major_formatter(pl.NullFormatter())
ax5.legend(loc=3,ncol=2,fontsize=15,columnspacing=0.5)
pl.title('(e) Indian sector',loc='left',fontsize=18)

ax6 = pl.subplot(336)
ax6.plot(years,arp_sp[:,0],color='r',ls='-',lw=2,label='AP1')
ax6.plot(years,arp_sp[:,1],color='r',ls='--',lw=2,label='AP2')
ax6.plot(years,arp_sp[:,2],color='r',ls=':',lw=2,label='AP3')
pl.xlim(1979,2014,35); ax6.set_xticks(pl.linspace(1979,2014,8))
pl.ylim(-0.2,0.2)
ax6.xaxis.set_major_formatter(pl.NullFormatter())
pl.yticks(fontsize=16)#; pl.xticks(fontsize=15,rotation=45)
ax6.grid()#ax3.tick_params(axis='x', pad=10); 
ax6.yaxis.set_major_formatter(pl.NullFormatter())
ax6.legend(loc=3,columnspacing=0.5,ncol=2,fontsize=14)
pl.title('(f) Pacific sector',loc='left',fontsize=18)

#pl.tight_layout()
#pl.savefig(sheddir+'arctic_splits.png')
###############################################################################

###############################################################################
#ax,fig = pl.subplots(1,3,figsize=(17,5))

ax7 = pl.subplot(337)
ax7.plot(years,soa_sp[:,0],color='darkgoldenrod',ls='-',lw=2,label='SA1')
ax7.plot(years,soa_sp[:,1],color='darkgoldenrod',ls='--',lw=2,label='SA2')
ax7.set_xticks(pl.linspace(1979,2014,8)); pl.xlim(1979,2014,35)
pl.ylim(-0.5,0.1); pl.ylabel('Sv',fontsize=22)
pl.yticks(fontsize=16); pl.xticks(fontsize=15,rotation=45)
ax7.grid()#ax1.tick_params(axis='x', pad=10);
ax7.legend(loc=3,fontsize=18,ncol=2,columnspacing=0.5)
pl.title('(g) Atlantic sector',loc='left',fontsize=18)

ax8 = pl.subplot(338)
ax8.plot(years,soi_sp[:,0],color='darkgoldenrod',ls='-',lw=2,label='SI1')
ax8.plot(years,soi_sp[:,1],color='darkgoldenrod',ls='--',lw=2,label='SI2')
ax8.set_xticks(pl.linspace(1979,2014,8)); pl.xlim(1979,2014,35)
pl.ylim(-0.5,0.1)
pl.yticks(fontsize=16); pl.xticks(fontsize=15,rotation=45)
ax8.grid()#ax2.tick_params(axis='x', pad=10); 
ax8.yaxis.set_major_formatter(pl.NullFormatter())
ax8.legend(loc=3,fontsize=18)
pl.title('(h) Indian sector',loc='left',fontsize=18)

ax9 = pl.subplot(339)
ax9.plot(years,sop_sp[:,0],color='darkgoldenrod',ls='-',lw=2,label='SP1')
ax9.plot(years,sop_sp[:,1],color='darkgoldenrod',ls='--',lw=2,label='SP2')
ax9.plot(years,sop_sp[:,2],color='darkgoldenrod',ls=':',lw=2,label='SP3')
ax9.set_xticks(pl.linspace(1979,2014,8)); pl.xlim(1979,2014,35)
pl.ylim(-0.5,0.1)
pl.yticks(fontsize=16); pl.xticks(fontsize=15,rotation=45)
ax9.grid()#ax2.tick_params(axis='x', pad=10); 
ax9.yaxis.set_major_formatter(pl.NullFormatter())
ax9.legend(loc=7,fontsize=18,ncol=2,columnspacing=0.5)
pl.title('(i) Pacific sector',loc='left',fontsize=18)

pl.tight_layout()
#pl.savefig(sheddir+'shed_splits_iv_all.png')
###############################################################################
#pl.close('all')

amr_ss = pl.genfromtxt(sheddir+'splits_Am_sscyc.csv',delimiter=',')
afr_ss = pl.genfromtxt(sheddir+'splits_AfMe_sscyc.csv',delimiter=',')
eaa_ss = pl.genfromtxt(sheddir+'splits_EAA_sscyc.csv',delimiter=',')
ara_ss = pl.genfromtxt(sheddir+'splits_ArA_sscyc.csv',delimiter=',')
ari_ss = pl.genfromtxt(sheddir+'splits_ArI_sscyc.csv',delimiter=',')
arp_ss = pl.genfromtxt(sheddir+'splits_ArP_sscyc.csv',delimiter=',')
soa_ss = pl.genfromtxt(sheddir+'splits_SOA_sscyc.csv',delimiter=',')
soi_ss = pl.genfromtxt(sheddir+'splits_SOI_sscyc.csv',delimiter=',')
sop_ss = pl.genfromtxt(sheddir+'splits_SOP_sscyc.csv',delimiter=',')

amr_ss = SeasonalCyc(amr_ss); afr_ss = SeasonalCyc(afr_ss)
eaa_ss = SeasonalCyc(eaa_ss); ara_ss = SeasonalCyc(ara_ss)
ari_ss = SeasonalCyc(ari_ss); arp_ss = SeasonalCyc(arp_ss)
soa_ss = SeasonalCyc(soa_ss); soi_ss = SeasonalCyc(soi_ss)
sop_ss = SeasonalCyc(sop_ss)

ssns = ['DJF','MAM','JJA','SON']

###############################################################################
ax,fig = pl.subplots(3,3,figsize=(12,10.5))

ind = pl.arange(4); width = 0.4
ax1 = pl.subplot(331)
ax1.bar(ind,amr_ss[:,0],width,color='b',ls='-',lw=2,label='AM1')
ax1.bar(ind+width,amr_ss[:,1],width,color='b',ec='k',hatch='/',ls='-',lw=2,label='AM2')
pl.ylim(-0.6,0.6); ax1.grid(axis='y')#pl.xlim(0,11)
ax1.set_yticks(pl.linspace(-0.6,0.6,7))
pl.xticks(ind+width); ax1.set_xticklabels(ssns); pl.xticks(fontsize=16)
pl.yticks(fontsize=15); pl.ylabel('Sv',fontsize=22,labelpad=-5)
ax1.legend(loc=0,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(a) Americas',fontsize=18,loc='left')
#ax1.annotate('Americas',(0.15,0.93),xycoords='figure fraction',fontsize=18)

width = 0.3
ax2 = pl.subplot(332)
ax2.bar(ind,afr_ss[:,0],width,ec='w',color='k',hatch='/',lw=0,label='AF1')
ax2.bar(ind+width,afr_ss[:,1],width,color='k',lw=2,label='AF2')
ax2.bar(ind+2*width,afr_ss[:,2],width,ec='w',color='k',hatch='|',lw=0,label='AF3')
pl.ylim(-0.6,0.6); ax2.grid(axis='y'); pl.xticks(fontsize=16)
pl.xticks(ind+1.5*width); ax2.set_xticklabels(ssns)
pl.yticks(fontsize=15)#; pl.ylabel('Sv',fontsize=22,labelpad=-5)
ax2.legend(loc=0,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(b) Africa',fontsize=18,loc='left')

width = 0.4
ax3 = pl.subplot(333)
ax3.bar(ind,eaa_ss[:,0],width,color='g',lw=2,label='EA1')
ax3.bar(ind+width,eaa_ss[:,1],width,color='g',ec='k',hatch='/',lw=2,label='EA2')
#ax3.plot(eaa_ss[:,2],color='g',ls=':',lw=2,label='7S-20S')
pl.ylim(-0.6,0.6); ax3.grid(axis='y')
ax3.set_yticks(pl.linspace(-0.6,0.6,7))
pl.yticks(fontsize=15); pl.xticks(fontsize=16)
pl.xticks(ind+width)#ax3.set_xticks(pl.linspace(0,11,12))
ax3.set_xticklabels(ssns)#; pl.ylabel('Sv',fontsize=22,labelpad=-5)
ax3.legend(loc=0,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(c) South-East Asia',fontsize=18,loc='left')

#pl.tight_layout()
#pl.savefig(sheddir+'amr_afr_eaa_sscyc_bars.png')
###############################################################################

###############################################################################
#ax,fig = pl.subplots(3,1,figsize=(8,10))

width=0.3
ax4 = pl.subplot(334)
ax4.bar(ind,ara_ss[:,0],width,color='r',hatch='/',ec='k',lw=2,label='AA1')
ax4.bar(ind+width,ara_ss[:,1],width,color='r',lw=2,label='AA2')
ax4.bar(ind+2*width,ara_ss[:,2],width,color='r',hatch='|',ec='k',lw=2,label='AA3')
pl.ylim(-0.3,0.2); ax4.grid(axis='y')
ax4.set_yticks(pl.linspace(-0.3,0.2,6)); pl.xticks(ind+width)
ax4.set_xticklabels(ssns); pl.xticks(fontsize=15)
pl.yticks(fontsize=15); pl.ylabel('Sv',fontsize=22,labelpad=-3)
ax4.legend(loc=3,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(d) Atlantic sector',fontsize=18,loc='left')
#ax1.annotate('Americas',(0.15,0.93),xycoords='figure fraction',fontsize=18)

ax5 = pl.subplot(335)
ax5.bar(ind,ari_ss[:,0],width,color='r',hatch='/',ec='k',lw=2,label='AI1')
ax5.bar(ind+width,ari_ss[:,1],width,color='r',lw=2,label='AI2')
ax5.bar(ind+2*width,ari_ss[:,2],width,color='r',hatch='|',ec='k',lw=2,label='AI3')
pl.ylim(-0.3,0.2); ax5.grid(axis='y')
ax5.set_yticks(pl.linspace(-0.3,0.2,6)); pl.xticks(ind+width)
ax5.set_xticklabels(ssns); pl.xticks(fontsize=15)
pl.yticks(fontsize=15)#; pl.ylabel('Sv',fontsize=22,labelpad=-3)
ax5.legend(loc=3,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(e) Indian sector',fontsize=18,loc='left')

ax6 = pl.subplot(336)
ax6.bar(ind,arp_ss[:,0],width,color='r',hatch='/',ec='k',lw=2,label='AP1')
ax6.bar(ind+width,arp_ss[:,1],width,color='r',lw=2,label='AP2')
ax6.bar(ind+2*width,arp_ss[:,2],width,color='r',hatch='|',ec='k',lw=2,label='AP3')
pl.ylim(-0.3,0.2); ax6.grid(axis='y')
ax6.set_yticks(pl.linspace(-0.3,0.2,6))
pl.yticks(fontsize=15); pl.xticks(ind+width); pl.xticks(fontsize=15)
ax6.set_xticklabels(ssns)#; pl.ylabel('Sv',fontsize=22,labelpad=-3)
ax6.legend(loc=3,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(f) Pacific sector',fontsize=18,loc='left')

#pl.tight_layout()
#pl.savefig(sheddir+'arctic_sscyc_bars.png')
###############################################################################

###############################################################################
#ax,fig = pl.subplots(3,1,figsize=(8,10))
#
dg = 'darkgoldenrod'; width = 0.4
ax7 = pl.subplot(337)
ax7.bar(ind,soa_ss[:,0],width,color=dg,lw=2,label='SA1')
ax7.bar(ind+width,soa_ss[:,1],width,color=dg,hatch='/',ec='k',lw=2,label='SA2')
pl.ylim(-0.6,0.6); ax7.grid(axis='y')
ax7.set_yticks(pl.linspace(-0.6,0.6,7)); pl.xticks(ind+width)
ax7.set_xticklabels(ssns); pl.xticks(fontsize=15)
pl.yticks(fontsize=15); pl.ylabel('Sv',fontsize=22,labelpad=-3)
ax7.legend(loc=2,fontsize=16,ncol=2,columnspacing=0.5)
pl.title('(g) Atlantic sector',fontsize=18,loc='left')

ax8 = pl.subplot(338)
ax8.bar(ind,soi_ss[:,0],width,color=dg,lw=2,label='SI1')
ax8.bar(ind+width,soi_ss[:,1],width,color=dg,lw=2,hatch='/',ec='k',label='SI2')
pl.ylim(-0.6,0.26); ax8.grid(axis='y'); pl.xticks(ind+width)
ax8.set_xticklabels(ssns); pl.xticks(fontsize=15)
ax8.set_yticks(pl.linspace(-0.6,0.6,7))
pl.yticks(fontsize=15)#; pl.ylabel('Sv',fontsize=22,labelpad=-3)
ax8.legend(loc=2,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(h) Indian sector',fontsize=18,loc='left')

width = 0.3
ax9 = pl.subplot(339)
ax9.bar(ind,sop_ss[:,0],width,color=dg,hatch='/',ec='k',lw=2,label='SP1')
ax9.bar(ind+width,sop_ss[:,1],width,color=dg,lw=2,label='SP2')
ax9.bar(ind+2*width,sop_ss[:,2],width,color=dg,hatch='|',ec='k',lw=2,label='SP3')
pl.ylim(-0.6,0.6); ax9.grid(axis='y')
ax9.set_yticks(pl.linspace(-0.6,0.6,7))
pl.yticks(fontsize=15); pl.xticks(fontsize=15)
ax9.set_xticks(ind)
ax9.set_xticklabels(ssns)#; pl.ylabel('Sv',fontsize=22,labelpad=-3)
ax9.legend(loc=2,ncol=2,fontsize=16,columnspacing=0.5)
pl.title('(i) Pacific sector',fontsize=18,loc='left')

pl.tight_layout()
#pl.savefig(sheddir+'shed_splits_sc_all.png')
###############################################################################