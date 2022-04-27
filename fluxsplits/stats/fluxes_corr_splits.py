# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:32:20 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from scipy.stats import pearsonr, linregress


def split_corrs(fullflux,splitflux,divQ,shed):
    """
    """
    if shed == 'amr':
        I = 0.; J = 2.
    elif shed == 'afr':
        I = 1.; J = 0.
    elif shed == 'eaa':
        I = 2.; J = 1.
    elif shed == 'ara':
        I = 3.; J = 0.
    elif shed == 'ari':
        I = 3.; J = 1.
    elif shed == 'arp':
        I = 3.; J = 2.
    elif shed == 'soa':
        I = 4.; J = 0.
    elif shed == 'soi':
        I = 4.; J = 1.
    elif shed == 'sop':
        I = 4.; J = 2.
    
    flx_corrs = pl.zeros([splitflux.shape[1],2])
    bas_cors1 = pl.zeros_like(flx_corrs)
    bas_cors2 = pl.zeros_like(flx_corrs)
    
    for i in range(splitflux.shape[1]):
        fr = pearsonr(fullflux,splitflux[:,i])
        flx_corrs[i,0] = fr[0]; flx_corrs[i,1] = fr[1]
        br1 = pearsonr(divQ[:,I],splitflux[:,i])
        br2 = pearsonr(divQ[:,J],splitflux[:,i])
        bas_cors1[i,0] = br1[0]; bas_cors1[i,1] = br1[1]
        bas_cors2[i,0] = br2[0]; bas_cors2[i,1] = br2[1]
    
    all_corrs = pl.zeros([splitflux.shape[1],flx_corrs.shape[1]*3])
    all_corrs[:,:2] = flx_corrs
    all_corrs[:,2:4] = bas_cors1
    all_corrs[:,4:] = bas_cors2
    
#    sheddir = '/home/np838619/Watershed/' + 'fluxsplits/stats/'
#    f = open(sheddir+shed+'_stats.csv','w')
#    f.write('flux_flux, , dq1_flux, , dq2_flux, \n')
#    f.write('r, sig, r, sig, r, sig\n')
#    pl.savetxt(f,all_corrs,fmt='%9.5f',delimiter=',')
#    f.close()
    
    return all_corrs#flx_corrs, bas_cors1, bas_cors2

def Covs(fluxes,divq,shed):
    if shed == 'amr':
        I = 0.
    elif shed == 'afr':
        I = 1.
    elif shed == 'eaa':
        I = 2.
    elif shed == 'ara':
        I = 3.
    elif shed == 'ari':
        I = 4
    elif shed == 'arp':
        I = 5
    elif shed == 'soa':
        I = 6
    elif shed == 'soi':
        I = 7
    elif shed == 'sop':
        I = 8
    
    atlcov = pl.cov(divq[:,0],fluxes[:,I])
    indcov = pl.cov(divq[:,1],fluxes[:,I])
    paccov = pl.cov(divq[:,2],fluxes[:,I])
    arccov = pl.cov(divq[:,3],fluxes[:,I])
    soucov = pl.cov(divq[:,4],fluxes[:,I])
    
    covs = pl.array([atlcov[0,1],indcov[0,1],paccov[0,1],arccov[0,1],soucov[0,1]])
    
    return covs

def ScatterPlot(fullflux,splitflux,colour):
    """
    """
    r = linregress(fullflux,splitflux)
    x = pl.linspace(-2,2,11); y = r[0]*x + r[1]
    pl.scatter(fullflux,splitflux,color=colour,edgecolors='k',linewidths=0.6,zorder=2)
    pl.plot(x,y,label='regression line',color=colour,linewidth=2,zorder=0)
    #pl.plot(x,x,ls='--',color='k',label='1:1 line')
    pl.xlim(-0.1,0.1); pl.ylim(-0.1,0.1)
    pl.axvline(x=0.1,color='k',ls=':'); pl.axvline(x=-0.1,color='k',ls=':')
    pl.axvline(x=0.2,color='k',ls=':'); pl.axvline(x=-0.2,color='k',ls=':')
    pl.axhline(y=0.1,color='k',ls=':'); pl.axhline(y=-0.1,color='k',ls=':')
    pl.axhline(y=0.2,color='k',ls=':'); pl.axhline(y=-0.2,color='k',ls=':')
#    pl.xlabel()
    pl.text(-0.09,0.08,'$r=$'+str(round(r[2],2)),fontsize=15)
    
    return r[0], r[3]

def MoveAxes(ax):
    """
    """
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    #ax.spines['bottom'].set_label('hello')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xticks([-0.2,-0.1,0,0.1,0.2])
    ax.set_yticks([-0.2,-0.1,0,0.1,0.2])

pl.close('all')
sheddir = '/home/np838619/Watershed/'

shedfluxes = pl.genfromtxt(sheddir+'variability/shedfluxes_7914_var.csv',
                                                          skip_header=1)
amr_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_Am_interann.csv',delimiter=',')
afr_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_AfMe_interann.csv',delimiter=',')
eaa_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_EAA_interann.csv',delimiter=',')
ara_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_ArA_interann.csv',delimiter=',')
ari_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_ArI_interann.csv',delimiter=',')
arp_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_ArP_interann.csv',delimiter=',')
soa_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_SOA_interann.csv',delimiter=',')
soi_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_SOI_interann.csv',delimiter=',')
sop_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_SOP_interann.csv',delimiter=',')
divQ = pl.genfromtxt(sheddir+'variability/divQ_7914_var.csv',skip_header=1)

amr_cors = split_corrs(shedfluxes[:,0],amr_splt,divQ,'amr')
afr_cors = split_corrs(shedfluxes[:,1],afr_splt,divQ,'afr')
eaa_cors = split_corrs(shedfluxes[:,2],eaa_splt,divQ,'eaa')
ara_cors = split_corrs(shedfluxes[:,3],ara_splt,divQ,'ara')
ari_cors = split_corrs(shedfluxes[:,4],ari_splt,divQ,'ari')
arp_cors = split_corrs(shedfluxes[:,5],arp_splt,divQ,'arp')
soa_cors = split_corrs(shedfluxes[:,6],soa_splt,divQ,'soa')
soi_cors = split_corrs(shedfluxes[:,7],soi_splt,divQ,'soi')
sop_cors = split_corrs(shedfluxes[:,8],sop_splt,divQ,'sop')

#cormat = pl.array([amrcor,afrcor,eaacor,aracor,aricor,arpcor,soacor,soicor,sopcor])
#sigmat = pl.array([amrsig,afrsig,eaasig,arasig,arisig,arpsig,soasig,soisig,sopsig])

#f = open(sheddir+'variability/flux_divq_corrs_sc.csv','w')
#f.write(' atl ind pac arc sou\n')
#pl.savetxt(f,cormat,fmt='%9.5f',delimiter=' ')
#f.close()

#f = open(sheddir+'variability/flux_divq_sigs_sc.csv','w')
#f.write(' atl ind pac arc sou\n')
#pl.savetxt(f,sigmat,fmt='%9.5f',delimiter=' ')
#f.close()

#amrcov = Covs(shedfluxes,divQ,'amr')
#afrcov = Covs(shedfluxes,divQ,'afr')
#eaacov = Covs(shedfluxes,divQ,'eaa')
#aracov = Covs(shedfluxes,divQ,'ara')
#aricov = Covs(shedfluxes,divQ,'ari')
#arpcov = Covs(shedfluxes,divQ,'arp')
#soacov = Covs(shedfluxes,divQ,'soa')
#soicov = Covs(shedfluxes,divQ,'soi')
#sopcov = Covs(shedfluxes,divQ,'sop')

#covmat = pl.array([amrcov,afrcov,eaacov,aracov,aricov,arpcov,soacov,soicov,sopcov])

#f = open(sheddir+'variability/flux_divq_covs_iv.csv','w')
#f.write(' atl ind pac arc sou\n')
#pl.savetxt(f,covmat,fmt='%9.5f',delimiter=' ')
#f.close()

flxstd = pl.std(shedfluxes,axis=0)
dvqstd = pl.std(divQ,axis=0)

slopes = pl.zeros([9]); rsigs = pl.zeros([9])

basin=0.
divQ = divQ - pl.mean(divQ,axis=0)
shedfluxes = shedfluxes - pl.mean(shedfluxes,axis=0)
afr_splt = afr_splt - pl.mean(afr_splt,axis=0)
qn = '$\mathbf{Q}\cdot\mathbf{\hat{n}}$'

fig, ax = pl.subplots(3,3,figsize=(10,10)) # 8,8 for 2x2!
###############################################################################
ax1 = pl.subplot(331)
MoveAxes(ax1); pl.title('(a) AF1 vs Africa '+qn, fontsize=15,loc='left',y=1.01)
slopes[0], rsigs[0] = ScatterPlot(shedfluxes[:,1],afr_splt[:,0],'gray')
###############################################################################
ax2 = pl.subplot(332)
MoveAxes(ax2); pl.title('(b) AF2 vs Africa '+qn, fontsize=15,loc='left',y=1.01)
slopes[1], rsigs[1] = ScatterPlot(shedfluxes[:,1],afr_splt[:,1],'gray')
###############################################################################
ax3 = pl.subplot(333)
MoveAxes(ax3); pl.title('(c) AF3 vs Africa '+qn, fontsize=15,loc='left',y=1.01)
slopes[2], rsigs[2] = ScatterPlot(shedfluxes[:,1],afr_splt[:,2],'gray')
###############################################################################
ax4 = pl.subplot(334)
MoveAxes(ax4); pl.title('(d) AF1 vs Indian $P-E$', fontsize=15,loc='left',y=1.01)
slopes[3], rsigs[3] = ScatterPlot(divQ[:,1],afr_splt[:,0],'b')
#################################################################################
ax5 = pl.subplot(335)
MoveAxes(ax5); pl.title('(e) AF2 vs Indian $P-E$', fontsize=15,loc='left',y=1.01)
slopes[4], rsigs[4] = ScatterPlot(divQ[:,1],afr_splt[:,1],'b')
################################################################################
ax6 = pl.subplot(336)
MoveAxes(ax6); pl.title('(f) AF3 vs Indian $P-E$', fontsize=15,loc='left',y=1.01)
slopes[5], rsigs[5] = ScatterPlot(divQ[:,1],afr_splt[:,2],'b')
###############################################################################
ax7 = pl.subplot(337)
MoveAxes(ax7); pl.title('(g) AF1 vs Atlantic $P-E$', fontsize=15,loc='left',y=1.01)
slopes[6], rsigs[6] = ScatterPlot(divQ[:,0],-afr_splt[:,0],'r')
##############################################################################
ax8 = pl.subplot(338)
MoveAxes(ax8); pl.title('(h) AF2 vs Atlantic $P-E$', fontsize=15,loc='left',y=1.01)
slopes[7], rsigs[7] = ScatterPlot(divQ[:,0],-afr_splt[:,1],'r')
##############################################################################
ax9 = pl.subplot(339)
MoveAxes(ax9); pl.title('(i) AF3 vs Atlantic $P-E$', fontsize=15,loc='left',y=1.01)
slopes[8], rsigs[8] = ScatterPlot(divQ[:,0],-afr_splt[:,2],'r')
##############################################################################

pl.tight_layout()#pl.subplots_adjust(left=0.04,right=0.96,top=0.90,bottom=0.06)
#pl.suptitle('Southern $-$div$\mathbf{Q}$ vs. $\overline{\mathbf{Q}}\cdot\hat{\mathbf{n}}$',fontsize=20)
#pl.savefig(sheddir+'fluxsplits/stats/afr_splits_regs.png')

#f = open(sheddir+'fluxsplits/stats/regress_afr_splits.csv','w')
#pl.savetxt(f,zip(slopes,rsigs),fmt='%9.5f',delimiter=' ')
#f.close()