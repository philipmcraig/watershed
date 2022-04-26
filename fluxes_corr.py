# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:32:20 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from scipy.stats import pearsonr, linregress
import matplotlib.gridspec as gridspec

def Corrs(fluxes,divq,shed):
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
    
    atlcor = pearsonr(divq[:,0],fluxes[:,I])
    indcor = pearsonr(divq[:,1],fluxes[:,I])
    paccor = pearsonr(divq[:,2],fluxes[:,I])
    arccor = pearsonr(divq[:,3],fluxes[:,I])
    soucor = pearsonr(divq[:,4],fluxes[:,I])
    
    cors = pl.array([atlcor[0],indcor[0],paccor[0],arccor[0],soucor[0]])
    sigs = pl.array([atlcor[1],indcor[1],paccor[1],arccor[1],soucor[1]])
    
    return cors, sigs

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

def ScatterPlot(divQ,flux):
    """
    """
    r = linregress(divQ,flux)
    x = pl.linspace(-2,2,11); y = r[0]*x + r[1]
    pl.scatter(divQ,flux,color='deeppink',edgecolors='k',linewidths=0.6,zorder=2)
    pl.plot(x,y,label='regression line',color='deeppink',linewidth=2,zorder=0)
    #pl.plot(x,x,ls='--',color='k',label='1:1 line')
    pl.xlim(-0.05,0.05); pl.ylim(-0.05,0.05)
    pl.axvline(x=0.025,color='k',ls=':'); pl.axvline(x=-0.025,color='k',ls=':')
    pl.axvline(x=0.05,color='k',ls=':'); pl.axvline(x=-0.05,color='k',ls=':')
    pl.axhline(y=0.025,color='k',ls=':'); pl.axhline(y=-0.025,color='k',ls=':')
    pl.axhline(y=0.05,color='k',ls=':'); pl.axhline(y=-0.05,color='k',ls=':')
#    pl.xlabel()
    pl.text(-0.047,0.036,'$r=$'+str(round(r[2],2)),fontsize=15)
    
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
    ax.set_xticks([-0.05,-0.025,0,0.025,0.05])
    ax.set_yticks([-0.05,-0.025,0,0.025,0.05])

pl.close('all')
sheddir = '/home/np838619/Watershed/'

shedfluxes = pl.genfromtxt(sheddir+'variability/shedfluxes_7914_var.csv',
                                                          skip_header=1)
sf_splt = pl.genfromtxt(sheddir+'fluxsplits/splits_Am_interann.csv',delimiter=',')
divQ = pl.genfromtxt(sheddir+'variability/divQ_7914_var.csv',skip_header=1)

amrcor, amrsig = Corrs(shedfluxes,divQ,'amr')
afrcor, afrsig = Corrs(shedfluxes,divQ,'afr')
eaacor, eaasig = Corrs(shedfluxes,divQ,'eaa')
aracor, arasig = Corrs(shedfluxes,divQ,'ara')
aricor, arisig = Corrs(shedfluxes,divQ,'ari')
arpcor, arpsig = Corrs(shedfluxes,divQ,'arp')
soacor, soasig = Corrs(shedfluxes,divQ,'soa')
soicor, soisig = Corrs(shedfluxes,divQ,'soi')
sopcor, sopsig = Corrs(shedfluxes,divQ,'sop')

cormat = pl.array([amrcor,afrcor,eaacor,aracor,aricor,arpcor,soacor,soicor,sopcor])
sigmat = pl.array([amrsig,afrsig,eaasig,arasig,arisig,arpsig,soasig,soisig,sopsig])

#f = open(sheddir+'variability/flux_divq_corrs_iv.csv','w')
#f.write(' atl ind pac arc sou\n')
#pl.savetxt(f,cormat,fmt='%9.5f',delimiter=' ')
#f.close()
#
#f = open(sheddir+'variability/flux_divq_sigs_iv.csv','w')
#f.write(' atl ind pac arc sou\n')
#pl.savetxt(f,sigmat,fmt='%9.5f',delimiter=' ')
#f.close()

amrcov = Covs(shedfluxes,divQ,'amr')
afrcov = Covs(shedfluxes,divQ,'afr')
eaacov = Covs(shedfluxes,divQ,'eaa')
aracov = Covs(shedfluxes,divQ,'ara')
aricov = Covs(shedfluxes,divQ,'ari')
arpcov = Covs(shedfluxes,divQ,'arp')
soacov = Covs(shedfluxes,divQ,'soa')
soicov = Covs(shedfluxes,divQ,'soi')
sopcov = Covs(shedfluxes,divQ,'sop')

covmat = pl.array([amrcov,afrcov,eaacov,aracov,aricov,arpcov,soacov,soicov,sopcov])

#f = open(sheddir+'variability/flux_divq_covs_iv.csv','w')
#f.write(' atl ind pac arc sou\n')
#pl.savetxt(f,covmat,fmt='%9.5f',delimiter=' ')
#f.close()

flxstd = pl.std(shedfluxes,axis=0)
dvqstd = pl.std(divQ,axis=0)

slopes = pl.zeros([9]); rsigs = pl.zeros([9])

basin=3.
divQ = divQ - pl.mean(divQ,axis=0)
shedfluxes = shedfluxes - pl.mean(shedfluxes,axis=0)

fig, ax = pl.subplots(3,1,figsize=(4,12)) # 8,8 for 2x2!
###############################################################################
#ax1 = pl.subplot(221)
#MoveAxes(ax1); pl.title('(a) Americas', fontsize=15,loc='left',y=1.01)
#slopes[0], rsigs[0] = ScatterPlot(divQ[:,basin],-shedfluxes[:,0])
###############################################################################
#ax2 = pl.subplot(222)
#MoveAxes(ax2); pl.title('(a) Africa', fontsize=15,loc='left',y=1.01)
#slopes[1], rsigs[1] = ScatterPlot(divQ[:,basin],shedfluxes[:,1])
###############################################################################
#ax3 = pl.subplot(222)
#MoveAxes(ax3); pl.title('(b) SE Asia', fontsize=15,loc='left',y=1.01)
#slopes[2], rsigs[2] = ScatterPlot(divQ[:,basin],shedfluxes[:,2])
###############################################################################
#ax4 = pl.subplot(311)
#MoveAxes(ax4); pl.title('(a) Arctic Atlantic', fontsize=15,loc='left',y=1.01)
#slopes[3], rsigs[3] = ScatterPlot(divQ[:,basin],shedfluxes[:,3])
#################################################################################
#ax5 = pl.subplot(312)
#MoveAxes(ax5); pl.title('(b) Arctic Indian', fontsize=15,loc='left',y=1.01)
#slopes[4], rsigs[4] = ScatterPlot(divQ[:,basin],shedfluxes[:,4])
################################################################################
#ax6 = pl.subplot(313)
#MoveAxes(ax6); pl.title('(c) Arctic Pacific', fontsize=15,loc='left',y=1.01)
#slopes[5], rsigs[5] = ScatterPlot(divQ[:,basin],shedfluxes[:,5])
###############################################################################
ax7 = pl.subplot(311)
MoveAxes(ax7); pl.title('(a) Southern Atlantic', fontsize=15,loc='left',y=1.01)
slopes[6], rsigs[6] = ScatterPlot(divQ[:,basin],-shedfluxes[:,6])
##############################################################################
ax8 = pl.subplot(312)
MoveAxes(ax8); pl.title('(b) Southern Indian', fontsize=15,loc='left',y=1.01)
slopes[7], rsigs[7] = ScatterPlot(divQ[:,basin],-shedfluxes[:,7])
##############################################################################
ax9 = pl.subplot(313)
MoveAxes(ax9); pl.title('(c) Southern Pacific', fontsize=15,loc='left',y=1.01)
slopes[8], rsigs[8] = ScatterPlot(divQ[:,basin],-shedfluxes[:,8])
##############################################################################

pl.subplots_adjust(left=0.04,right=0.96,top=0.90,bottom=0.06)
pl.suptitle('Southern $-$div$\mathbf{Q}$ vs. $\overline{\mathbf{Q}}\cdot\hat{\mathbf{n}}$',fontsize=20)
#pl.savefig(sheddir+'variability/divq_vs_flux_sou_2x2_rmvmn.png')

#f = open(sheddir+'variability/regress_sou.csv','w')
#pl.savetxt(f,zip(slopes,rsigs),fmt='%9.5f',delimiter=' ')
#f.close()


pl.figure(figsize=(8,8))
gs = gridspec.GridSpec(2, 4)
ig = [gs[0,:2],gs[0,2:],gs[1,1:3]]#[0,0,1]; iy = [:2,2:,1:3]

ax1 = pl.subplot(ig[0])
MoveAxes(ax1); pl.title('(a) Atlantic sector',fontsize=15,loc='left',y=1.01)
slopes[3], rsigs[3] = ScatterPlot(divQ[:,basin],shedfluxes[:,3])

ax2 = pl.subplot(ig[1])
MoveAxes(ax2); pl.title('(b) Indian sector',fontsize=15,loc='left',y=1.01)
slopes[4], rsigs[4] = ScatterPlot(divQ[:,basin],shedfluxes[:,4])

ax3 = pl.subplot(ig[2])
MoveAxes(ax3); pl.title('(c) Pacific sector',fontsize=15,loc='left',y=1.01)
slopes[5], rsigs[5] = ScatterPlot(divQ[:,basin],shedfluxes[:,5])

pl.suptitle('Arctic $-$div$\mathbf{Q}$ vs. $\overline{\mathbf{Q}}\cdot\hat{\mathbf{n}}$',fontsize=20)
pl.subplots_adjust(left=0.05,right=0.95,wspace=0.5)
#pl.savefig(sheddir+'variability/divq_vs_flux_arc_2x2_rmvmn.png')