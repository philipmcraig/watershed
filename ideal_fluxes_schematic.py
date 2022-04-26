# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 12:09:17 2018

@author: np838619
"""

from __future__ import division
import pylab as pl

def TidyUp(axx):
    """
    """
    axx.spines['right'].set_color('none')
    axx.spines['top'].set_color('none')
    axx.spines['left'].set_color('none')
    axx.spines['bottom'].set_color('none')
    axx.xaxis.set_major_formatter(pl.NullFormatter())
    axx.yaxis.set_major_formatter(pl.NullFormatter())
    pl.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        left='off',
        top='off',         # ticks along the top edge are off
        right='off',
        labelbottom='off') # labels along the bottom edge are off

    return None

def Vlines(axx):
    """
    """
    axx.axvline(x=0.05,color='k',lw=3); axx.axvline(x=0.95,color='k',lw=3)
    axx.axvline(x=0.35,color='g',lw=3); axx.axvline(x=0.65,color='b',lw=3)

    return None

def Oceans(axx):
    """
    """
    axx.annotate('Indian Ocean',(0.1,0.9),fontsize=18)
    axx.annotate('Pacific Ocean',(0.4,0.9),fontsize=18)
    axx.annotate('Atlantic Ocean',(0.7,0.9),fontsize=18)

pl.close('all')

fig,ax = pl.subplots(3,1,figsize=(10,10))

ax1 = pl.subplot(311)
TidyUp(ax1); Vlines(ax1); Oceans(ax1)
ax1.arrow(0.075,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax1.arrow(0.375,0.5,-0.05,0,fc="g", ec="g", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax1.arrow(0.675,0.5,-0.05,0,fc="b", ec="b", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax1.arrow(0.975,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax1.text(0.09,0.47,'0.36',color='k',fontsize=20)
ax1.text(0.39,0.47,'0.50',color='g',fontsize=20)
ax1.text(0.69,0.47,'0.24',color='b',fontsize=20)
ax1.text(0.99,0.47,'0.36',color='k',fontsize=20)
ax1.text(0.175,0.25,'-0.18',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
ax1.text(0.475,0.25,'-0.65',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
ax1.text(0.775,0.25,'-0.28',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
pl.title('(a) Zonally averaged moisture fluxes',fontsize=20)

ax2 = pl.subplot(312)
TidyUp(ax2); Vlines(ax2); Oceans(ax2)
ax2.arrow(0.075,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax2.arrow(0.375,0.5,-0.05,0,fc="g", ec="g", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax2.arrow(0.675,0.5,-0.05,0,fc="b", ec="b", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax2.arrow(0.975,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax2.text(0.09,0.47,'0.16',color='k',fontsize=20)
ax2.text(0.39,0.47,'0.50',color='g',fontsize=20)
ax2.text(0.69,0.47,'0.23',color='b',fontsize=20)
ax2.text(0.99,0.47,'0.16',color='k',fontsize=20)
ax2.text(0.175,0.25,'0.02',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
ax2.text(0.475,0.25,'-0.66',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
ax2.text(0.775,0.25,'-0.47',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
pl.title('(b) Change African and American fluxes',fontsize=20)


#ax2.arrow(0.075,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
#         head_length=0.01,width=0.02,head_starts_at_zero=False)
#ax2.arrow(0.375,0.5,-0.05,0,fc="g", ec="g", linewidth = 2, head_width=0.06,
#         head_length=0.01,width=0.02,head_starts_at_zero=False)
#ax2.arrow(0.675,0.5,-0.05,0,fc="b", ec="b", linewidth = 2, head_width=0.06,
#         head_length=0.01,width=0.02,head_starts_at_zero=False)
#ax2.arrow(0.975,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
#         head_length=0.01,width=0.02,head_starts_at_zero=False)
#ax2.text(0.09,0.47,'0.17',color='k',fontsize=20)
#ax2.text(0.39,0.47,'0.50',color='g',fontsize=20)
#ax2.text(0.69,0.47,'0.23',color='b',fontsize=20)
#ax2.text(0.99,0.47,'0.17',color='k',fontsize=20)
#ax2.text(0.175,0.25,'0.02',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
#ax2.text(0.475,0.25,'-0.66',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
#ax2.text(0.775,0.25,'-0.47',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
#pl.title('(b) Change South-East Asia',fontsize=20)

ax3 = pl.subplot(313)
TidyUp(ax3); Vlines(ax3); Oceans(ax3)

ax3.arrow(0.075,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax3.arrow(0.325,0.5,0.05,0,fc="g", ec="g", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax3.arrow(0.675,0.5,-0.05,0,fc="b", ec="b", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax3.arrow(0.975,0.5,-0.05,0,fc="k", ec="k", linewidth = 2, head_width=0.06,
         head_length=0.01,width=0.02,head_starts_at_zero=False)
ax3.text(0.09,0.47,'0.16',color='k',fontsize=20)
ax3.text(0.39,0.47,'0.16',color='g',fontsize=20)
ax3.text(0.69,0.47,'0.23',color='b',fontsize=20)
ax3.text(0.99,0.47,'0.16',color='k',fontsize=20)
ax3.text(0.175,0.25,'-0.64',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
ax3.text(0.475,0.25,'0.00',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
ax3.text(0.775,0.25,'-0.47',bbox={'facecolor':'white', 'alpha':1.0, 'pad':10},fontsize=20)
pl.title('(c) ERA-Interim/Change South-East Asian flux',fontsize=20)

#ax4 = pl.subplot(414)
#TidyUp(ax4); Vlines(ax4); Oceans(ax4)


pl.subplots_adjust(top=0.95,bottom=0.04,left=0.06,right=0.94,hspace=0.27)
fig.set_edgecolor('k')

pl.savefig('/home/users/qx911590/np838619/Watershed/ideal_fluxes_new.png',dpi=350)