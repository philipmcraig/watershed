# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 14:29:12 2018

@author: np838619
"""

pl.close('all')

F = pl.array([shedfluxes[:,2],-shedfluxes[:,0],divQ[:,2]])
T = ['(a) SE Asia $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
     '(b) Americas $\mathbf{Q}\cdot\mathbf{\hat{n}}$',
#    '(c) Southern Atlantic $\mathbf{Q}\cdot\mathbf{\hat{n}}$'
   '(c) Pacific $P-E$']
proj = ccrs.PlateCarree(central_longitude=-165)
ext = [90,300,-40,70]
lons,lats = pl.meshgrid(eralon,eralat)

pl.figure(figsize=(12,8))
gs = gridspec.GridSpec(2, 4)
ig = [gs[0,:2],gs[0,2:],gs[1,1:3]]#[0,0,1]; iy = [:2,2:,1:3]
for i in range(3):
    x,y = Regression(F[i]*10,pac_tcdq,eralat,eralon)
    axx = pl.subplot(ig[i],projection=proj)
    cs = PlotStats(axx,x,y,lons,lats,proj,ext,levels,norm,False)
    pl.title(T[i],fontsize=18,loc='left')
    gx = axx.gridlines(draw_labels=True)
    if i == 0:
        gx.ylabels_right = False; gx.xlabels_top = False
    elif i == 1:
        gx.ylabels_left = False; gx.xlabels_top = False
    elif i == 2:
        gx.xlabels_top = False
    gx.xlocator = mticker.FixedLocator([-150,-90,-60,90,150,210])
    gx.ylocator = mticker.FixedLocator([-40,-20,0,20,40,60,80])
    gx.xformatter = LONGITUDE_FORMATTER; gx.yformatter = LATITUDE_FORMATTER
    gx.xlabel_style = {'size': 13}; gx.ylabel_style = {'size': 13}

f = pl.gcf()
colax = f.add_axes([0.13,0.08,0.76,0.03])
cb=pl.colorbar(cs,cax=colax,orientation='horizontal')
cb.set_label('mm/day/dSv',fontsize=18)
cb.ax.tick_params(labelsize=14)
pl.subplots_adjust(top=0.95,bottom=0.15,left=0.05,right=0.95,wspace=0.24)
#pl.savefig(sheddir+'variability/sss_reg_maps/ERA_tcdq/pac_tcdq_panel.png')