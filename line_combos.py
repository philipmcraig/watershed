# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 10:31:27 2017
Code to create one line for the Arctic Ocean watershed, a combination of the
Bering Strait to Greenland, Greenland to Scandinavia, and Scandinavia to Bering
Strait lines.

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap

def ArcticLine(string):
    """Function to read files as I need to do it 3 times.
    Input should be a string, either 'clicks', 'traj_release'
    read 3 files for each part of the Arctic line
    Make an empty array, stick the lines in and remove repeated end points
    Write new file here?
    """
    sheddir = '/home/np838619/Watershed/shed_defs/'
    BG = pl.genfromtxt(sheddir+'BSGl_'+string+'.txt',skip_header=5) # Bering Strait to Greenland
    GS = pl.genfromtxt(sheddir+'GlSc_'+string+'.txt',skip_header=5) # Greenland to Scandinavia
    SB = pl.genfromtxt(sheddir+'ScBS_'+string+'.txt',skip_header=5) # Scandinavia to Bering Strait lines
    
    if string == 'clicks':
        points = BG.shape[0] + GS.shape[0] + SB.shape[0]
        arc_line = pl.zeros([points-2,2])
        BG_end = BG.shape[0]; GS_end = BG.shape[0]+GS.shape[0] 
        arc_line[:BG_end] = BG[:] # all points of BG
        arc_line[BG_end:GS_end-2] = GS[1:-1] # first point of GS removed
        arc_line[GS_end-2:] = SB[:]
    elif string == 'traj_release':
        arc_line = []
        for pt in range(BG.shape[0]):
            arc_line.append(BG[pt])
        for pt in range(GS.shape[0]):
            arc_line.append(GS[pt])
        for pt in range(SB.shape[0]):
            arc_line.append(SB[pt])
        arc_line = pl.asarray(arc_line)
    
    for i in range(arc_line.shape[0]):
        if arc_line[i,0] > 180.:
            arc_line[i,0] = arc_line[i,0] - 360.
    
    f = open(sheddir+'Ar_'+string+'.txt','w')
    if string == 'clicks':
        f.write('Longitude, Latitude co-ordinates defining the Arctic Ocean ' + 
                'drainage divide. \nUsed to release back trajectories. \n \n' 
            + str(len(arc_line)) + '  points defined in this file. \n \n')
    else:
        f.write('Longitude, Latitude co-ordinates defining the Arctic Ocean' + 
        'Arctic Ocean drainage divide. \nUsed to release back trajectories. \n \n' + 
        str(len(arc_line)) + '  points defined in this file. \n \n')
    pl.savetxt(f,arc_line,fmt='%9.5f')
    f.close()
    
    return arc_line

def SouthernLine(string):
    """Function to read 5 files
    Input: 'clicks', 'traj_release'
    """
    sheddir = '/home/np838619/Watershed/shed_defs/'
    IO = pl.genfromtxt(sheddir+'SOI_'+string+'.txt',skip_header=5) # Indian Ocean
    OZ = pl.genfromtxt(sheddir+'Oz_'+string+'.txt',skip_header=5) # Australia
    PO = pl.genfromtxt(sheddir+'SOP_'+string+'.txt',skip_header=5) # Pacific Ocean
    LP = pl.genfromtxt(sheddir+'LP_'+string+'.txt',skip_header=5) # River Plate
    AO = pl.genfromtxt(sheddir+'SOA_'+string+'.txt',skip_header=5) # Atlantic Ocean
    
    if string == 'clicks':
        points = IO.shape[0] + OZ.shape[0] + PO.shape[0] + LP.shape[0] + AO.shape[0]
        sou_line = pl.zeros([points-4,2])
        IO_end = IO.shape[0]; OZ_end = IO_end + OZ.shape[0]; PO_end = OZ_end + PO.shape[0]
        LP_end = PO_end + LP.shape[0]; AO_end = LP_end + AO.shape[0]
        sou_line[:IO_end] = IO[:]
        sou_line[IO_end:OZ_end-1] = OZ[1:]
        sou_line[OZ_end-1:PO_end-2] = PO[1:]
        sou_line[PO_end-2:LP_end-4] = LP[1:-1]
        sou_line[LP_end-4:] = AO[:]
    elif string == 'traj_release':
        sou_line = []
        for pt in range(IO.shape[0]):
            sou_line.append(IO[pt])
        for pt in range(OZ.shape[0]):
            sou_line.append(OZ[pt])
        for pt in range(PO.shape[0]):
            sou_line.append(PO[pt])
        for pt in range(LP.shape[0]):
            sou_line.append(LP[pt])
        for pt in range(AO.shape[0]):
            sou_line.append(AO[pt])
        sou_line = pl.asarray(sou_line)
    
    for i in range(sou_line.shape[0]):
        if sou_line[i,0] > 180.:
            sou_line[i,0] = sou_line[i,0] - 360.
    
    f = open(sheddir+'SO_'+string+'.txt','w')
    if string == 'clicks':
        f.write('Longitude, Latitude co-ordinates defining the Southern Ocean ' +
                'drainage divide. \nUsed to release back trajectories. \n \n' 
            + str(len(sou_line)) + '  points defined in this file. \n \n')
    else:
        f.write('Longitude, Latitude co-ordinates defining the Southern Ocean ' 
        'drainage divide. \nUsed to release back trajectories. \n \n' + 
        str(len(sou_line)) + '  points defined in this file. \n \n')
    pl.savetxt(f,sou_line,fmt='%9.5f')
    f.close()
    
    return sou_line
    
def RR2(interp_pts):
    repeat = []
    for r in range(1,interp_pts.shape[0]):
        if interp_pts[r,0] == interp_pts[r-1,0] and interp_pts[r,1] == interp_pts[r-1,1]:
            repeat.append(r)
    
    RP2 = []#; L2 = []
    for r in range(interp_pts.shape[0]):
        if r not in repeat:
            RP2.append(interp_pts[r])
            #L2.append(labels[r])
    rlspts = pl.asarray(RP2)#; labs = pl.asarray(L2)
    
    return rlspts#, labs

pl.close('all')
sheddir = '/home/np838619/Watershed/shed_defs/'
# Call function 3 times, each using a different string
cl = ArcticLine('clicks')#SouthernLine('clicks')#
tj = ArcticLine('traj_release')#SouthernLine('traj_release')#
tjn = RR2(tj)#ArcticLine('traj_release_new')#

f = open(sheddir+'Ar_'+'traj_release_new'+'.txt','w')
f.write('Longitude, Latitude co-ordinates defining the Southern Ocean ' +
                'drainage divide. \nUsed to release back trajectories. \n \n' 
            + str(len(tjn)) + '  points defined in this file. \n \n')
pl.savetxt(f,tjn,fmt='%9.5f')
f.close()

m = Basemap(projection='nplaea',boundinglat=25,lon_0=0.,round=True)
m.drawcoastlines(); pl.tight_layout()#; m.drawparallels([-35])
#m.plot(cl[:,0],cl[:,1],'b',latlon=True)
#for i in range(len(cl)-1):
#    m.drawgreatcircle(cl[i,0],cl[i,1],cl[i+1,0],cl[i+1,1],color='b')
m.plot(cl[:,0],cl[:,1],'ro',latlon=True)
m.plot(tj[:,0],tj[:,1],'b',latlon=True)
pl.savefig(sheddir+'Ar_clicks.png')#;pl.close()
pl.figure(2)
m.drawcoastlines(); pl.tight_layout()
m.plot(tjn[:,0],tjn[:,1],'go',latlon=True)
pl.savefig(sheddir+'Ar_interp_pts.png')