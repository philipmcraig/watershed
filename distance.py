# -*- coding: utf-8 -*-
"""
Created on Sun Dec 11 14:21:46 2016

@author: np838619
"""

import pylab as pl

def Vicenty(pt1,pt2):
    """
    """
    pt1 = pl.radians(pt1); pt2 = pl.radians(pt2)
    
    phi1 = pt1[1]; phi2 = pt2[1]
    lam1 = pt1[0]; lam2 = pt2[0]; dlam = lam2 - lam1
    
    top = (pl.cos(phi2)*pl.sin(dlam))**2 +\
        (pl.cos(phi1)*pl.sin(phi2) - pl.sin(phi1)*pl.cos(phi2)*pl.cos(dlam))**2
    
    bot = pl.sin(phi1)*pl.sin(phi2) + pl.cos(phi1)*pl.cos(phi2)*pl.cos(dlam)
    
    dsig = pl.arctan(pl.sqrt(top)/bot)
    
    R = 6.37e6
    
    d = R*dsig
    
    return d

def Haversine(pt1,pt2):
    """
    """
    pt1 = pl.radians(pt1); pt2 = pl.radians(pt2)
    
    phi1 = pt1[1]; phi2 = pt2[1]; dphi = phi2 - phi1
    lam1 = pt1[0]; lam2 = pt2[0]; dlam = lam2 - lam1
    
    dsig = 2*pl.arcsin(pl.sqrt((pl.sin(dphi/2)**2 + pl.cos(phi1)*pl.cos(phi2)*pl.sin(dlam/2)**2)))
    
    R = 6.37e6
    
    d = R*dsig
    
    return d