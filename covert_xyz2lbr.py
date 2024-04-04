import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import pprint as pp #prints 3D list pretty
from scipy.integrate import trapz, cumtrapz
from numpy import pi, sqrt, cos, sin, arctan2, linspace, log, exp, sinh, arcsinh
from astropy.io import fits
import sys

#Conversions from (x,y,z) to (l,b) from Mattia Sormani#

###################################
# xyz are cartesian coordinates centered at the GC
# XYZ are cartesian coordinates centered at the Sun position.
    #X points to the GC, Y points in the -(X \times z) direction,
    #Z points is obtained from the right-hand rule from the previous two.
# lbr are spherical coordinates corresponding to XYZ centered at the Sun position
    #(i.e., the usual Galactic coordinates)
# the Sun position is assumed to be (xsun,ysun,zsun)
# the Sun velocity is assumed to be (vxsun,vysun,vzsun)
# the following functions convert back and forth between the above coordinates
###################################
pi = np.pi
sin=np.sin
cos=np.cos
tan = np.tan
def XYZhat(xsun,ysun,zsun):
 # returns: Xhat, Yhat, Zhat components in the xyz frame
    rsun = np.sqrt(xsun**2+ysun**2+zsun**2)
    Rsun = np.sqrt(xsun**2+ysun**2)
    sintheta = Rsun/rsun
    costheta = zsun/rsun
    sinphi   = ysun/Rsun
    cosphi   = xsun/Rsun
    Xhat = np.array([-cosphi*sintheta,-sinphi*sintheta,-costheta])
    Yhat = np.array([+sinphi,         -cosphi,         +0.0     ])
    Zhat = np.array([-cosphi*costheta,-sinphi*costheta,+sintheta])
    return Xhat, Yhat, Zhat

def rblhat(l,b,r):
 # returns: rhat, bhat, lhat components in the xyz frame
    theta = pi/2 - b
    rhat  = [+sin(theta)*cos(l),+sin(theta)*sin(l),+cos(theta)]
    bhat  = [-cos(theta)*cos(l),-cos(theta)*sin(l),+sin(theta)]
    lhat  = [-sin(l),cos(l),0]
    return rhat, bhat, lhat


def xyz2XYZ(x,y,z,vx,vy,vz,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
    Xhat,Yhat,Zhat = XYZhat(xsun,ysun,zsun)
    Deltax  = x-xsun
    Deltay  = y-ysun
    Deltaz  = z-zsun
    Deltavx = vx-vxsun
    Deltavy = vy-vysun
    Deltavz = vz-vzsun
    X  = Deltax *Xhat[0] + Deltay *Xhat[1] + Deltaz *Xhat[2]
    Y  = Deltax *Yhat[0] + Deltay *Yhat[1] + Deltaz *Yhat[2]
    Z  = Deltax *Zhat[0] + Deltay *Zhat[1] + Deltaz *Zhat[2]
    vX = Deltavx*Xhat[0] + Deltavy*Xhat[1] + Deltavz*Xhat[2]
    vY = Deltavx*Yhat[0] + Deltavy*Yhat[1] + Deltavz*Yhat[2]
    vZ = Deltavx*Zhat[0] + Deltavy*Zhat[1] + Deltavz*Zhat[2]
    return X,Y,Z,vX,vY,vZ

def XYZ2xyz(X,Y,Z,vX,vY,vZ,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
    Xhat,Yhat,Zhat = XYZhat(xsun,ysun,zsun)
    Deltax  = X *Xhat[0] + Y *Yhat[0] + Z *Zhat[0]
    Deltay  = X *Xhat[1] + Y *Yhat[1] + Z *Zhat[1]
    Deltaz  = X *Xhat[2] + Y *Yhat[2] + Z *Zhat[2]
    Deltavx = vX*Xhat[0] + vY*Yhat[0] + vZ*Zhat[0]
    Deltavy = vX*Xhat[1] + vY*Yhat[1] + vZ*Zhat[1]
    Deltavz = vX*Xhat[2] + vY*Yhat[2] + vZ*Zhat[2]
    x  = Deltax + xsun
    y  = Deltay + ysun
    z  = Deltaz + zsun
    vx = Deltavx + vxsun
    vy = Deltavy + vysun
    vz = Deltavz + vzsun
    return x,y,z,vx,vy,vz


def XYZ2lbr(X,Y,Z,vX,vY,vZ):
    r     = np.sqrt(X**2+Y**2+Z**2)
    l     = np.arctan2(Y,X)
    theta = np.arccos(Z/r)
    b     = pi/2 - theta
    rhat, bhat, lhat = rblhat(l,b,r)
    vr = vX*rhat[0] + vY*rhat[1] + vZ*rhat[2]
    vb = vX*bhat[0] + vY*bhat[1] + vZ*bhat[2]
    vl = vX*lhat[0] + vY*lhat[1] + vZ*lhat[2]
    return l,b,r,vl,vb,vr

def lbr2XYZ(l,b,r,vl,vb,vr):
    theta = pi/2 - b
    rhat, bhat, lhat = rblhat(l,b,r)
    X  = r*sin(theta)*cos(l)
    Y  = r*sin(theta)*sin(l)
    Z  = r*cos(theta)
    vX = vr*rhat[0] + vl*lhat[0] + vb*bhat[0]
    vY = vr*rhat[1] + vl*lhat[1] + vb*bhat[1]
    vZ = vr*rhat[2] + vl*lhat[2] + vb*bhat[2]
    return X,Y,Z,vX,vY,vZ

def xyz2lbr(x,y,z,vx,vy,vz,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
    X,Y,Z,vX,vY,vZ = xyz2XYZ(x,y,z,vx,vy,vz,xsun,ysun,zsun,vxsun,vysun,vzsun)
    l,b,r,vl,vb,vr = XYZ2lbr(X,Y,Z,vX,vY,vZ)
    return l,b,r,vl,vb,vr

def lbr2xyz(l,b,r,vl,vb,vr,xsun=-8.0,ysun=0.0,zsun=0.0,vxsun=0.0,vysun=2.2,vzsun=0.0):
    X,Y,Z,vX,vY,vZ = lbr2XYZ(l,b,r,vl,vb,vr)
    x,y,z,vx,vy,vz = XYZ2xyz(X,Y,Z,vX,vY,vZ,xsun,ysun,zsun,vxsun,vysun,vzsun)
    return x,y,z,vx,vy,vz



