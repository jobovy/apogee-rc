import numpy
from galpy.util import bovy_coords
from galpy.potential import SteadyLogSpiralPotential
def simulate_vlos_spiral(alpha=-7.,m=2.,gamma=0.7853981633974483,omegas=0.65,
                         xmin=None,xmax=None,ymin=None,ymax=None,
                         dx=None):
    """Perform a simple simulation of the effect of spiral structure on the velocity field; use the WKB approximation"""
    #Set up galpy spiral
    sp= SteadyLogSpiralPotential(alpha=alpha,m=m,gamma=gamma,omegas=omegas)
    #Now evaluate on the grid
    nx= int(numpy.ceil((xmax-xmin)/dx))
    ny= int(numpy.ceil((ymax-ymin)/dx))
    xs= numpy.arange(xmin+dx/2.,xmax-dx/2.+0.000001,dx)
    ys= numpy.arange(ymin+dx/2.,ymax-dx/2.+0.000001,dx)
    xv,yv= numpy.meshgrid(xs,ys,indexing='xy')
    rs= numpy.sqrt(xv**2.+yv**2.)/8.
    phis= numpy.arctan2(yv,xv)
    (d,l)= bovy_coords.rphi_to_dl_2d(rs,phis)
    cospl= numpy.cos(phis+l)
    sinpl= numpy.sin(phis+l)
    vlos= numpy.empty((nx,ny))
    Delta= 2.-m**2.*(omegas-1.)**2.
    for ii in range(nx):
        for jj in range(ny):
            pot= sp(rs[ii,jj],phi=phis[ii,jj])
            vr= m*(omegas-1.)/Delta*alpha/1.*reduction_factor(alpha/1.)*pot
            vt= -2.*-0.5/Delta*alpha/1.*reduction_factor(alpha/1.)*pot
            vlos[ii,jj]= -vr*cospl[ii,jj]+vt*sinpl[ii,jj]
    return vlos

def reduction_factor(k):
    return 1.
