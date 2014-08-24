import numpy
from galpy.util import bovy_coords
from galpy.potential import SteadyLogSpiralPotential
def simulate_vlos_spiral(alpha=-7.,m=2.,gamma=0.7853981633974483,omegas=0.65,
                         xmin=None,xmax=None,ymin=None,ymax=None,
                         dx=None,returnvrvt=False):
    """Perform a simple simulation of the effect of spiral structure on the velocity field; use the WKB approximation; if returnvrvt == True, return (vr,vt), otherwise return vlos"""
    #Set up galpy spiral
    sp= SteadyLogSpiralPotential(alpha=alpha,m=m,gamma=gamma,omegas=omegas)
    #Now evaluate on the grid
    nx= int(numpy.ceil((xmax-xmin)/dx))
    ny= int(numpy.ceil((ymax-ymin)/dx))
    xs= numpy.arange(xmin+dx/2.,xmax-dx/2.+0.000001,dx)
    ys= numpy.arange(ymin+dx/2.,ymax-dx/2.+0.000001,dx)
    xv,yv= numpy.meshgrid(xs,ys,indexing='ij')
    rs= numpy.sqrt(xv**2.+yv**2.)/8.
    phis= numpy.arctan2(yv,xv)
    (d,l)= bovy_coords.rphi_to_dl_2d(rs,phis)
    cospl= numpy.cos(phis+l)
    sinpl= numpy.sin(phis+l)
    vlos= numpy.empty((nx,ny))
    vr= numpy.empty((nx,ny))
    vt= numpy.empty((nx,ny))
    Delta= 2.-m**2.*(omegas-1.)**2.
    for ii in range(nx):
        for jj in range(ny):
            pot= sp(rs[ii,jj],phi=phis[ii,jj])
            potsin= sp._amp*sp._A/sp._alpha\
                *numpy.sqrt(1.-(sp(rs[ii,jj],phi=phis[ii,jj])/sp._A/sp._amp*sp._alpha)**2.)
            vr[ii,jj]= m*(omegas-1.)/Delta*rs[ii,jj]*alpha/1.\
                *reduction_factor(alpha/1.)*pot
            vt[ii,jj]= 2.*-0.5/Delta*rs[ii,jj]*alpha/1.\
                *reduction_factor(alpha/1.)*potsin
            vlos[ii,jj]= -vr[ii,jj]*cospl[ii,jj]+vt[ii,jj]*sinpl[ii,jj]
    if returnvrvt:
        return (vr,vt)
    else:
        return vlos

def reduction_factor(k):
    return 1.
