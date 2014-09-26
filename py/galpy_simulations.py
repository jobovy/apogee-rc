import os
import cPickle as pickle
import numpy
from galpy.util import bovy_coords
def read_sim(savefilename):
    savefile= open(savefilename,'rb')
    smass= pickle.load(savefile)
    meanvr= pickle.load(savefile)
    meanvt= pickle.load(savefile)
    sigmar2= pickle.load(savefile)
    sigmat2= pickle.load(savefile)
    sigmart= pickle.load(savefile)
    vertexdev= pickle.load(savefile)
    surfmass_init= pickle.load(savefile)
    meanvt_init= pickle.load(savefile)
    sigmar2_init= pickle.load(savefile)
    sigmat2_init= pickle.load(savefile)
    savefile.close()
    return (smass,meanvr,meanvt,sigmar2,sigmat2,sigmart,vertexdev,
            surfmass_init,meanvt_init,sigmar2_init,sigmat2_init)

def read_dmeanvrvt(savefile):
    saved= read_sim(savefile)
    return (saved[1][:,:,0],saved[2][:,:,0]-saved[8])

def vlos(savefile):
    dvr,dvt= read_dmeanvrvt(savefile)
    resx, resy= dvr.shape[0], dvr.shape[1]
    #Project onto the line-of-sight
    from plot_psd import _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
    xgrid= numpy.linspace((_RCXMIN-8.)/8.+_RCDX/8./2.,
                          (_RCXMAX-8.)/8.-_RCDX/8./2.,
                          resx)
    ygrid= numpy.linspace(_RCYMIN/8.+_RCDX/8./2.,
                          _RCYMAX/8.-_RCDX/8./2.,
                          resy)
    xv,yv= numpy.meshgrid(xgrid,ygrid,indexing='ij')
    rs= numpy.sqrt((1.+xv)**2.+yv**2.)
    phis= numpy.arctan2(yv,1.+xv)
    (d,l)= bovy_coords.rphi_to_dl_2d(rs,phis)
    cospl= numpy.cos(phis+l)
    sinpl= numpy.sin(phis+l)
    vlos= -cospl*dvr+sinpl*dvt
    return vlos[:,:]
    
def vlos_altrect(savefile):
    dvr,dvt= read_dmeanvrvt(savefile)
    resx, resy= dvr.shape[0], dvr.shape[1]
    #Project onto the line-of-sight
    from plot_psd import _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
    xgrid= numpy.linspace((_RCXMIN-2.25-8.)/8.+_RCDX/8./2.,
                          (_RCXMAX+2.25-8.)/8.-_RCDX/8./2.,
                          resx)
    ygrid= numpy.linspace(_RCYMIN/8.-2.25/8.+_RCDX/8./2.,
                          _RCYMAX/8.+2.25/8.-_RCDX/8./2.,
                          resy)
    xv,yv= numpy.meshgrid(xgrid,ygrid,indexing='ij')
    rs= numpy.sqrt((1.+xv)**2.+yv**2.)
    phis= numpy.arctan2(yv,1.+xv)
    (d,l)= bovy_coords.rphi_to_dl_2d(rs,phis)
    cospl= numpy.cos(phis+l)
    sinpl= numpy.sin(phis+l)
    vlos= -cospl*dvr+sinpl*dvt
    return vlos[:,:]
    
def vlos_polar(savefile):
    dvr,dvt= read_dmeanvrvt(savefile)
    resx, resy= dvr.shape[0], dvr.shape[1]
    #Project onto the line-of-sight
    xgrid= numpy.linspace(0.,2.*numpy.pi*(1.-1./resy/2.),resx)
    ygrid= numpy.linspace(0.5,2.,resy)
    phis,rs= numpy.meshgrid(xgrid,ygrid,indexing='ij')
    (d,l)= bovy_coords.rphi_to_dl_2d(rs,phis)
    cospl= numpy.cos(phis+l)
    sinpl= numpy.sin(phis+l)
    vlos= -cospl*dvr+sinpl*dvt
    return vlos[:,:]
    
def vlos_elliptical(res=19,cp=0.05,sp=0.,p=0.,beta=0.,xgrid=None,ygrid=None):
    sr= 31.4
    dvramp= 1./(1.-beta)*(1.+0.5*p)*(1.-8.*(sr/220.)**2.)
    dvtamp= 1./(1.-beta)*(1.+0.25*p*(1.+beta))*(1.-2.5*(sr/220.)**2.)
    #Project onto the line-of-sight
    if xgrid is None:
        from plot_psd import _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
        xgrid= numpy.linspace((_RCXMIN-8.)/8.+_RCDX/8./2.,
                              (_RCXMAX-8.)/8.-_RCDX/8./2.,
                              res)
        ygrid= numpy.linspace(_RCYMIN/8.+_RCDX/8./2.,
                              _RCYMAX/8.-_RCDX/8./2.,
                              res)
    xv,yv= numpy.meshgrid(xgrid,ygrid,indexing='ij')
    rs= numpy.sqrt((1.+xv)**2.+yv**2.)
    phis= numpy.arctan2(yv,1.+xv)
    (d,l)= bovy_coords.rphi_to_dl_2d(rs,phis)
    cospl= numpy.cos(phis+l)
    sinpl= numpy.sin(phis+l)
    dvr= -dvramp*(numpy.sin(2.*phis)*cp-numpy.cos(2.*phis)*sp)
    dvt= -dvtamp*(numpy.cos(2.*phis)*cp+numpy.sin(2.*phis)*sp)
    vlos= (-cospl*dvr+sinpl*dvt)*rs**(p-2.*beta)
    return vlos
    
