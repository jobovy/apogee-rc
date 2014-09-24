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
    return (saved[1],saved[2]-saved[8])

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
    return vlos[:,:,0]
    
