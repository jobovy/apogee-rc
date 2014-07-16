import sys
import re
import os, os.path
import numpy
from scipy import interpolate
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter
import fitsio
import apogee.tools.read as apread
import pixelize_sample
import bovy_psd
from plot_2dkinematics import dvlosgal
import hackGCS
import readAndHackHoltz
from simple_spiral_simulation import simulate_vlos_spiral
from plot_psd import _ADDLLOGGCUT, _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
#Parameters of the pixelizations
def plot_psd2d(plotfilename):
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    print "Using %i stars for low-Z 2D kinematics analysis" % numpy.sum(indx)
    data= data[indx]
    #Get residuals
    dx= _RCDX
    binsize= .8#.765
    pix= pixelize_sample.pixelXY(data,
                                 xmin=_RCXMIN,xmax=_RCXMAX,
                                 ymin=_RCYMIN,ymax=_RCYMAX,
                                 dx=dx,dy=dx)
    resv= pix.plot(lambda x: dvlosgal(x,vtsun=220.+22.6),
                   returnz=True,justcalc=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd2d= bovy_psd.psd2d(resv)
    print psd2d
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(psd2d.T,origin='lower',cmap='jet',
                          interpolation='nearest',
                          colorbar=True,shrink=0.78)
    bovy_plot.bovy_end_print(plotfilename)
    if True:
        alpha= -12.5 #-12.5
        spvlos= simulate_vlos_spiral(alpha=alpha,
                                     gamma=1.2,
                                     omegas=.65,
                                     xmin=_RCXMIN,xmax=_RCXMAX,
                                     ymin=_RCYMIN,ymax=_RCYMAX,
                                     dx=_RCDX)
        potscale= 1.35
        spvlos+= numpy.random.normal(size=spvlos.shape)*resvunc/220./potscale
        print numpy.arctan(2./alpha)/numpy.pi*180., numpy.sqrt(0.035/numpy.fabs(alpha)/2.)*potscale*220., numpy.sqrt(0.035/numpy.fabs(alpha))*potscale*220.
        simpsd2d= bovy_psd.psd2d(spvlos*220.*potscale)
        bovy_plot.bovy_print()
        bovy_plot.bovy_dens2d(simpsd2d.T,#(spvlos*220.*potscale).T,
                              origin='lower',cmap='jet', 
                              interpolation='nearest',
                              colorbar=True,shrink=0.78)
        fileparts= re.split(r'\.',plotfilename)
        nparts= len(fileparts)
        simfilename= ''
        for ii in range(nparts):
            if ii == nparts-2: simfilename+= fileparts[ii]+'_simpsd2d.'
            elif ii == nparts-1: simfilename+= fileparts[ii]
            else: simfilename+= fileparts[ii]+'.'
        bovy_plot.bovy_end_print(simfilename)
        bovy_plot.bovy_print()
        bovy_plot.bovy_dens2d((spvlos*220.*potscale).T,
                              origin='lower',cmap='jet', 
                              interpolation='nearest',
                              colorbar=True,shrink=0.78)
        fileparts= re.split(r'\.',plotfilename)
        nparts= len(fileparts)
        simfilename= ''
        for ii in range(nparts):
            if ii == nparts-2: simfilename+= fileparts[ii]+'_sim.'
            elif ii == nparts-1: simfilename+= fileparts[ii]
            else: simfilename+= fileparts[ii]+'.'
        bovy_plot.bovy_end_print(simfilename)
        return None

if __name__ == '__main__':
    plot_psd2d(sys.argv[1])
