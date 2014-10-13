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
import galpy_simulations
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
    tmax= numpy.unravel_index(numpy.argmax(psd2d),psd2d.shape)
    tmax0= float(psd2d.shape[0]/2-tmax[0])/psd2d.shape[0]*2
    tmax1= float(tmax[1]-psd2d.shape[1]/2)/psd2d.shape[1]*2
    kmax= 1./_RCDX
    print tmax0*kmax, tmax1*kmax
    #kmax= numpy.amax(numpy.fft.fftfreq(resv.shape[0]*2,_RCDX))
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(psd2d.T,origin='lower',cmap='jet',
                          interpolation='nearest',
                          xrange=[-kmax,kmax],
                          yrange=[-kmax,kmax],
                          xlabel=r'$k_x\,(\mathrm{kpc}^{-1})$',
                          ylabel=r'$k_y\,(\mathrm{kpc}^{-1})$')
    bovy_plot.bovy_end_print(plotfilename)
    if True:
        spvlos= galpy_simulations.vlos('../sim/bar_rect_alpha0.015_hivres.sav')[1::2,1::2]
        potscale= 1.
k        spvlos+= numpy.random.normal(size=spvlos.shape)*resvunc/220./potscale
        simpsd2d= bovy_psd.psd2d(spvlos*220.*potscale)
        bovy_plot.bovy_print()
        bovy_plot.bovy_dens2d(simpsd2d.T,
                              origin='lower',cmap='jet', 
                              interpolation='nearest',
                              xrange=[-kmax,kmax],
                              yrange=[-kmax,kmax],
                              xlabel=r'$k_x\,(\mathrm{kpc}^{-1})$',
                              ylabel=r'$k_y\,(\mathrm{kpc}^{-1})$')
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
                              xrange=[_RCXMIN,_RCXMAX],
                              yrange=[_RCYMIN,_RCYMAX],
                              xlabel=r'$X_{\mathrm{GC}}\,(\mathrm{kpc})$',
                              ylabel=r'$Y_{\mathrm{GC}}\,(\mathrm{kpc})$',
                              vmin=-16.,vmax=16.,
                              zlabel=r'$\Delta V^{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
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
    numpy.random.seed(17) #17, went to 30 / 23 good w/ diff gamma
    plot_psd2d(sys.argv[1])
