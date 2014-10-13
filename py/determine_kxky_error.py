#Determine the error on (kx,ky) by using the best-fitting bar as a noiseless model, adding noise, and seeing how the peak shifts
import numpy
import apogee.tools.read as apread
from galpy.util import bovy_plot
import bovy_psd
import pixelize_sample
import galpy_simulations
from plot_psd import _ADDLLOGGCUT, _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
_HIVRESSTR= '_hivres'
def determine_kxky_error():
    #Load fiducial bar model
    spvlos= galpy_simulations.vlos('../sim/bar_rect_alpha0.015%s.sav' % _HIVRESSTR)[1::2,1::2]
    #spvlos= galpy_simulations.vlos('../sim/spiral_rect_omegas0.33_alpha-14%s.sav' % _HIVRESSTR)[1::2,1::2]*0.85
    #spvlos= galpy_simulations.vlos('../sim/spiral_rect_omegas0.33_alpha-7%s.sav' % _HIVRESSTR)[1::2,1::2]*0.25
    psd2d= bovy_psd.psd2d(spvlos)
    kmax= 1./_RCDX
    #Calculate maximum
    psd2d[psd2d.shape[0]/2-1:psd2d.shape[0]/2+2,psd2d.shape[1]/2-1:psd2d.shape[1]/2+2]= 0.
    tmax= numpy.unravel_index(numpy.argmax(psd2d),psd2d.shape)
    tmax0= float(psd2d.shape[0]/2-tmax[0])/psd2d.shape[0]*2
    tmax1= float(tmax[1]-psd2d.shape[1]/2)/psd2d.shape[1]*2
    print tmax0*kmax, tmax1*kmax
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(psd2d.T,origin='lower',cmap='jet',
                          interpolation='nearest',
                          xrange=[-kmax,kmax],
                          yrange=[-kmax,kmax],
                          xlabel=r'$k_x\,(\mathrm{kpc}^{-1})$',
                          ylabel=r'$k_y\,(\mathrm{kpc}^{-1})$')
    bovy_plot.bovy_end_print('/Users/bovy/Desktop/test.png')
    #Read the data for the noise
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    print "Using %i stars for low-Z 2D kinematics analysis" % numpy.sum(indx)
    data= data[indx]
    #Get residuals
    dx= _RCDX
    pix= pixelize_sample.pixelXY(data,
                                 xmin=_RCXMIN,xmax=_RCXMAX,
                                 ymin=_RCYMIN,ymax=_RCYMAX,
                                 dx=dx,dy=dx)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    #Now do 1000 MC simulations to determine the error on kmax
    nmc= 1000
    kxmax= numpy.zeros(nmc)
    kymax= numpy.zeros(nmc)
    for ii in range(nmc):
        newresv= spvlos+numpy.random.normal(size=spvlos.shape)*resvunc/220.
        simpsd2d= bovy_psd.psd2d(newresv*220.)
        simpsd2d[simpsd2d.shape[0]/2-1:simpsd2d.shape[0]/2+2,simpsd2d.shape[1]/2-1:simpsd2d.shape[1]/2+2]= 0.
        tmax= numpy.unravel_index(numpy.argmax(simpsd2d),psd2d.shape)
        tmax0= float(psd2d.shape[0]/2-tmax[0])/psd2d.shape[0]*2
        tmax1= float(tmax[1]-psd2d.shape[1]/2)/psd2d.shape[1]*2
        kmax= 1./_RCDX
        kxmax[ii]= tmax0*kmax
        kymax[ii]= tmax1*kmax
    print numpy.mean(kxmax), numpy.std(kxmax)
    print numpy.mean(kymax), numpy.std(kymax)
    return None

if __name__ == '__main__':
    determine_kxky_error()
