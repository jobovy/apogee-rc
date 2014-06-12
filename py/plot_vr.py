#Plot vlos in the direction of l=0 and l=180
import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import apogee.tools.read as apread
import pixelize_sample
from plot_2dkinematics import vlosgal, dvlosgal
from plot_psd import _ADDLLOGGCUT, \
    _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
def plot_vr(plotfilename):
    #Read the APOGEE-RC data and pixelate it
    #APOGEE-RC
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    data= data[indx]
    #Get velocity field
    xmin, xmax= 5.5, 13.5
    dx= 1.
    pix= pixelize_sample.pixelXY(data,
                                 xmin=xmin,xmax=xmax,
                                 ymin=-dx/2.,ymax=dx/2.,
                                 dx=dx,dy=dx)
#                                 dx=_RCDX,dy=_RCDX)
    vr= pix.plot(lambda x: dvlosgal(x),
                 returnz=True,justcalc=True)
    vrunc= pix.plot(lambda x: dvlosgal(x),
                    func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                    returnz=True,justcalc=True)
    sr= pix.plot(lambda x: dvlosgal(x),
                 func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x))),
                 returnz=True,justcalc=True)
    srunc= pix.plot(lambda x: dvlosgal(x),
                    func=disperror,
                    returnz=True,justcalc=True)
    #print numpy.median(vr.flatten()[numpy.array([True,True,False,True,True,True,True,True,True],dtype='bool')])
    print vr.flatten()
    print vrunc.flatten()
    print sr.flatten()
    print srunc.flatten()
    rs= numpy.arange(xmin+dx/2.,xmax-dx/2.+0.00001,dx)
    print rs
    bovy_plot.bovy_print()
    srAxes= pyplot.axes([0.1,0.5,0.8,0.4])
    vrAxes= pyplot.axes([0.1,0.1,0.8,0.4])
    pyplot.sca(srAxes)
    pyplot.errorbar(rs,sr.flatten(),yerr=srunc.flatten(),
                    marker='o',ls='none',ms=6.,color='k')
    pyplot.xlim(0.,15.)
    pyplot.ylim(9.5,49.)
    #srAxes.set_yscale('log')
    bovy_plot._add_ticks(yticks=False)
    bovy_plot._add_axislabels(r'$ $',
                              r'$\sigma_R\,(\mathrm{km\,s}^{-1})$')
    nullfmt   = NullFormatter()         # no labels
    srAxes.xaxis.set_major_formatter(nullfmt)
    pyplot.sca(vrAxes)
    pyplot.errorbar(rs,vr.flatten(),yerr=vrunc.flatten(),
                    marker='o',ls='none',ms=6.,color='k')
    pyplot.plot([0.,20.],numpy.median(vr.flatten())*numpy.ones(2),'k--')
    bovy_plot._add_ticks()
    pyplot.xlim(0.,15.)
    pyplot.ylim(-14.,14.)
    bovy_plot._add_axislabels(r'$R\,(\mathrm{kpc})$',
                              r'$\langle V_R\rangle\,(\mathrm{km\,s}^{-1})$')
    bovy_plot.bovy_end_print(plotfilename)
    return None

def disperror(x):
    gaussdisp= 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))
    varerr2= 2.*gaussdisp**4./(len(x)-1.)
    return 0.5*varerr2**0.5/gaussdisp

if __name__ == '__main__':
    plot_vr(sys.argv[1])
