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
_EXT='png'
_ADDLLOGGCUT= True
_ADDGCS= True
_ADDRED= False
_ADDRAVE= True
_NNOISE= 1000
_PLOTBAND= False
_SUBTRACTERRORS= 1.
def plot_psd():
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    print "Using %i stars for low-Z 2D kinematics analysis" % numpy.sum(indx)
    data= data[indx]
    #Get residuals
    dx= 0.75
    binsize= .8#.765
    pix= pixelize_sample.pixelXY(data,
                                 xmin=5.,xmax=12.5,
                                 ymin=-3.,ymax=4.5,
                                 dx=dx,dy=dx)
    resv= pix.plot(lambda x: dvlosgal(x),returnz=True,justcalc=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    scale= 3.25/numpy.median(numpy.sqrt(noisepsd))
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    ks= psd1d[0][0:-3]
    if _ADDGCS:
        xrange=[.03,110.]
    else:
        xrange= [0.,1.]
    yrange= [0.,10.]
    bovy_plot.bovy_print(fig_width=7.5)
    bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][0:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                        'ko',lw=2.,
                        zorder=12,
                        xlabel=r'$k\,(\mathrm{kpc}^{-1})$',
                        ylabel=r'$\sqrt{P_k}\,(\mathrm{km\,s}^{-1})$',
                        semilogx=_ADDGCS,
                        xrange=xrange,yrange=yrange)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][0:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*psd1d[2][0:-3]/numpy.sqrt(psd1d[1][0:-3]),
                    marker='None',ls='none',color='k')
    if _PLOTBAND:
        bovy_plot.bovy_plot(ks,
                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
                            '--',lw=2.,zorder=10,
                            color='0.85',overplot=True)
    if _PLOTBAND:
        bovy_plot.bovy_plot(ks,
                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
                            '-',lw=8.,zorder=9,
                            color='0.65',overplot=True)
    if _ADDGCS:
        ks_gcs, psd_gcs, e_psd_gcs= plot_psd_gcs()
    if _ADDRAVE:
        ks_rave, psd_rave, e_psd_rave= plot_psd_rave()
    if _ADDRED:
        plot_psd_red()
    interpks= list(ks[:-5])
    interppsd= list((scale*numpy.sqrt(psd1d[1][0:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)))[:-5])
    interppsd_w= (scale*0.5*psd1d[2][0:-3]/numpy.sqrt(psd1d[1][0:-3]))[:-5]
    interppsd_w[:8]*= 0.025 #fiddling to get a decent fit
    interppsd_w[3:5]*= 0.001
    interppsd_w= list(interppsd_w)
    interpks.append(0.025)
    interppsd.append(10.**-5.)
    interppsd_w.append(0.00001)
    if _ADDGCS:
        interpks.extend(ks_gcs)
        interppsd.extend(psd_gcs)
        interppsd_w.extend(e_psd_gcs)
    if _ADDRAVE:
        interpks.extend(ks_rave[5:])
        interppsd.extend(psd_rave[5:])
        interppsd_w.extend(e_psd_rave[5:])
    interpks= numpy.array(interpks)
    sortindx= numpy.argsort(interpks)
    interpks= interpks[sortindx]
    interppsd= numpy.array(interppsd)[sortindx]
    interppsd_w= interppsd/numpy.array(interppsd_w)[sortindx]
    interpindx= True-numpy.isnan(interppsd)
    #interpspec= interpolate.InterpolatedUnivariateSpline(interpks[interpindx],
    interpspec= interpolate.UnivariateSpline(numpy.log(interpks[interpindx]),
                                             numpy.log(interppsd[interpindx]/3.),
                                             w=interppsd_w,
                                             k=3,s=len(interppsd_w)*0.8)
    pks= numpy.linspace(interpks[0],interpks[-1],201)
    bovy_plot.bovy_plot(pks,
                        3.*numpy.exp(interpspec(numpy.log(pks))),
                        'k-',overplot=True)
    def my_formatter(x, pos):
        return r'$%g$' % x
    def my_formatter2(x, pos):
        return r'$%g$' % (1./x)
    major_formatter = FuncFormatter(my_formatter)
    major_formatter2 = FuncFormatter(my_formatter2)
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(major_formatter)
    ax2= pyplot.twiny()
    xmin, xmax= ax.xaxis.get_view_interval()
    ax2.set_xscale('log')
    ax2.xaxis.set_view_interval(xmin,xmax,ignore=True)
    ax2.set_xlabel('$\mathrm{Approximate\ scale}\,(\mathrm{kpc})$',
                   fontsize=12.,ha='center',x=0.5)
    ax2.xaxis.set_major_formatter(major_formatter2)
    bovy_plot.bovy_end_print('/Users/bovy/Desktop/test.png')
    return None

def plot_psd_gcs():
    data= hackGCS.hackGCS()
    dx= 0.025
    binsize= .8
    pix= pixelize_sample.pixelXY(data,xmin=-0.075,xmax=0.075,
                                 ymin=-0.075,ymax=0.075,
                                 dx=dx,dy=dx)
    resv= pix.plot('VVel',returnz=True,justcalc=True)
    resvunc= pix.plot('VVel',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    scale= 3.25/numpy.median(numpy.sqrt(noisepsd))
    print scale
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    ks= psd1d[0][0:-3]
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][0:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),'kx',mew=2.,
                        overplot=True)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][0:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*psd1d[2][0:-3]/numpy.sqrt(psd1d[1][0:-3]),
                    marker='None',ls='none',color='k')
    if False:
        interpindx= True-numpy.isnan(psd1d[1][0:-3])
        interpspec= interpolate.InterpolatedUnivariateSpline(ks[interpindx],
                                                             scale*numpy.sqrt(psd1d[1][0:-3])[interpindx],
                                                             k=3)
        pks= numpy.linspace(ks[0],ks[-1],201)
        bovy_plot.bovy_plot(pks,interpspec(pks),'k-',overplot=True)
    if _PLOTBAND:
        bovy_plot.bovy_plot(ks,
                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
                            '-',lw=8.,zorder=9,
                            color='0.65',overplot=True)
    return (ks,
            scale*numpy.sqrt(psd1d[1][0:-3]
                             -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
            scale*0.5*psd1d[2][0:-3]/numpy.sqrt(psd1d[1][0:-3]))

def plot_psd_rave():
    data= fitsio.read('/Users/bovy/data/rave/ravedr4_rc.fits')
    dx= .25
    binsize= 0.8#.735
    pix= pix= pixelize_sample.pixelXY(data,xmin=6.75,xmax=8.75,
                                      ymin=-1.75,ymax=0.25,
                                      dx=dx,dy=dx)
    resv= pix.plot(lambda x: dvlosgal(x,vtsun=230.),returnz=True,justcalc=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    scale= 3.25/numpy.median(numpy.sqrt(noisepsd))
    print scale
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    ks= psd1d[0][0:-3]
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][0:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),'k+',mew=2.,
                        overplot=True)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][0:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*psd1d[2][0:-3]/numpy.sqrt(psd1d[1][0:-3]),
                    marker='None',ls='none',color='k')
    if False:
        interpindx= True-numpy.isnan(psd1d[1][0:-3])
        interpspec= interpolate.InterpolatedUnivariateSpline(ks[interpindx],
                                                             scale*numpy.sqrt(psd1d[1][0:-3])[interpindx],
                                                             k=3)
        pks= numpy.linspace(ks[0],ks[-1],201)
        bovy_plot.bovy_plot(pks,interpspec(pks),'k-',overplot=True)
    if _PLOTBAND:
        bovy_plot.bovy_plot(ks,
                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
                            '-',lw=8.,zorder=9,
                            color='0.65',overplot=True)
        #bovy_plot.bovy_plot(numpy.tile(ks,(_NNOISE,1)).T,
        #                    scale*numpy.sqrt(noisepsd).T,
        #                    '-',zorder=0,alpha=0.5,
        #                    color='0.45',overplot=True)
    return (ks,
            scale*numpy.sqrt(psd1d[1][0:-3]
                             -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
            scale*0.5*psd1d[2][0:-3]/numpy.sqrt(psd1d[1][0:-3]))

def plot_psd_red():
    data= readAndHackHoltz.readAndHackHoltz()
    dx= 1.
    binsize= .735
    pix= pixelize_sample.pixelXY(data,xmin=5.,xmax=13,
                                 ymin=-3.,ymax=7.,
                                 dx=dx,dy=dx)
    resv= pix.plot(lambda x: dvlosgal(x),returnz=True,justcalc=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    scale= 3.25/numpy.median(numpy.sqrt(noisepsd))
    print scale
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    ks= psd1d[0][0:-3]
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][0:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),'ko',mew=2.,
                        mfc='none',overplot=True)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][0:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*psd1d[2][0:-3]/numpy.sqrt(psd1d[1][0:-3]),
                    marker='None',ls='none',color='k')
    if False:
        interpindx= True-numpy.isnan(psd1d[1][0:-3])
        interpspec= interpolate.InterpolatedUnivariateSpline(ks[interpindx],
                                                             scale*numpy.sqrt(psd1d[1][0:-3])[interpindx],
                                                             k=3)
        pks= numpy.linspace(ks[0],ks[-1],201)
        bovy_plot.bovy_plot(pks,interpspec(pks),'k-',overplot=True)
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    if _PLOTBAND:
        bovy_plot.bovy_plot(ks,
                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
                            '-',lw=8.,zorder=9,
                            color='0.65',overplot=True)
    return None

if __name__ == '__main__':
    plot_psd()
