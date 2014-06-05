import sys
import numpy
from scipy import optimize
import fitsio
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
import asymmetricDriftModel
import pixelize_sample
_EXT='png'
_ADDLLOGGCUT= True
_JACKERRS= True
_DEGTORAD= numpy.pi/180.
_VRSUN= -10.5
_VTSUN= 242.
_VZSUN= 7.25
def plot_2dkinematics(basesavefilename,datafilename=None):
    #Plot 2D los field
    #First read the sample
    if not datafilename is None:
        data= fitsio.read(datafilename)
    else:
        data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    print "Using %i stars for low-Z 2D kinematics analysis" % numpy.sum(indx)
    data= data[indx]
    pix= pixelize_sample.pixelXY(data)
    vmin, vmax= -75., 75.
    bovy_plot.bovy_print()
    pix.plot('VHELIO_AVG',
             zlabel=r'$\mathrm{median}\ V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
             vmin=vmin,vmax=vmax)
    bovy_plot.bovy_text(r'$\mathrm{typical\ uncertainty\!:}\ 3\,\mathrm{km\,s}^{-1}$',
                        bottom_left=True,size=18.)
    bovy_plot.bovy_text(r'$|Z| < 250\,\mathrm{pc}$',top_right=True,size=18.)
    bovy_plot.bovy_end_print(basesavefilename+'_XY.'+_EXT)
    #R,phi
    pix= pixelize_sample.pixelXY(data,rphi=True,
                                 ymin=-22.5,ymax=37.5,
#                                 dx=0.3,dy=2.)
#                                 dx=3.,dy=20.)
                                 dx=1.,dy=5.)
    bovy_plot.bovy_print()
    pix.plot('VHELIO_AVG',
             zlabel=r'$\mathrm{median}\ V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
             vmin=vmin,vmax=vmax)
    #bovy_plot.bovy_text(r'$|Z| < 250\,\mathrm{pc}$',bottom_left=True,size=18.)
    bovy_plot.bovy_end_print(basesavefilename+'_RPHI.'+_EXT)
    #Plot the dispersion / sqrt(n)
    pix= pixelize_sample.pixelXY(data,
                                 dx=1.,dy=1.)
    bovy_plot.bovy_print()
    pix.plot('VHELIO_AVG',
             func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
             zlabel=r'$\delta\ V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
             vmin=0.,vmax=10.)
    #bovy_plot.bovy_text(r'$|Z| < 250\,\mathrm{pc}$',bottom_left=True,size=18.)
    bovy_plot.bovy_end_print(basesavefilename+'_RPHI_vlosunc.'+_EXT)
    #Plot the dispersion
    bovy_plot.bovy_print()
    pix.plot('VHELIO_AVG',
             func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x))),
             zlabel=r'$\mathrm{MAD}\ V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
             vmin=0.,vmax=40.)
    #bovy_plot.bovy_text(r'$|Z| < 250\,\mathrm{pc}$',bottom_left=True,size=18.)
    bovy_plot.bovy_end_print(basesavefilename+'_RPHI_vlosdisp.'+_EXT)
    #Now plot the los velocity corrected for the Solar motion
    #XY
    vmin, vmax= -250., 250.
    bovy_plot.bovy_print()
    resv= pix.plot(lambda x: vlosgal(x),
                   zlabel=r'$\mathrm{median}\ V_{\mathrm{los,rot}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
    bovy_plot.bovy_end_print(basesavefilename+'_Vrot_XY.'+_EXT)
    pix= pixelize_sample.pixelXY(data,rphi=True,
                                 ymin=-22.5,ymax=37.5,
                                 dx=1.,dy=5.)
    #R,phi
    vmin, vmax= -250., 250.
    bovy_plot.bovy_print()
    resv= pix.plot(lambda x: vlosgal(x),
                   zlabel=r'$\mathrm{median}\ V_{\mathrm{los,rot}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
    #Plot tangent point
    rs= numpy.linspace(6.,8.,101)
    bovy_plot.bovy_plot(rs,numpy.arccos(rs/8.)*180./numpy.pi,'k--',lw=2.,
                        overplot=True)
    bovy_plot.bovy_plot(rs,-numpy.arccos(rs/8.)*180./numpy.pi,'k--',lw=2.,
                        overplot=True)
    bovy_plot.bovy_end_print(basesavefilename+'_Vrot_RPHI.'+_EXT)
    #Now plot the residuals wrt the Bovy et al. (2012) disk model
    #R,phi
    vmin, vmax= -20., 20.
    bovy_plot.bovy_print()
    resv= pix.plot(lambda x: dvlosgal(x,beta=0.,vc=218.),
                   zlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
    medindx= True-numpy.isnan(resv)
    meanresv= numpy.mean(resv[medindx])
    sigresv= numpy.std(resv[medindx])
    bovy_plot.bovy_text(r'$\mathrm{Residual} = %.1f \pm %.1f / %i\,\mathrm{km\,s}^{-1}$' % (meanresv,sigresv,round(numpy.sqrt(numpy.sum(medindx)))),
                        top_left=True,size=16.)
    bovy_plot.bovy_end_print(basesavefilename+'_dVBovy12_RPHI.'+_EXT)
    #Plot a histogram of the residuals
    bovy_plot.bovy_print()
    bovy_plot.bovy_hist(resv[medindx],color='k',bins=11,xrange=[-25.,25.],
                        yrange=[0.,15.],
                        histtype='step',normed=False,
                        xlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$')
    xs= numpy.linspace(-25.,25.,1001)
    bovy_plot.bovy_plot(xs,
                        50./11.*numpy.sum(medindx)/numpy.sqrt(2.*numpy.pi)/sigresv\
                            *numpy.exp(-0.5*(xs-meanresv)**2./sigresv**2.),
                        'k-',lw=2,overplot=True)
    bovy_plot.bovy_end_print(basesavefilename+'_dVBovy12_hist.'+_EXT) 
   #Now plot the residuals wrt the Bovy et al. (2012) disk model w/ Vc-240 and standard solar motion
    #R,phi
    vmin, vmax= -20., 20.
    bovy_plot.bovy_print()
    resv= pix.plot(lambda x: dvlosgal(x,beta=0.,vc=240.,vtsun=252.),
                   zlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
    medindx= True-numpy.isnan(resv)
    #bovy_plot.bovy_text(r'$\mathrm{Residual} = %.1f \pm %.1f / %i\,\mathrm{km\,s}^{-1}$' % (numpy.median(resv[medindx]),numpy.median(numpy.fabs(resv[medindx]-numpy.median(resv[medindx])))*1.4826,round(numpy.sqrt(numpy.sum(medindx)))),
#                        top_left=True,size=16.)
    bovy_plot.bovy_end_print(basesavefilename+'_dVBovy12Vc240VsSBD_RPHI.'+_EXT)
    #Now plot the residuals wrt the Bovy et al. (2012) disk model w/ Vc-220 and standard solar motion
    #R,phi
    vmin, vmax= -20., 20.
    bovy_plot.bovy_print()
    resv= pix.plot(lambda x: dvlosgal(x,beta=0.,vc=220.,vtsun=232.),
                   zlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
    medindx= True-numpy.isnan(resv)
#    bovy_plot.bovy_text(r'$\mathrm{Residual} = %.1f \pm %.1f / %i\,\mathrm{km\,s}^{-1}$' % (numpy.median(resv[medindx]),numpy.median(numpy.fabs(resv[medindx]-numpy.median(resv[medindx])))*1.4826,round(numpy.sqrt(numpy.sum(medindx)))),
#                        top_left=True,size=16.)
    bovy_plot.bovy_end_print(basesavefilename+'_dVBovy12Vc220VsSBD_RPHI.'+_EXT)
    #FFT
    pix= pixelize_sample.pixelXY(data,
                                 xmin=4.,xmax=14,
                                 ymin=-5.,ymax=5.)
    resv= pix.plot(lambda x: dvlosgal(x),
                   zlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      zlabel=r'$\delta\ V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                      vmin=0.,vmax=10.,returnz=True)
    pixsm= pixelize_sample.pixelXY(data,
                                   xmin=7.5,xmax=11.,
                                   ymin=-3.,ymax=3.,
                                   dx=0.35,dy=0.35)
    resvsm= pixsm.plot(lambda x: dvlosgal(x),
                   zlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
    resvsmunc= pixsm.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      zlabel=r'$\delta\ V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                      vmin=0.,vmax=10.,returnz=True)
    import psds
    psd2d= psds.power_spectrum(resv,
                               binsize=1.,radbins=7,wavenumber=True)
    psd2dsm= psds.power_spectrum(resvsm,
                                 binsize=1.,radbins=7,wavenumber=True)
    #Simulate the noise
    newresv= numpy.random.normal(size=resv.shape)*resvunc
    newresvsm= numpy.random.normal(size=resvsm.shape)*resvsmunc
    psd2dnoise= psds.power_spectrum(newresv,
                                    binsize=1.,radbins=7,wavenumber=True)
    psd2dsmnoise= psds.power_spectrum(newresvsm,
                                      binsize=1.,radbins=7,wavenumber=True)
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(psd2d[0][1:],numpy.sqrt(psd2d[1]/numpy.sqrt(numpy.prod(resv.shape)))[1:],'ko-',
                        xlabel=r'$k\,(\mathrm{kpc}^{-1})$',
                        xrange=[0.,.5],
                        yrange=[1.,30.],
                        semilogy=False)
    bovy_plot.bovy_plot(psd2dnoise[0][1:],numpy.sqrt(psd2dnoise[1]/numpy.sqrt(numpy.prod(resv.shape)))[1:],'-',color='0.5',overplot=True,lw=2.)
    #bovy_plot.bovy_plot(psd2dsm[0][1:]/0.35,numpy.sqrt(0.35*psd2dsm[1]/numpy.sqrt(numpy.prod(resvsm.shape)))[1:],'bo-',
    #                    overplot=True)
    #bovy_plot.bovy_plot(psd2dsmnoise[0][1:]/.35,numpy.sqrt(psd2dsmnoise[1]/numpy.sqrt(numpy.prod(resvsm.shape)))[1:],'-',color='0.5',overplot=True,lw=2.)
    bovy_plot.bovy_end_print(basesavefilename+'_FFTPSD.'+_EXT)
    return None

def dvlosgal(data,beta=0.,vc=218.,vtsun=_VTSUN):
    l= data['GLON']*_DEGTORAD
    sinl= numpy.sin(data['GLON']*_DEGTORAD)
    cosl= numpy.cos(data['GLON']*_DEGTORAD)
    sinb= numpy.sin(data['GLAT']*_DEGTORAD)
    cosb= numpy.cos(data['GLAT']*_DEGTORAD)
    vlosgal= data['VHELIO_AVG']/cosb-_VRSUN*cosl+vtsun*sinl\
        +_VZSUN*sinb/cosb\
        -(vc*(data['RC_GALR']/8.)**beta\
              -asymmetricDriftModel.va(data['RC_GALR']/8.,31.4/vc,
                                       vc=(data['RC_GALR']/8.)**beta,hR=3./8.,
                                       hs=33.3)*vc)*numpy.sin(l+data['RC_GALPHI'])
    return vlosgal

def vlosgal(data,vtsun=_VTSUN):
    sinl= numpy.sin(data['GLON']*_DEGTORAD)
    cosl= numpy.cos(data['GLON']*_DEGTORAD)
    sinb= numpy.sin(data['GLAT']*_DEGTORAD)
    cosb= numpy.cos(data['GLAT']*_DEGTORAD)
    vlosgal= data['VHELIO_AVG']/cosb-_VRSUN*cosl+vtsun*sinl\
        +_VZSUN*sinb/cosb
    return vlosgal

def linfit(x,slope,zeropoint):
    return slope*(x-8.)+zeropoint

if __name__ == '__main__':
    if len(sys.argv) > 2:
        plot_2dkinematics(sys.argv[1],sys.argv[2])
    else:
        plot_2dkinematics(sys.argv[1])
