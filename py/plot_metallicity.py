import sys
import numpy
from scipy import optimize
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
import pixelize_sample
_EXT='ps'
def plot_metallicity(basesavefilename):
    #First read the sample
    data= apread.rcsample()
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.05)*(data['METALS'] > -1000.)
    data= data[indx]
    #First do the metallicity gradient
    rs= numpy.arange(0.5,18.5001,.85)
    fehs= numpy.zeros_like(rs)+numpy.nan
    sigfehs= numpy.zeros_like(rs)+numpy.nan
    efehs= numpy.zeros_like(rs)+numpy.nan
    for ii in range(len(rs)-1):
        tindx= (data['RC_GALR'] > rs[ii])\
            *(data['RC_GALR'] <= rs[ii+1])
        if numpy.sum(tindx) < 20: continue
        fehs[ii]= numpy.median(data['METALS'][tindx])
        sigfehs[ii]= numpy.std(data['METALS'][tindx])
        efehs[ii]= numpy.std(data['METALS'][tindx])/numpy.sqrt(numpy.sum(tindx))
    #Plot
    bovy_plot.bovy_print(fig_width=6.)
    bovy_plot.bovy_plot(rs+0.5*(rs[1]-rs[0]),fehs,'ko',
                        xlabel=r'$R\,(\mathrm{kpc})$',
                        ylabel=r'$\mathrm{median\ [Fe/H]}\,(\mathrm{dex})$',
                        xrange=[0.,15.],
                        yrange=[-.6,0.5],zorder=10)
    pyplot.errorbar(rs+0.5*(rs[1]-rs[0]),fehs,yerr=efehs,marker='None',
                    ls='none',color='k',zorder=9)
    #bovy_plot.bovy_plot(data['RC_GALR'],data['METALS'],'k,',overplot=True)
    #FIT
    indx= True-numpy.isnan(fehs)
    bestfit= optimize.curve_fit(linfit,(rs+0.5*(rs[1]-rs[0]))[indx],
                                fehs[indx],sigma=efehs[indx],
                                p0=(-0.1,0.))
    bovy_plot.bovy_plot(rs+0.5*(rs[1]-rs[0]),
                        bestfit[0][0]*((rs+0.5*(rs[1]-rs[0])-8.))+bestfit[0][1],
                        'k--',overplot=True)
    bovy_plot.bovy_text(r'$|Z| < 50\,\mathrm{pc}$',top_right=True,size=18.)
    bovy_plot.bovy_text(r'$[\mathrm{Fe/H}] = %.2f\,(R-8\,\mathrm{kpc}) + %.2f$' % (bestfit[0][0],bestfit[0][1]),
                        bottom_left=True,size=18.)
    bovy_plot.bovy_end_print(basesavefilename+'_radialgradient.'+_EXT)
    #Now plot azimuthal stuff
    #First read the sample again
    data= apread.rcsample()
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    data= data[indx]
    pix= pixelize_sample.pixelXY(data)
    bovy_plot.bovy_print()
    pix.plot('METALS',
             zlabel=r'$\mathrm{median\ [Fe/H]}\,(\mathrm{dex})$',
             vmin=-0.4,vmax=0.3)
    bovy_plot.bovy_text(r'$\mathrm{typical\ uncertainty\!:}\ 0.02\,\mathrm{dex}$',
                        bottom_left=True,size=18.)
    bovy_plot.bovy_text(r'$|Z| < 250\,\mathrm{pc}$',top_right=True,size=18.)
    bovy_plot.bovy_end_print(basesavefilename+'_XY.'+_EXT)
    #R,phi
    pix= pixelize_sample.pixelXY(data,rphi=True,
                                 ymin=-22.5,ymax=37.5,dy=5.)
    bovy_plot.bovy_print()
    pix.plot('METALS',
             zlabel=r'$\delta\,\mathrm{median\ [Fe/H]}\,(\mathrm{dex})$',
             vmin=-0.1,vmax=0.1,submediany=True)
    #bovy_plot.bovy_text(r'$|Z| < 250\,\mathrm{pc}$',bottom_left=True,size=18.)
    bovy_plot.bovy_end_print(basesavefilename+'_RPHI.'+_EXT)

def linfit(x,slope,zeropoint):
    return slope*(x-8.)+zeropoint

if __name__ == '__main__':
    plot_metallicity(sys.argv[1])
