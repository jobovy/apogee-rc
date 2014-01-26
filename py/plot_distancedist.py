import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
_EXT='ps'
_ADDLLOGGCUT= True
def plot_distancedist(basesavefilename):
    #First read the sample
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Histogram
    bovy_plot.bovy_print(fig_width=6.)
    bovy_plot.bovy_hist(data['RC_DIST'],
                        histtype='step',
                        color='k',
                        lw=1.5,
                        bins=51,
                        xrange=[0.,10.],
                        xlabel=r'$\mathrm{distance}\,(\mathrm{kpc})$',
                        ylabel=r'$\mathrm{Number\ of\ stars}$',zorder=10)
    ax= pyplot.gca()
    ax2= ax.twinx()
    pyplot.sca(ax2)
    bovy_plot.bovy_hist(data['RC_DIST'],
                        histtype='step',
                        cumulative=True,
                        color='k',
                        lw=1.,normed=True,
                        overplot=True,bins=len(data))
    ax2.set_xlim(0.,10.)
    ax2.set_ylim(0.,1.)
    bovy_plot._add_ticks()
    ax2.set_ylabel(r'$\mathrm{cumulative\ distribution}$')
    bovy_plot.bovy_end_print(basesavefilename+'_hist.'+_EXT)
    #RZ
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(data['RC_GALR'],data['RC_GALZ'],
                        'k.',ms=2.,
                        xrange=[0.,16.],
                        yrange=[-3.,3.],
                        xlabel=r'$R\,(\mathrm{kpc})$',
                        ylabel=r'$Z\,(\mathrm{kpc})$',
                        onedhists=True)
    bovy_plot.bovy_end_print(basesavefilename+'_RZ.'+_EXT)
    #XY
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(data['RC_GALR']*numpy.cos(data['RC_GALPHI']),
                        data['RC_GALR']*numpy.sin(data['RC_GALPHI']),
                        'k.',
                        xrange=[0.,16.],
                        yrange=[5.,-5.],
                        xlabel=r'$X_{\mathrm{GC}}\,(\mathrm{kpc})$',
                        ylabel=r'$Y_{\mathrm{GC}}\,(\mathrm{kpc})$',
                        onedhists=True,ms=2.)
    bovy_plot.bovy_end_print(basesavefilename+'_XY.'+_EXT)

if __name__ == '__main__':
    plot_distancedist(sys.argv[1])
