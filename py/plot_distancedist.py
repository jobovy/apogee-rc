import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
_EXT='png'
_ADDLLOGGCUT= False
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
    bovy_plot.bovy_plot(sorted(data['RC_DIST']),
                        numpy.linspace(1./len(data),1.,len(data)),
                        'k-',
                        overplot=True)
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
    xarr, dx=1.5, -1.
    from matplotlib.patches import FancyArrowPatch
    _legendsize= 16
    arr= FancyArrowPatch(posA=(xarr+0.05,0.),
                         posB=(xarr+dx*10./8.,0.),
                         arrowstyle='->', 
                         connectionstyle='arc3,rad=%4.2f' % (0.), 
                         shrinkA=2.0, shrinkB=2.0,
                         mutation_scale=20.0, 
                         mutation_aspect=None,fc='k')
    ax = pyplot.gca()
    ax.add_patch(arr)
    bovy_plot.bovy_text(xarr+7.*dx/8.,-0.45,r'$\mathrm{GC}$',
                        size=_legendsize)
    arr= FancyArrowPatch(posA=(1.5,-0.05),
                         posB=(1.5,.75),
                         arrowstyle='->', 
                         connectionstyle='arc3,rad=%4.2f' % (0.), 
                         shrinkA=2.0, shrinkB=2.0,
                         mutation_scale=20.0, 
                                 mutation_aspect=None,fc='k')
    ax = pyplot.gca()
    ax.add_patch(arr)
    bovy_plot.bovy_text(1.59,0.2,r'$\mathrm{NGP}$',
                        size=_legendsize)
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
    xarr, dx= 2.2, -2.
    arr= FancyArrowPatch(posA=(xarr,0.),
                         posB=(xarr+dx,0.),
                         arrowstyle='->', 
                         connectionstyle='arc3,rad=%4.2f' % (0.), 
                         shrinkA=2.0, shrinkB=2.0,
                         mutation_scale=20.0, 
                         mutation_aspect=None,fc='k')
    ax = pyplot.gca()
    ax.add_patch(arr)
    bovy_plot.bovy_text(xarr+7.*dx/8.,-0.25,r'$\mathrm{GC}$',
                        size=_legendsize)
    xcen, ycen, dr, t= -2., 0., 4., 14.*numpy.pi/180.
    arr= FancyArrowPatch(posA=(xcen+dr*numpy.cos(t),
                               ycen+dr*numpy.sin(t)),
                         posB=(xcen+dr*numpy.cos(-t),
                               ycen+dr*numpy.sin(-t)),
                         arrowstyle='<-', 
                         connectionstyle='arc3,rad=%4.2f' % (2.*t), 
                         shrinkA=2.0, shrinkB=2.0,
                         mutation_scale=20.0, 
                         mutation_aspect=None,fc='k')
    ax.add_patch(arr)
    bovy_plot.bovy_end_print(basesavefilename+'_XY.'+_EXT)

if __name__ == '__main__':
    plot_distancedist(sys.argv[1])
