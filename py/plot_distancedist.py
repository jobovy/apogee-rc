import sys
import numpy
from scipy import optimize
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
import pixelize_sample
_EXT='ps'
def plot_distancedist(basesavefilename):
    #First read the sample
    data= apread.rcsample()
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
