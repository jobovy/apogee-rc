import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
def plot_magcorr(plotfilename):
    corrs= -1.61\
        -numpy.array([-1.615557,#solar
                       -1.653314,#avg'ed over Z dist
                       -1.652762,# exp. SFH
                       -1.652846])# Kroupa
    print corrs
    bovy_plot.bovy_print(fig_width=7.)
    bovy_plot.bovy_plot([1,2,3,4],corrs,'ko',ms=10.,
                        xrange=[0,5],
                        yrange=[-0.1,0.1],
                        ylabel=r'$\mathrm{magnitude\ offset}$')
    #Put labels and rotate them
    pyplot.xticks([1,2,3,4],
                  [r'$\mathrm{solar\ metallicity}$',
                   r"$\mathrm{avg'ed\ over\ metal.}$"+'\n'+r"$\mathrm{dist.}$",
                   r"$\mathrm{avg'ed\ over\ metal.}$"+'\n'+r"$\mathrm{dist.,\ exp.\ SFH}$",
                   r"$\mathrm{avg'ed\ over\ metal.}$"+'\n'+r"$\mathrm{dist.,\ Kroupa\,(2003)\ IMF}$"],size=16.,
                  rotation=45.)
    bovy_plot.bovy_end_print(plotfilename)

if __name__ == '__main__':
    plot_magcorr(sys.argv[1])
