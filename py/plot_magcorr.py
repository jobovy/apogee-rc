import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
def plot_magcorr(plotfilename):
    corrs= -1.61\
        -numpy.array([-1.668586,#solar
                       -1.649471,#avg'ed over Z dist
                       -1.643215,# exp. SFH
                       -1.649527])# Kroupa
    """
        -numpy.array([-1.615557,#solar
                       -1.653314,#avg'ed over Z dist
                       -1.652762,# exp. SFH
                       -1.652846])# Kroupa
    """
    print corrs
    bovy_plot.bovy_print(fig_width=7.,
                         text_fontsize=20.,
                         legend_fontsize=24.,
                         xtick_labelsize=18.,
                         ytick_labelsize=18.,
                         axes_labelsize=24.)
    bovy_plot.bovy_plot([1,2,3,4],corrs,'ko',ms=10.,
                        xrange=[0,5],
                        yrange=[0.,0.08],
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
