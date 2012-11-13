import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
from galpy.util import bovy_plot
import rcmodel
def plot_vs_jkz(parser):
    options,args= parser.parse_args()
    if options.basti:
        zs= numpy.array([0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,
                         0.01,0.0198,0.03,0.04])
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
        zs= numpy.arange(0.0005,0.03005,0.005)
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')
        plotthis= pickle.load(savefile)
        jks= pickle.load(savefile)
        zs= pickle.load(savefile)
        savefile.close()
    else:
        njks= 31
        jks= numpy.linspace(0.5,0.75,njks)
        plotthis= numpy.zeros((njks,len(zs)))
        for ii in range(len(zs)):
            print zs[ii]
            rc= rcmodel.rcmodel(Z=zs[ii],loggmin=1.8,loggmax=2.8,
                                band=options.band,basti=options.basti,
                                imfmodel=options.imfmodel)
            for jj in range(njks):
                if options.type == 'mode':
                    try:
                        plotthis[jj,ii]= rc.mode(jks[jj])
                    except ValueError:
                        plotthis[jj,ii]= numpy.nan
                elif options.type == 'sig':
                    try:
                        plotthis[jj,ii]= rc.sigmafwhm(jks[jj])
                    except ValueError:
                        plotthis[jj,ii]= numpy.nan
        #Save
        savefile= open(args[0],'wb')
        pickle.dump(plotthis,savefile)
        pickle.dump(jks,savefile)
        pickle.dump(zs,savefile)
        savefile.close()
    #Plot
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='gray',
                          xrange=[jks[0],jks[-1]],
                          yrange=[zs[0],zs[-1]],
                          vmin=0.,vmax=0.8,
                          xlabel=r'$J-Ks$',
                          ylabel=r'$Z$',
                          interpolation='nearest',
                          colorbar=True,
                          shrink=0.78,
                          zlabel=r'$\mathrm{FWHM} / 2\sqrt{2\,\ln 2}\ [\mathrm{mag}]$')
    bovy_plot.bovy_end_print(options.outfilename)
    return None

def get_options():
    usage = "usage: %prog [options] savefilename"
    parser = OptionParser(usage=usage)
    #Data options
    parser.add_option("-o",dest='outfilename',default=None,
                      help="Name for plot file")
    parser.add_option("-t",dest='type',default='sig',
                      help="type of plot ('sig' for sigma, 'mode' for mode)")
    parser.add_option("-b",dest='band',default='Ks',
                      help="Band to use")
    parser.add_option("--imfmodel",dest='imfmodel',
                      default='lognormalChabrier2001',
                      help="IMF model to use")
    parser.add_option("--basti",action="store_true", dest="basti",
                      default=False,
                      help="If set, use BaSTI isochrones")
    return parser
    
if __name__ == '__main__':
    parser= get_options()
    plot_vs_jkz(parser)
