import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
from galpy.util import bovy_plot, multi, save_pickles
import multiprocessing
from matplotlib import pyplot
import rcmodel
_NEW= True
if _NEW:
    loggmax= 'custom'
else:
    loggmax= 2.8
def plot_vs_jkz(parser):
    options,args= parser.parse_args()
    if options.basti:
        zs= numpy.array([0.004,0.008,0.01,0.0198,0.03,0.04])
    elif options.parsec:
        zs= numpy.arange(0.0005,0.06005,0.0005)
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')
        plotthis= pickle.load(savefile)
        jks= pickle.load(savefile)
        zs= pickle.load(savefile)
        savefile.close()
    else:
        njks= 101
        jks= numpy.linspace(0.5,0.8,njks)
        plotthis= numpy.zeros((njks,len(zs)))
        funcargs= (zs,options,njks,jks)
        multOut= multi.parallel_map((lambda x: indiv_calc(x,
                                                          *funcargs)),
                                    range(len(zs)),
                                    numcores=numpy.amin([64,len(zs),
                                                         multiprocessing.cpu_count()]))
        for ii in range(len(zs)):
            plotthis[:,ii]= multOut[ii]
        #Save
        save_pickles(args[0],plotthis,jks,zs)
    #Plot
    if options.type == 'sig':
        if options.band.lower() == 'age':
            if options.relative:
                raise NotImplementedError("relative age not implemented yet")
            else:
                vmin, vmax= 0.,.5
                zlabel= r'$\mathrm{FWHM} / 2\sqrt{2\,\ln 2}$'
        else:
            if options.relative:
                vmin, vmax= 0.8,1.2
                zlabel= r'$\mathrm{FWHM}/\mathrm{FWHM}_{\mathrm{fiducial}}$'
            else:
                vmin, vmax= 0., 0.4
                zlabel= r'$\mathrm{FWHM} / 2\sqrt{2\,\ln 2}$'
    elif options.type == 'mode':
        if options.band.lower() == 'age':
            if options.relative:
                raise NotImplementedError("relative age not implemented yet")
            else:
                vmin, vmax= 0.,1.
                zlabel= r'$\Delta\displaystyle\arg\!\max_{\substack{\log_{10}\mathrm{Age}}}{p(\log_{10}\mathrm{Age}|[J-K_s]_0)}$'
        else:
            if options.relative:
                vmin, vmax= -0.05,0.05
                zlabel= r'$\Delta\displaystyle\arg\!\max_{\substack{K_s}}{p(M_{K_s}|[J-K_s]_0)}$'
            else:
                vmin, vmax= -1.8, -1.5
                if options.band.lower() == 'h':
                    zlabel= r'$\displaystyle\arg\!\max_{\substack{H}}{p(M_{H}|[J-K_s]_0)}$'
                else:
                    zlabel= r'$\displaystyle\arg\!\max_{\substack{K_s}}{p(M_{K_s}|[J-K_s]_0)}$'
    if options.basti:#Remap the Zs
        zs= numpy.array([0.004,0.008,0.01,0.0198,0.03,0.04])
        regularzs= numpy.arange(0.0005,0.04005,0.0005)
        njks= len(jks)
        regularplotthis= numpy.zeros((njks,len(regularzs)))
        for jj in range(len(regularzs)):
            #Find z
            thisindx= numpy.argmin(numpy.fabs(regularzs[jj]-zs))
            for ii in range(njks):
                regularplotthis[ii,jj]= plotthis[ii,thisindx]
        zs= regularzs
        plotthis= regularplotthis
    if options.relative and os.path.exists(options.infilename):
        savefile= open(options.infilename,'rb')
        plotthisrel= pickle.load(savefile)
        savefile.close()
        if options.basti:
            plotthisrel= plotthisrel[:,:80]
        elif not options.parsec:
            plotthisrel= plotthisrel[:,:60]
        if options.type == 'mode':
            plotthis-= plotthisrel
        elif options.type == 'sig':
            plotthis/= plotthisrel
    bovy_plot.bovy_print()
    if options.type == 'sig':
        plotthis[numpy.isnan(plotthis)]= vmax
    if options.relative:
        #Only plot between the cuts
        for ii in range(plotthis.shape[0]):
            indx= zs >= rcmodel.jkzcut(jks[ii],upper=True)
            indx+= zs <= rcmodel.jkzcut(jks[ii],upper=False)
            plotthis[ii,indx]= numpy.nan
    bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='jet',
                          xrange=[jks[0],jks[-1]],
                          yrange=[zs[0],zs[-1]],
                          vmin=vmin,vmax=vmax,
                          xlabel=r'$(J-K_s)_0$',
                          ylabel=r'$Z$',
                          interpolation='nearest',
                          colorbar=True,
                          shrink=0.78,
                          zlabel=zlabel)
    #Overplot cuts
    bovy_plot.bovy_plot(jks,rcmodel.jkzcut(jks),
                        'w--',lw=2.,overplot=True)
    bovy_plot.bovy_plot(jks,rcmodel.jkzcut(jks,upper=True),
                        'w--',lw=2.,overplot=True)
    if options.basti:
        pyplot.annotate(r'$\mathrm{BaSTI}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif not options.parsec:
        pyplot.annotate(r'$\mathrm{Padova}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif options.imfmodel == 'kroupa2003':
        pyplot.annotate(r'$\mathrm{Kroupa\ (2003)\ IMF}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif 'expsfh' in args[0]:
        pyplot.annotate(r'$\mathrm{p(\mathrm{Age}) \propto e^{\mathrm{Age}/(8\,\mathrm{Gyr})}}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif not options.eta is None:
        pyplot.annotate(r'$\eta_{\mathrm{Reimers}} = %.1f$' % options.eta,
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif False:
        pyplot.annotate(r'$\mathrm{Padova}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    bovy_plot.bovy_end_print(options.outfilename)
    return None

def indiv_calc(ii,zs,options,njks,jks):
    print zs[ii]
    rc= rcmodel.rcmodel(Z=zs[ii],loggmin=1.8,loggmax=loggmax,
                        band=options.band,basti=options.basti,
                        imfmodel=options.imfmodel,
                        expsfh=options.expsfh,
                        parsec=options.parsec,
                        eta=options.eta)
    out= numpy.empty(njks)
    for jj in range(njks):
        if options.type == 'mode':
            try:
                out[jj]= rc.mode(jks[jj])
            except ValueError:
                out[jj]= numpy.nan
        elif options.type == 'sig':
            try:
                out[jj]= rc.sigmafwhm(jks[jj])
            except ValueError:
                out[jj]= numpy.nan
    return out

def get_options():
    usage = "usage: %prog [options] savefilename"
    parser = OptionParser(usage=usage)
    #Data options
    parser.add_option("-o",dest='outfilename',default=None,
                      help="Name for plot file")
    parser.add_option("-i",dest='infilename',default=None,
                      help="Name of the file that holds the surface that will be subtracted")
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
    parser.add_option("--expsfh",action="store_true", dest="expsfh",
                      default=False,
                      help="If set, use an exponentially-declining SFH")
    parser.add_option("--parsec",action="store_true", dest="parsec",
                      default=False,
                      help="If set, PARSEC")
    parser.add_option("-r","--relative",action="store_true", dest="relative",
                      default=False,
                      help="If set, plot relative surface")
    parser.add_option("--feh",dest='feh',default=0.,type='float',
                      help="Metallicity range to use")
    parser.add_option("--zsolar",dest='zsolar',default=0.017,type='float',
                      help="Solar Z")
    parser.add_option("--allapogee",action="store_true", dest="allapogee",
                      default=False,
                      help="If set, plot astro_sampling for all of APOGEE")
    parser.add_option("--redapogee",action="store_true", dest="redapogee",
                      default=False,
                      help="If set, plot astro_sampling for APOGEE J-Ks > 0.8")
    parser.add_option("--coarseage",action="store_true", dest="coarseage",
                      default=False,
                      help="If set, use a coarse age distribution for the average RC mass")
    parser.add_option("--remapz",action="store_true", dest="remapz",
                      default=False,
                      help="If set, remap Z onto [Fe/H] for astro_sampling")
    parser.add_option("--eta",dest='eta',default=None,type='float',
                      help="Mass-loss efficiency parameter")
    return parser
    
if __name__ == '__main__':
    parser= get_options()
    plot_vs_jkz(parser)
