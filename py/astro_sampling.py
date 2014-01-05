import os, os.path
import cPickle as pickle
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import rcmodel
from plot_vs_jkz import get_options
def astro_sampling(parser):
    options,args= parser.parse_args()
    if options.basti:
        zs= numpy.array([0.004,0.008,0.01,0.0198,0.03,0.04])
    elif options.parsec:
        zs= numpy.arange(0.0005,0.06005,0.0005)
        #zs= numpy.arange(0.0005,0.06005,0.005)
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
        #zs= numpy.arange(0.0005,0.03005,0.005)
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')
        plotthis= pickle.load(savefile)
        zs= pickle.load(savefile)
        lages= pickle.load(savefile)
        savefile.close()
        dlages= (lages[1]-lages[0])
        if options.type == 'massperrc':
            savefile= open(args[1],'rb')
            plotthis/= pickle.load(savefile)
            savefile.close()
    else:
        nages= 31
        if options.type == 'omega' or options.type == 'numfrac' \
                or options.coarseage:
            nages= 16
        lages= numpy.linspace(-1.,1.,nages)
        dlages= (lages[1]-lages[0])
        plotthis= numpy.zeros((len(zs),nages))
        for ii in range(len(zs)):
            print zs[ii]
            if options.allapogee or options.redapogee:
                rc= rcmodel.rcmodel(Z=zs[ii],loggmax=3.5,
                                    band=options.band,basti=options.basti,
                                    imfmodel=options.imfmodel,
                                    parsec=options.parsec)
            else:
                rc= rcmodel.rcmodel(Z=zs[ii],loggmin=1.8,loggmax=2.8,
                                    band=options.band,basti=options.basti,
                                    imfmodel=options.imfmodel,
                                    parsec=options.parsec)
            for jj in range(nages):
                jk= rc._jks
                aindx= (rc._lages <= lages[jj]+dlages)\
                    *(rc._lages > lages[jj]-dlages)
                if options.allapogee:
                    aindx*= (jk > 0.5)
                elif options.redapogee:
                    aindx*= (jk > 0.8)
                else:
                    aindx*= (jk < 0.8)*(jk > 0.5)\
                        *(zs[ii] <= rcmodel.jkzcut(jk,upper=True))\
                        *(zs[ii] >= rcmodel.jkzcut(jk))\
                        *(zs[ii] <= 0.06)
                if options.type == 'omega':
                    try:
                        plotthis[ii,jj]= numpy.mean(rc._massweights[aindx])
                    except ValueError:
                        plotthis[ii,jj]= numpy.nan
                elif options.type == 'numfrac':
                    try:
                        plotthis[ii,jj]= numpy.mean(rc._weights[aindx])
                    except ValueError:
                        plotthis[ii,jj]= numpy.nan
                elif options.type == 'mass':
                    try:
                        plotthis[ii,jj]= numpy.sum(rc._masses[aindx]*rc._weights[aindx])/numpy.sum(rc._weights[aindx])
                    except ValueError:
                        plotthis[ii,jj]= numpy.nan
        #Save
        savefile= open(args[0],'wb')
        pickle.dump(plotthis,savefile)
        pickle.dump(zs,savefile)
        pickle.dump(lages,savefile)
        savefile.close()
    #Plot
    #Fist cut out youngest ages, since they are irrelevant
    aindx= lages > numpy.log10(0.8)
    lages= lages[aindx]
    plotthis= plotthis[:,aindx]
    if options.type == 'mass':
        vmin, vmax= 0.5, 2.3
        vmin2, vmax2= 0.5, 2.
        zlabel= r'$\langle M_{\mathrm{RC}} \rangle \,(M_\odot)$'
        cmap= 'gist_yarg'
    elif options.type == 'omega':
        vmin, vmax= 0.,.03
        vmin2, vmax2= 0.,.015
        if options.allapogee:
            vmin, vmax= 0.,.035
            zlabel= r'$\mathrm{Mass\ fraction\ in}\ (J-K_s)_0 > 0.5\ \mathrm{giants\ (\%)}$'
        elif options.redapogee:
            vmin, vmax= 0.,.005
            vmin2, vmax2= 0.,.003
            zlabel= r'$\mathrm{Mass\ fraction\ in}\ (J-K_s)_0 > 0.8\ \mathrm{giants\ (\%)}$'
        else:
            zlabel= r'$\mathrm{Mass\ fraction\ in\ RC\ stars\ (\%)}$'
        cmap= 'gist_yarg'
        plotthis*= 100.
    elif options.type == 'numfrac':
        vmin, vmax= 0.,0.005
        vmin2, vmax2= 0.,0.004
        zlabel= r'$\mathrm{Number\ fraction\ in\ RC\ stars\ (\%)}$'
        cmap= 'gist_yarg'
        plotthis*= 100.
    elif options.type == 'massperrc':
        vmin, vmax= 0.,100.
        vmin2, vmax2= 0.,100.
        zlabel= r'$\mathrm{Stellar\ population\ mass\ per\ RC\ star}$'
        cmap= 'gist_yarg'
    print numpy.nanmin(plotthis), numpy.nanmax(plotthis)
    if options.basti:#Remap the Zs
        zs= numpy.array([0.004,0.008,0.01,0.0198,0.03,0.04])
        regularzs= numpy.arange(0.0005,0.03005,0.0005)/0.019*0.0198
        regularplotthis= numpy.zeros((nages,len(regularzs)))
        for jj in range(len(regularzs)):
            #Find z
            thisindx= numpy.argmin(numpy.fabs(regularzs[jj]-zs))
            for ii in range(nages):
                regularplotthis[ii,jj]= plotthis[ii,thisindx]
        zs= regularzs
        plotthis= regularplotthis
    bovy_plot.bovy_print(fig_height=7.,fig_width=6.)
    fig= pyplot.gcf()
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.6
    axBottom= pyplot.axes([left,bottom,width,height])
    fig.sca(axBottom)
    xlimits= [zs[0],zs[-1]]
    ylimits= [lages[0]-dlages,lages[-1]+dlages]
    bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap=cmap,
                          xrange=xlimits,
                          yrange=ylimits,
                          vmin=vmin,vmax=vmax,
                          interpolation='nearest',
                          colorbar=True,
                          shrink=0.7,
                          zlabel=zlabel,
                          overplot=True)
    extent= xlimits+ylimits
    pyplot.axis(extent)
    bovy_plot._add_axislabels(r'$Z$',
                              r'$\log_{10}\,\mathrm{Age} / 1\,\mathrm{Gyr}$')
    bovy_plot._add_ticks()
    left, bottom, width, height= 0.1, 0.68, 0.64, 0.2
    axTop= pyplot.axes([left,bottom,width,height])
    fig.sca(axTop)
    #Plot the average over SFH
    mtrend= numpy.zeros(len(zs))
    exppage= 10.**lages*numpy.exp((10.**(lages+2.))/800.) #e.g., Binney (2010)
    exexppage= 10.**lages*numpy.exp((10.**(lages+2.))/100.) #e.g., Binney (2010)
    page= 10.**lages
    mtrend= numpy.sum(page*plotthis,axis=1)/numpy.sum(page)
    expmtrend= numpy.sum(exppage*plotthis,axis=1)/numpy.sum(exppage)
    exexpmtrend= numpy.sum(exexppage*plotthis,axis=1)/numpy.sum(exexppage)
    pyplot.plot(zs,mtrend,'k-')
    pyplot.plot(zs,expmtrend,'k--')
    pyplot.plot(zs,exexpmtrend,'k-.')
    pyplot.ylim(vmin2,vmax2)
    pyplot.xlim(xlimits[0],xlimits[1])
    nullfmt   = NullFormatter()         # no labels
    thisax= pyplot.gca()
    thisax.xaxis.set_major_formatter(nullfmt)
    bovy_plot._add_ticks()
    if options.type == 'mass':
        pyplot.ylabel(zlabel)
    if options.basti:
        pyplot.annotate(r'$\mathrm{BaSTI}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif options.imfmodel == 'kroupa2003':
        pyplot.annotate(r'$\mathrm{Padova, Kroupa\ (2003)\ IMF}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif 'expsfh' in args[0]:
        pyplot.annotate(r'$\mathrm{Padova, p(\mathrm{Age}) \propto e^{\mathrm{Age}/(8\,\mathrm{Gyr})}}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    elif options.parsec:
        pass
        #pyplot.annotate(r'$\mathrm{PARSEC}$',
        #                (0.5,1.08),xycoords='axes fraction',
        #                horizontalalignment='center',
        #                verticalalignment='top',size=16.)
    else:
        pyplot.annotate(r'$\mathrm{Padova}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=16.)
    bovy_plot.bovy_end_print(options.outfilename)
    return None

if __name__ == '__main__':
    parser= get_options()
    astro_sampling(parser)
