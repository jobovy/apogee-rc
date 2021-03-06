import os, os.path
import cPickle as pickle
import numpy
from scipy import interpolate
from galpy.util import bovy_plot, multi
import isodist
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import rcmodel
from plot_vs_jkz import get_options
_CUTLOWAGE= True
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
        if options.type == 'mass' and len(args) == 3:
            #Also load mass_coarseage and omega
            savefile= open(args[1],'rb')
            masscoarse= pickle.load(savefile)
            savefile.close()
            savefile= open(args[2],'rb')
            omega= pickle.load(savefile)
            savefile.close()
    else:
        nages= 31
        if options.type == 'omega' or options.type == 'numfrac' \
                or options.coarseage:
            nages= 16
        lages= numpy.linspace(-1.,1.,nages)
        dlages= (lages[1]-lages[0])
        plotthis= numpy.zeros((len(zs),nages))
        multOut= multi.parallel_map(lambda x: _calc_one(zs[x],options,nages,lages,dlages),
                                    range(len(zs)),
                                    numcores=32)
        for ii in range(len(zs)):
            plotthis[ii,:]= multOut[ii]
        #Save
        savefile= open(args[0],'wb')
        pickle.dump(plotthis,savefile)
        pickle.dump(zs,savefile)
        pickle.dump(lages,savefile)
        savefile.close()
    #Plot
    #Fist cut out youngest ages, since they are irrelevant
    if _CUTLOWAGE:
        aindx= lages > numpy.log10(0.8)
        lages= lages[aindx]
        plotthis= plotthis[:,aindx]
    if options.type == 'mass':
        vmin, vmax= 0.5, 2.3
        vmin2, vmax2= 0.5, 2.
        zlabel= r'$\langle M_{\mathrm{RC}} \rangle \,(M_\odot)$'
        #cmap= 'gist_yarg'
        cmap= 'jet'
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
        #cmap= 'gist_yarg'
        cmap= 'jet'
        plotthis*= 100.
    elif options.type == 'numfrac':
        vmin, vmax= 0.,0.005
        vmin2, vmax2= 0.,0.004
        zlabel= r'$\mathrm{Number\ fraction\ in\ RC\ stars\ (\%)}$'
        #cmap= 'gist_yarg'
        cmap= 'jet'
        plotthis*= 100.
    elif options.type == 'massperrc':
        vmin, vmax= 0.,50000.
        vmin2, vmax2= 0.,25000.
        zlabel= r'$\mathrm{Stellar\ population\ mass\ per\ RC\ star}\,(M_\odot)$'
        #cmap= 'gist_yarg'
        cmap= 'jet'
        if options.redapogee:
            vmin, vmax= 0.,100000.
            vmin2, vmax2= 0.,200000.
            zlabel= r'$\mathrm{Mass\ fraction\ in}\ (J-K_s)_0 > 0.8\ \mathrm{giants\ (\%)}$'
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
    if options.remapz:
        zs= zs[:-1]
        plotthis= plotthis[:-1,:]
        fehs= numpy.linspace(-1.05,
                              isodist.Z2FEH(zs[-1],zsolar=0.017),len(zs))
        fehzs= isodist.FEH2Z(fehs,zsolar=0.017)
        new_plotthis= numpy.empty_like(plotthis)
        for ii in range(plotthis.shape[1]):
            goodz= True-numpy.isnan(plotthis[:,ii])
            tip= interpolate.InterpolatedUnivariateSpline(zs[goodz],
                                                          plotthis[goodz,ii],
                                                          k=3)
            new_plotthis[:,ii]= tip(fehzs)
            try:
                new_plotthis[fehs < numpy.nanmax(isodist.Z2FEH(zs[True-goodz],zsolar=0.017)),ii]= numpy.nan
            except ValueError: continue
        plotthis= new_plotthis
        xlabel= r'$[\mathrm{Fe/H}]\,(\mathrm{dex})$'
    else:
        xlabel= r'$Z$'
    bovy_plot.bovy_print(fig_height=7.,fig_width=6.)
    fig= pyplot.gcf()
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.6
    axBottom= pyplot.axes([left,bottom,width,height])
    fig.sca(axBottom)
    if options.remapz:
        xlimits= [fehs[0],fehs[-1]]
    else:
        xlimits= [zs[0],zs[-1]]
    ylimits= [lages[0]-dlages,lages[-1]+dlages]
    bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap=cmap,
                          xrange=xlimits,
                          yrange=ylimits,
                          vmin=vmin,vmax=vmax,
                          interpolation='nearest',
                          colorbar=True,
                          shrink=.9,
                          zlabel=zlabel,
                          overplot=True)
    extent= xlimits+ylimits
    pyplot.axis(extent)
    bovy_plot._add_axislabels(xlabel,
                              r'$\log_{10}\,\mathrm{Age} / 1\,\mathrm{Gyr}$')
    bovy_plot._add_ticks()
    left, bottom, width, height= 0.1, 0.68, 0.64, 0.2
    axTop= pyplot.axes([left,bottom,width,height])
    fig.sca(axTop)
    #Plot the average over SFH
    lages= numpy.linspace(-1.,1.,16)
    if _CUTLOWAGE:
        aindx= lages > numpy.log10(0.8)
        lages= lages[aindx]
        if options.type == 'mass':
            omega= omega[:,aindx]
            masscoarse= masscoarse[:,aindx]
    mtrend= numpy.zeros(len(zs))
    exppage= 10.**lages*numpy.exp((10.**(lages+2.))/800.) #e.g., Binney (2010)
    exexppage= 10.**lages*numpy.exp((10.**(lages+2.))/100.) #e.g., Binney (2010)
    page= 10.**lages
    if options.type == 'massperrc':
        mtrend= 1./(numpy.sum(page*1./plotthis,axis=1)/numpy.sum(page))
        expmtrend= 1./(numpy.sum(exppage*1./plotthis,axis=1)/numpy.sum(exppage))
        exexpmtrend= 1./(numpy.sum(exexppage*1./plotthis,axis=1)/numpy.sum(exexppage))
    elif options.type == 'mass' and len(args) == 3:
        if options.remapz:
            omega= omega[:-1,:]
            masscoarse= masscoarse[:-1,:]
        mtrend= numpy.nansum(page*omega,axis=1)/numpy.nansum(page*omega/masscoarse,axis=1)
        expmtrend= numpy.nansum(exppage*omega,axis=1)/numpy.nansum(exppage*omega/masscoarse,axis=1)
        exexpmtrend= numpy.nansum(exexppage*omega,axis=1)/numpy.nansum(exexppage*omega/masscoarse,axis=1)
    else:
        mtrend= numpy.sum(page*plotthis,axis=1)/numpy.sum(page)
        expmtrend= numpy.sum(exppage*plotthis,axis=1)/numpy.sum(exppage)
        exexpmtrend= numpy.sum(exexppage*plotthis,axis=1)/numpy.sum(exexppage)
    if options.remapz:
        zs= fehs
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

def _calc_one(z,options,nages,lages,dlages):
    print z
    if options.allapogee or options.redapogee:
        rc= rcmodel.rcmodel(Z=z,loggmax=3.5,
                            band=options.band,basti=options.basti,
                            imfmodel=options.imfmodel,
                            parsec=options.parsec,eta=options.eta)
    else:
        rc= rcmodel.rcmodel(Z=z,loggmin=1.8,loggmax='custom',
                            band=options.band,basti=options.basti,
                            imfmodel=options.imfmodel,
                            parsec=options.parsec,eta=options.eta)
    out= numpy.zeros(nages)
    for jj in range(nages):
        jk= rc._jks
        aindx= (rc._lages <= lages[jj]+dlages)\
            *(rc._lages > lages[jj]-dlages)
        if options.allapogee:
            aindx*= (jk > 0.5)
        elif options.redapogee:
            aindx*= (jk > 0.8)
        else:
            rcd= rcmodel.rcdist('../../rcdist-apogee/data/rcmodel_mode_jkz_ks_parsec_newlogg.sav')
            predH= numpy.array([rcd(j,z) for j in jk])
            predH= numpy.reshape(predH,len(jk))
            aindx*= (jk < 0.8)*(jk > 0.5)\
                *(z <= rcmodel.jkzcut(jk,upper=True))\
                *(z >= rcmodel.jkzcut(jk))\
                *(z <= 0.06)\
                *(rc._sample[:,1] > (predH-0.4))\
                *(rc._sample[:,1] < (predH+0.4))\
                *(rc._sample[:,1] > -3.)\
                *(rc._loggs[:,0] <= 3.5)
        if options.type == 'omega':
            try:
                out[jj]= numpy.mean(rc._massweights[aindx])
            except ValueError:
                out[jj]= numpy.nan
        elif options.type == 'numfrac':
            try:
                out[jj]= numpy.mean(rc._weights[aindx])
            except ValueError:
                out[jj]= numpy.nan
        elif options.type == 'mass':
            try:
                out[jj]= numpy.sum(rc._masses[aindx]*rc._weights[aindx])/numpy.sum(rc._weights[aindx])
            except ValueError:
                out[jj]= numpy.nan
    return out
                
if __name__ == '__main__':
    parser= get_options()
    astro_sampling(parser)
