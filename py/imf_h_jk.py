#Plot the IMF vs. H and J-K, sensibly, also plot logg etc.
import sys
import os, os.path
import copy
import math
import numpy
from optparse import OptionParser
import isodist, isodist.imf
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter, MultipleLocator
from matplotlib.patches import FancyArrowPatch
OUTDIR= os.path.join(os.getenv('HOME'),'Desktop')
#OUTDIR= '../tex/'
OUTEXT= 'png'
_ADDEXTINCT= False
_ADDCOLORCUT= True
def imf_h_jk(plotfile,Z=None,dwarf=False,log=False,h=12.,basti=False,
             dartmouth=False,kroupa=False,loggash=False):
    #Read isochrones
    if basti:
        zs= numpy.array([0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,
                         0.01,0.0198,0.03,0.04])
    elif dartmouth:
        zs= isodist.FEH2Z(numpy.array([-2.5,-2.,-1.5,-1.,-0.5,0.,0.2,0.3,0.5]))
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
    if Z is None:
        Zs= zs
    elif not basti and not dartmouth:
        if Z < 0.0015:
            Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
        elif Z < 0.01:
            Zs= [Z-0.001,Z-0.0005,Z,Z+0.0005,Z+0.001] #build up statistics
        else:
            Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
    else:
        Zs= [Z]
    if basti:
        p= isodist.BastiIsochrone(Z=Zs)
    elif dartmouth:
        p= isodist.DartmouthIsochrone(feh=isodist.Z2FEH(Zs),onlyold=True)
    else:
        p= isodist.PadovaIsochrone(Z=Zs)
    #Get relevant data
    sample= []
    weights= []
    quants= []
    for logage in p.logages():
        for z in Zs:
            thisiso= p(logage,z)
            if basti: mini= thisiso['M_ini']
            elif dartmouth: mini= thisiso['M']
            else: mini= thisiso['M_ini']
            if basti:
                int_IMF= isodist.imf.lognormalChabrier2001(thisiso['M_ini'],
                                                           int=True)
                dmpm= numpy.roll(int_IMF,-1)-int_IMF
            elif dartmouth:
                int_IMF= isodist.imf.lognormalChabrier2001(thisiso['M'],
                                                           int=True)
                dmpm= numpy.roll(int_IMF,-1)-int_IMF
            else:
                if kroupa:
                    int_IMF= isodist.imf.kroupa2003(thisiso['M_ini'],int=True)
                    dmpm= numpy.roll(int_IMF,-1)-int_IMF
                else:
                    dmpm= numpy.roll(thisiso['int_IMF'],-1)-thisiso['int_IMF']
            for ii in range(1,len(mini)-1):
                if basti:
                    JK= thisiso['J'][ii]-thisiso['K'][ii]
                    H= thisiso['H'][ii]
                else:
                    JK= thisiso['J'][ii]-thisiso['Ks'][ii]
                    H= thisiso['H'][ii]
                if JK < 0.: # or thisiso['logg'][ii] > 3.5:
                    continue
                if dmpm[ii] > 0.: 
                    if basti:
                        if loggash:
                            sample.append([thisiso['J'][ii]-thisiso['K'][ii],
                                           thisiso['logg'][ii]])
                        else:
                            sample.append([thisiso['J'][ii]-thisiso['K'][ii],
                                           thisiso['H'][ii]])
                    else:
                        if loggash:
                            sample.append([thisiso['J'][ii]-thisiso['Ks'][ii],
                                           thisiso['logg'][ii]])
                        else:
                            sample.append([thisiso['J'][ii]-thisiso['Ks'][ii],
                                           thisiso['H'][ii]])
                    if dartmouth:
                        if logage > numpy.log10(5.)+9.:
                            weights.append(2.*dmpm[ii]) #Dartmouth are linearly spaced, but spacing is bigger at > 5 Gyr
                        else:
                            weights.append(dmpm[ii]) #Dartmouth are linearly spaced?
                    else:
                        weights.append(dmpm[ii]*10**(logage-7.))
                    #weights.append(dmpm[ii]*10**(logage-7.)*numpy.exp((10.**(logage-7.))/800.))
                    if 'logg' in options.type:
                        quants.append(thisiso['logg'][ii])
                    elif 'teff' in options.type:
                        quants.append(10.**(thisiso['logTe'][ii]))
                    elif 'M_ini' in options.type:
                        quants.append(thisiso['M_ini'][ii])
                    elif 'logage' in options.type:
                        quants.append(thisiso['logage'][ii]-9.)
                else: 
                    continue #no use in continuing here
    #Form array
    sample= numpy.array(sample)
    weights= numpy.array(weights)
    quants= numpy.array(quants)
    #Histogram
    if loggash: yrange=[1.5,3.4]
    elif dwarf: yrange= [2.,9.]
    else: yrange= [-11.,2.]
    if dwarf:
        hist, edges= numpy.histogramdd(sample,weights=weights,bins=51,
                                       range=[[0.,1.6],yrange])
        if not options.type == 'dens':
            quanthist, edges= numpy.histogramdd(sample,weights=weights*quants,
                                                bins=51,
                                                range=[[0.,1.6],yrange])
            if 'sigma' in options.type:
                quant2hist, edges= numpy.histogramdd(sample,weights=weights*quants**2.,
                                                     bins=51,
                                                     range=[[0.,1.6],yrange])
    else:
        hist, edges= numpy.histogramdd(sample,weights=weights,bins=49,
#                                       range=[[0.5,0.8],yrange])
                                       range=[[0.3,1.6],yrange])
        if not options.type == 'dens':
            quanthist, edges= numpy.histogramdd(sample,weights=weights*quants,
                                                bins=49,
                                                range=[[0.3,1.6],yrange])
            if 'sigma' in options.type:
                quant2hist, edges= numpy.histogramdd(sample,weights=weights*quants**2.,
                                                     bins=49,
                                                     range=[[0.3,1.6],yrange])
    if not options.type == 'dens':
        #Determine normalized hist
        normhist= numpy.zeros_like(hist)
        for ii in range(len(hist[:,0])):
            normhist[ii,:]= hist[ii,:]/numpy.nanmax(hist[ii,:])*numpy.nanmax(hist)
        indx= (normhist > 10.**-2.*numpy.nanmax(normhist))
        if 'sigma' in options.type:
            hist[indx]= numpy.sqrt((quant2hist[indx]/hist[indx]-(quanthist[indx]/hist[indx])**2.))
        else:
            hist[indx]= quanthist[indx]/hist[indx]
        hist[True-indx]= numpy.nan
        for ii in range(len(hist[:,0])):
            rev= copy.copy(hist[ii,::-1]) #reverse, but in one go does not always work
            hist[ii,:]= rev
    else:
        #Normalize each J-K
        for ii in range(len(hist[:,0])):
            hist[ii,:]/= numpy.nanmax(hist[ii,:])/numpy.nanmax(hist)
            rev= copy.copy(hist[ii,::-1]) #reverse, but in one go does not always work
            hist[ii,:]= rev
    #Plot
    bovy_plot.bovy_print()
    if log:
        hist= numpy.log(hist)
    if options.type == 'loggmean':
        colorbar= True
        zlabel= r'$\mathrm{average}\ \log g$'
        cmap= 'jet'
        zrange= [0.,4.5]
    elif options.type == 'loggsigma':
        colorbar= True
        zlabel= r'$\sigma_{\log g}$'
        cmap= 'jet'
        zrange= [0.,0.1]
    elif options.type == 'teffmean':
        colorbar= True
        zlabel= r'$\mathrm{average}\ T_{\mathrm{eff}}\ [\mathrm{K}]$'
        cmap= 'jet'
        zrange= [3000.,6000.]
    elif options.type == 'teffsigma':
        colorbar= True
        zlabel= r'$\sigma_{T_{\mathrm{eff}}}\ [\mathrm{K}]$'
        cmap= 'jet'
        zrange= [0.,100.]
    elif options.type == 'M_inimean':
        colorbar= True
        zlabel= r'$\mathrm{average}\ M_{\mathrm{init}}\ [\mathrm{M_\odot}]$'
        cmap= 'jet'
        zrange= [0.,4.]
    elif options.type == 'M_inisigma':
        colorbar= True
        zlabel= r'$\sigma_{M_{\mathrm{init}}}\ [\mathrm{M_\odot}]$'
        cmap= 'jet'
        zrange= [0.,0.2]
    elif options.type == 'logagemean':
        colorbar= True
        zlabel= r'$\mathrm{average}\ \log_{10} \mathrm{Age} / \mathrm{1\,Gyr}$'
        cmap= 'jet'
        zrange= [-1.,1.]
    elif options.type == 'logagesigma':
        colorbar= True
        zlabel= r'$\sigma_{\log_{10} \mathrm{Age} / \mathrm{1\,Gyr}}$'
        cmap= 'jet'
        zrange= [0.,.2]
    else:
        colorbar= True
        zlabel= None
        cmap= 'gist_yarg'
        zrange= [numpy.nanmin(hist),numpy.nanmax(hist)]
    shrink=.78
    if loggash:
        ylabel= r'$\log g$'
    else:
        ylabel=r'$M_H\ [\mathrm{mag}]$'
    bovy_plot.bovy_dens2d(hist.T,origin='lower',cmap=cmap,
                          xrange=[edges[0][0],edges[0][-1]],
                          yrange=[edges[1][-1],edges[1][0]],
                          aspect=(edges[0][-1]-edges[0][0])/float(edges[1][-1]-edges[1][0]),
                          xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                          ylabel=ylabel,
                          colorbar=colorbar,zlabel=zlabel,
                          shrink=shrink,vmin=zrange[0],vmax=zrange[1],
                          interpolation='nearest')
    if _ADDEXTINCT and not loggash:
        #Add extinction arrow
        djk= 0.45
        dh= 1.55/1.5*djk
        ax=pyplot.gca()
        ax.add_patch(FancyArrowPatch((1.,-2.),(1+djk,-2+dh),
                                     arrowstyle='->',mutation_scale=20,fill=True,
                                     lw=1.25))
        bovy_plot.bovy_text(1.05,-2.05,r'$\mathrm{extinction}$',
                            rotation=-math.atan(1.5/1.55*1.3/13.)/math.pi*180.,
                            size=14.)
    if _ADDCOLORCUT and not loggash:
        ax=pyplot.gca()
        #Add color cut
        bovy_plot.bovy_plot([0.5,0.5],[-20.,20.],'--',color='0.6',overplot=True)
        ax.add_patch(FancyArrowPatch((0.5,-6.),(0.7,-6.),
                                     arrowstyle='->',mutation_scale=20,fill=True,
                                     lw=1.25,ls='dashed',color='0.6'))
        bovy_plot.bovy_text(0.43,-8.,r'$\mathrm{APOGEE\ color\ cut}$',rotation=90.,
                            size=14.)
    if _ADDEXTINCT:
        #Add twin y axis
        ax= pyplot.gca()
        def my_formatter(x, pos):
            """distance in kpc for m=h"""
            xs= 10.**((h-x)/5.-2.)
            return r'$%.0f$' % xs
        def my_formatter2(x, pos):
            """distance in kpc for m=h"""
            xs= 10.**((h-x)/5.+1.)
            return r'$%.0f$' % xs
        ax2= pyplot.twinx()
        if dwarf:
            major_formatter = FuncFormatter(my_formatter2)
        else:
            major_formatter = FuncFormatter(my_formatter)
        ax2.yaxis.set_major_formatter(major_formatter)
        ystep= ax.yaxis.get_majorticklocs()
        ystep= ystep[1]-ystep[0]
        ax2.yaxis.set_minor_locator(MultipleLocator(ystep/5.))
        ax2.yaxis.tick_right()
        ax2.yaxis.set_label_position('right')
        ymin, ymax= ax.yaxis.get_view_interval()
        ax2.yaxis.set_view_interval(ymin,ymax,ignore=True)
        if dwarf:
            ax2.set_ylabel('$\mathrm{distance\ for}\ H_0\ =\ %.1f\ [\mathrm{pc}]$' % h)
        else:
            ax2.set_ylabel('$\mathrm{distance\ for}\ H_0\ =\ %.1f\ [\mathrm{kpc}]$' % h)
        xstep= ax.xaxis.get_majorticklocs()
        xstep= xstep[1]-xstep[0]
        ax2.xaxis.set_minor_locator(MultipleLocator(xstep/5.))
    if Z is None or not options.type == 'dens':
        bovy_plot.bovy_end_print(plotfile)
    else:
#        bovy_plot.bovy_text(r'$Z\ =\ %.3f$' % Z,top_right=True,size=14.)
        bovy_plot.bovy_text(r'$[\mathrm{M/H}]\ =\ %.2f$' % (isodist.Z2FEH(Z)),
                                                           top_right=True,size=14.)
        bovy_plot.bovy_end_print(plotfile)
    return None

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-o",dest='plotfile',default=None,
                      help="Name of file for plot")
    parser.add_option("-t","--type",dest='type',default='dens',
                      help="Type of plot ('dens' for density, 'loggmean')")
    parser.add_option("-Z",dest='Z',type='float',default=None,
                      help="Metallicity Z")
    parser.add_option("--ho",dest='ho',type='float',default=12.,
                      help="H_0 to use for distance scale")
    parser.add_option("--dwarf",action="store_true", 
                      dest="dwarf",
                      default=False,
                      help="Show the dwarf part of the isochrone")
    parser.add_option("--log",action="store_true", 
                      dest="log",
                      default=False,
                      help="Use a logarithmic grayscale")
    parser.add_option("--loggash",action="store_true", 
                      dest="loggash",
                      default=False,
                      help="Use logg rather than M_H as the y-axis")
    return parser

if __name__ == '__main__':
    parser= get_options()
    (options,args)= parser.parse_args()
    imf_h_jk(options.plotfile,Z=options.Z,dwarf=options.dwarf,log=options.log,
             h=options.ho,loggash=options.loggash)
