import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
import isodist
from galpy.util import bovy_plot
from matplotlib import pyplot
import rcmodel
_NEW= True
if _NEW:
    loggmax='custom'
else:
    loggmax= 2.8
def plot_mx_jkz(parser):
    options,args= parser.parse_args()
    rc= rcmodel.rcmodel(Z=options.Z,loggmin=1.8,loggmax=loggmax,
                        basti=options.basti,band=options.band,
                        parsec=options.parsec,
                        expsfh=options.expsfh,
                        imfmodel=options.imfmodel,
                        eta=options.eta,afe=options.afe)
    #Calculate mode and hm
    njks= 101
    jks= numpy.linspace(0.5,0.8,njks)
    modes= numpy.array([rc.mode(jk) for jk in jks])
    hms= numpy.zeros((njks,2))
    for ii in range(njks):
        try:
            minhm, maxhm= rc.sigmafwhm(jks[ii],straight=True)
        except ValueError:
            minhm, maxhm= numpy.nan, numpy.nan
        hms[ii,0]= minhm
        hms[ii,1]= maxhm
    #Now plot
    bovy_plot.bovy_print(text_fontsize=20.,
                         legend_fontsize=24.,
                         xtick_labelsize=18.,
                         ytick_labelsize=18.,
                         axes_labelsize=24.)
    rc.plot(nbins=101,conditional=True)
    bovy_plot.bovy_plot(jks,modes,'w-',lw=2.,overplot=True)
    bovy_plot.bovy_plot(jks,hms[:,0],'-',lw=2.,color='0.85',overplot=True)
    bovy_plot.bovy_plot(jks,hms[:,1],'-',lw=2.,color='0.85',overplot=True)
    #Overplot the cuts in J-Ks for this Z
    bovy_plot.bovy_plot([rcmodel.zjkcut(options.Z),rcmodel.zjkcut(options.Z)],
                        [0.,-3.],'k--',lw=2.,overplot=True)
    bovy_plot.bovy_plot([rcmodel.zjkcut(options.Z,upper=True),rcmodel.zjkcut(options.Z,upper=True)],
                        [0.,-3.],'k--',lw=2.,overplot=True)
    if options.Z < 0.01: zstr= r'$Z = %.3f$' % options.Z
    else: zstr= r'$Z = %.2f$' % options.Z
    if options.afe:
        bovy_plot.bovy_text(zstr+'\n'+r'$[\alpha/\mathrm{Fe}] = 0.4$',
                            bottom_right=True,size=20.)
    else:
        bovy_plot.bovy_text(zstr,
                            bottom_right=True,size=20.)
    if options.basti and options.Z < 0.01:
        pyplot.annotate(r'$\mathrm{BaSTI}\rightarrow$',(0.5,1.08),
                        xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=24.)
        #bovy_plot.bovy_text(r'$\mathrm{BaSTI}$',title=True,size=20.)
    elif options.parsec and options.Z < 0.01:
        pyplot.annotate(r'$\mathrm{PARSEC}\rightarrow$',(0.5,1.08),
                        xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=24.)
        #bovy_plot.bovy_text(r'$\mathrm{PARSEC}$',title=True,size=20.)
    elif options.Z < 0.01:
        pyplot.annotate(r'$\mathrm{Padova}\rightarrow$',(0.5,1.08),
                        xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=24.)
        #bovy_plot.bovy_text(r'$\mathrm{Padova}$',title=True,size=20.)
    bovy_plot.bovy_end_print(options.outfilename)
    return None

def get_options():
    usage = "usage: %prog [options] savefilename"
    parser = OptionParser(usage=usage)
    #Data options
    parser.add_option("-o",dest='outfilename',default=None,
                      help="Name for plot file")
    parser.add_option("-z",dest='Z',default=0.019,type='float',
                      help="Metallicity to use")
    parser.add_option("-b",dest='band',default='Ks',
                      help="Band to use")
    parser.add_option("--imfmodel",dest='imfmodel',
                      default='lognormalChabrier2001',
                      help="IMF model to use")
    parser.add_option("--basti",action="store_true", dest="basti",
                      default=False,
                      help="If set, use BaSTI isochrones")
    parser.add_option("--parsec",action="store_true", dest="parsec",
                      default=False,
                      help="If set, PARSEC")
    parser.add_option("--expsfh",action="store_true", dest="expsfh",
                      default=False,
                      help="If set, use an exponentially-declining SFH")
    parser.add_option("--eta",dest='eta',default=None,type='float',
                      help="Mass-loss efficiency parameter")
    parser.add_option("--afe",action="store_true", dest="afe",
                      default=False,
                      help="If set, use alpha-enhanced isochrones (works only for Basti, for now")
    return parser
    
if __name__ == '__main__':
    parser= get_options()
    plot_mx_jkz(parser)
    
