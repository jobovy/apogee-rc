import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
import isodist
from galpy.util import bovy_plot
import rcmodel
def plot_mx_jkz(parser):
    options,args= parser.parse_args()
    rc= rcmodel.rcmodel(Z=options.Z,loggmin=1.8,loggmax=2.8,
                        basti=options.basti,band=options.band,
                        imfmodel=options.imfmodel)
    #Calculate mode and hm
    njks= 101
    jks= numpy.linspace(0.5,0.75,njks)
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
    if options.basti:
        zsolar= 0.0198
    rc.plot(nbins=101,conditional=True)
    bovy_plot.bovy_plot(jks,modes,'w-',lw=2.,overplot=True)
    bovy_plot.bovy_plot(jks,hms[:,0],'-',lw=2.,color='0.85',overplot=True)
    bovy_plot.bovy_plot(jks,hms[:,1],'-',lw=2.,color='0.85',overplot=True)
    bovy_plot.bovy_text(r'$[\mathrm{M/H}]\ =\ %.2f$' % (isodist.Z2FEH(options.Z)),
                        bottom_left=True,size=14.)
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
    return parser
    
if __name__ == '__main__':
    parser= get_options()
    plot_mx_jkz(parser)
    
