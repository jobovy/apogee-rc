import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
from scipy import maxentropy
import isodist
import rcmodel
def localzdist(z,zsolar=0.019):
    feh= isodist.Z2FEH(z,zsolar=zsolar)
    logfehdist= maxentropy.logsumexp([numpy.log(0.9)-numpy.log(0.14)-0.5*(feh+0.13)**2./0.14**2.,
numpy.log(0.1)-numpy.log(0.36)-0.5*(feh+0.41)**2./0.36**2.])
    return logfehdist-numpy.log(z)
                                      
def calc_avg_rcmks(parser):
    options,args= parser.parse_args()
    rc= rcmodel.rcmodel(Z=zs[ii],loggmin=1.8,loggmax=2.8,
                        band=options.band,basti=options.basti,
                        imfmodel=options.imfmodel,
                        parsec=options.parsec)
    njks= 101
    nmks= 101
    jks= numpy.linspace(0.5,0.75,njks)
    mks= numpy.linspace(0.5,3.,nmks)
    if options.basti:
        zs= numpy.array([0.004,0.008,0.01,0.0198,0.03,0.04])
        zsolar= 0.0198
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
        zs= numpy.array([0.01,0.02])
        zsolar= 0.019
    logpz= localzdist(zs,zsolar=zsolar)
    logmkp= numpy.zeros((len(zs),njks,nmks))
    logp= numpy.zeros((len(zs),njks,nmks))
    for ii in range(len(zs)):
        for jj in range(njks):
            for kk in range(nmks):
                logmkp[ii,jj,kk]= logpz[ii]+rc(jks[jj],mks[kk])+numpy.log(mks[kk])
                logp[ii,jj,kk]= logpz[ii]+rc(jks[jj],mks[kk])
    avgmk= numpy.exp(maxentropy.logsumexp(logmkp.flatten())\
                         -maxentropy.logsumexp(logp.flatten()))
    print "Average mk: %f" % avgmk
    return avgmk

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
    parser.add_option("--parsec",action="store_true", dest="parsec",
                      default=False,
                      help="If set, PARSEC")
    return parser
    
if __name__ == '__main__':
    parser= get_options()
    plot_vs_jkz(parser)
