import os, os.path
import cPickle as pickle
from optparse import OptionParser
import numpy
from scipy import maxentropy
from galpy.util import multi, save_pickles
import isodist
import rcmodel
def localzdist(z,zsolar=0.019):
    #From 2 Gaussian XD fit to Casagrande et al. (2011)
    feh= isodist.Z2FEH(z,zsolar=zsolar)
    logfehdist= maxentropy.logsumexp([numpy.log(0.8)-numpy.log(0.15)-0.5*(feh-0.016)**2./0.15**2.,
                                      numpy.log(0.2)-numpy.log(0.22)-0.5*(feh+0.15)**2./0.22**2.])
    return logfehdist-numpy.log(z)
                                      
def calc_avg_rcmks(parser):
    options,args= parser.parse_args()
    njks= 101
    nmks= 101
    jks= numpy.linspace(0.5,0.8,njks)
    mks= numpy.linspace(-0.5,-3.,nmks)
    if options.basti:
        zs= numpy.array([0.004,0.008,0.01,0.0198,0.03,0.04])
        zsolar= 0.019
    elif options.parsec:
        zs= numpy.arange(0.0005,0.06005,0.0005)
#        zs= numpy.array([0.01,0.02])
        zsolar= 0.019
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
#        zs= numpy.array([0.01,0.02])
        zsolar= 0.019
    if not os.path.exists(options.outfilename):
        logpz= localzdist(zs,zsolar=zsolar)
        logmkp= numpy.zeros((len(zs),njks,nmks))
        logp= numpy.zeros((len(zs),njks,nmks))      
        funcargs= (zs,options,njks,jks,nmks,mks,logpz)
        multOut= multi.parallel_map((lambda x: indiv_calc(x,
                                                          *funcargs)),
                                    range(len(zs)),
                                    numcores=numpy.amin([64,len(zs),
                                                         multiprocessing.cpu_count()]))
        for ii in range(len(zs)):
            logmkp[ii,:,:]= out[ii][0,:,:]
            logp[ii,:,:]= out[ii][1,:,:]
        save_pickles(options.outfilename,logmkp,logp)
    else:
        savefile= open(options.outfilename,'rb')
        logmkp= pickle.load(savefile)
        logp= pickle.load(savefile)
        savefile.close()
    avgmk= numpy.exp(maxentropy.logsumexp(logmkp.flatten())\
                         -maxentropy.logsumexp(logp.flatten()))
    print "Average mk: %f" % (-avgmk)
    return -avgmk

def indiv_calc(x,zs,options.njks,jks,nmks,mks,logpz):
    print zs[ii]
    out= numpy.empty((2,njks,nmks))
    rc= rcmodel.rcmodel(Z=zs[ii],loggmin=1.8,loggmax=2.8,
                        band=options.band,basti=options.basti,
                        imfmodel=options.imfmodel,
                        expsfh=options.expsfh,
                        parsec=options.parsec,
                        stage=options.stage)
    for jj in range(njks):
        for kk in range(nmks):
            out[0,jj,kk]= logpz[ii]+rc(jks[jj],mks[kk])+numpy.log(-mks[kk])
            out[1,jj,kk]= logpz[ii]+rc(jks[jj],mks[kk])
    return out

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
    parser.add_option("--expsfh",action="store_true", dest="expsfh",
                      default=False,
                      help="If set, use an exponentially-declining SFH")
    parser.add_option("--parsec",action="store_true", dest="parsec",
                      default=False,
                      help="If set, PARSEC")
    parser.add_option("-s",dest='stage',default=None,type='int',
                      help="Evolutionary stage to use")
    return parser
    
if __name__ == '__main__':
    parser= get_options()
    calc_avg_rcmks(parser)
