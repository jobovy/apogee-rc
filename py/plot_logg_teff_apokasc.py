import os, os.path
import cPickle as pickle
import copy
import numpy
from galpy.util import bovy_plot, save_pickles
import isodist
import apogee.tools.read as apread
import rcmodel
from plot_vs_jkz import get_options
_DEBUG= False
def plot_logg_teff_apokasc(parser):
    options,args= parser.parse_args()
    #Setup Zs
    if os.path.exists(args[0]):
        savefile= open(args[0],'rb')
        outhist= pickle.load(savefile)
        hists= pickle.load(savefile)
        edgess= pickle.load(savefile)
        data= pickle.load(savefile)
        savefile.close()
    else:
        zs= numpy.arange(0.0005,0.03005,0.0005)
        if _DEBUG:
            zs= numpy.arange(0.0005,0.03005,0.005)
        fehs= isodist.Z2FEH(zs,zsolar=0.017)#0.017 for APOGEE analysis
        zs= zs[numpy.fabs(fehs-options.feh) < 0.2]
        if _DEBUG:
            zs= [zs[numpy.fabs(fehs-options.feh) < 0.2][0]]
        fehs= isodist.Z2FEH(zs,zsolar=0.017)   
        #Load the RC models for each feh individually
        rcms= []
        hists= []
        edgess= []
        for z in zs:
            print z
            trc= rcmodel.rcmodel(Z=z,loggmin=1.,loggmax=3.5,
                                 band=options.band,basti=options.basti,
                                 imfmodel=options.imfmodel,
                                 parsec=options.parsec)
            rcms.append(trc)
            sample= numpy.vstack([trc._teffs[:,0],trc._loggs[:,0]]).T
            weights= trc._weights
            hist, edges= numpy.histogramdd(sample,weights=weights,bins=31,
                                           range=[[4450.,5200.],[1.,3.5]])
            hists.append(hist)
            edgess.append(edges)
        #Load APOKASC data
        data= apread.apokasc()
        indx= (data['KASC_RG_LOGG_SCALE_2'] > 1.)\
            *(data['KASC_RG_LOGG_SCALE_2'] < 3.5)\
            *(data['METALS'] > options.feh-0.2)\
            *(data['METALS'] <= options.feh+0.2)
        print "Using %i APOKASC objects ..." % (numpy.sum(indx))
        data= data[indx]
        ndata= numpy.sum(indx)
        #Stack predictions
        outhist= numpy.zeros_like(hists[0])
        for ii in range(ndata):
            zindx= numpy.argmin(numpy.fabs(fehs-data['METALS'][ii]))
            outhist+= hists[zindx]
        save_pickles(args[0],outhist,hists,edgess,data)
    #Normalize each Teff
    for ii in range(len(outhist[:,0])):
        outhist[ii,:]/= numpy.nanmax(outhist[ii,:])/numpy.nanmax(outhist)
        rev= copy.copy(outhist[ii,::-1]) #reverse, but in one go does not always work
        outhist[ii,:]= rev
    if False:
        #Reload apokasc data
        data= apread.apokasc()
        indx= (data['KASC_RG_LOGG_SCALE_2'] > 1.)\
            *(data['KASC_RG_LOGG_SCALE_2'] < 3.5)\
            *(data['METALS'] > options.feh-0.2)\
            *(data['METALS'] <= options.feh+0.2)
        data= data[indx]
    #Plot everything
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(outhist.T,origin='lower',cmap='gist_yarg',
                          xrange=[edgess[0][0][0],edgess[0][0][-1]],
                          yrange=[edgess[0][1][-1],edgess[0][1][0]],
                          aspect=(edgess[0][0][-1]-edgess[0][0][0])/float(edgess[0][1][-1]-edgess[0][1][0]),
                          xlabel=r'$T_{\mathrm{eff}}\,(\mathrm{K})$',
                          ylabel=r'$\mathrm{Seismic}\ \log g$',
                          shrink=0.78,
                          interpolation='nearest')
    #Overplot APOKASC data
    noseismo= data['SEISMO EVOL'] == 'UNKNOWN'
    if numpy.sum(noseismo) > 0:
        bovy_plot.bovy_plot(data['TEFF'][noseismo],
                            data['KASC_RG_LOGG_SCALE_2'][noseismo],'bo',
                            overplot=True,
                            mec='none',ms=3.)
    clumpseismo= data['SEISMO EVOL'] == 'CLUMP'
    if numpy.sum(clumpseismo) > 0:
        bovy_plot.bovy_plot(data['TEFF'][clumpseismo],
                            data['KASC_RG_LOGG_SCALE_2'][clumpseismo],'yo',
                            overplot=True,
                            mec='none',ms=4.5)
    noclumpseismo= (data['SEISMO EVOL'] == 'RGB') \
        + (data['SEISMO EVOL'] == 'DWARF/SUBGIANT')
    if numpy.sum(noclumpseismo) > 0:
        bovy_plot.bovy_plot(data['TEFF'][noclumpseismo],
                            data['KASC_RG_LOGG_SCALE_2'][noclumpseismo],'ro',
                            overplot=True,
                            mec='none',ms=3.)
    bovy_plot.bovy_text(r'$%.1f < [\mathrm{M/H}] \leq %.1f$' % (options.feh-0.2,options.feh+0.2),
                        top_left=True,size=14.)
    bovy_plot.bovy_end_print(options.outfilename)
    return None
    
if __name__ == '__main__':
    parser= get_options()
    plot_logg_teff_apokasc(parser)
    
