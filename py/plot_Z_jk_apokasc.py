import os, os.path
import cPickle as pickle
import copy
import numpy
from galpy.util import bovy_plot, save_pickles
import isodist
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import apogee.tools.read as apread
import rcmodel
from plot_vs_jkz import get_options
_DEBUG= True
_PREDICT= True
def plot_Z_jk_apokasc(parser):
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
        if _PREDICT:
            zs= numpy.arange(0.0005,0.06005,0.0005)
            if _DEBUG:
                zs= numpy.arange(0.0005,0.06005,0.005)
            #Load the RC models for each feh individually
            rcms= []
            hists= []
            edgess= []
            for z in zs:
                print z
                trc= rcmodel.rcmodel(Z=z,loggmin=1.8,loggmax=2.8,
                                     band=options.band,basti=options.basti,
                                     imfmodel=options.imfmodel,
                                     parsec=options.parsec)
                rcms.append(trc)
                sample= numpy.vstack([trc._sample[:,0],
                                      z*numpy.ones(trc._sample.shape[0])]).T
                weights= trc._weights
                hist, edges= numpy.histogramdd(sample,weights=weights,
                                               bins=12,#12*(10-_DEBUG*9),
                                               range=[[0.5,0.8],[0.,0.06]])
                hists.append(hist)
                edgess.append(edges)
        #Load APOKASC data
        data= apread.apokasc()
        indx= (data['KASC_RG_LOGG_SCALE_2'] > 2.25)\
            *(data['KASC_RG_LOGG_SCALE_2'] < 2.65)
        print "Using %i APOKASC objects ..." % (numpy.sum(indx))
        data= data[indx]
        if _PREDICT:
            #Stack predictions
            outhist= numpy.zeros_like(hists[0])
            for ii in range(len(hists)):
                outhist+= hists[ii]
            save_pickles(args[0],outhist,hists,edgess,data)
    if _PREDICT:
        #Normalize each color
        pass
        #        for ii in range(len(outhist[:,0])):
#            outhist[ii,:]/= numpy.nanmax(outhist[ii,:])/numpy.nanmax(outhist)
        #for ii in range(len(outhist[0,:])):
        #    outhist[:,ii]/= numpy.nanmax(outhist[:,ii])/numpy.nanmax(outhist)
    #Plot everything
    bovy_plot.bovy_print()
    if _PREDICT:
        bovy_plot.bovy_dens2d(outhist.T,origin='lower',cmap='gist_yarg',
                              xrange=[edgess[0][0][0],edgess[0][0][-1]],
                              yrange=[edgess[0][1][0],edgess[0][1][-1]],
                              aspect=(edgess[0][0][-1]-edgess[0][0][0])/float(edgess[0][1][-1]-edgess[0][1][0]),
                              xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                              ylabel=r'$Z$',
                              shrink=0.78,
                              interpolation='nearest')
    #Overplot APOKASC data
        #Load APOKASC data
        data= apread.apokasc()
        indx= (data['KASC_RG_LOGG_SCALE_2'] > 1.8)\
            *(data['KASC_RG_LOGG_SCALE_2'] < 2.8)
        print "Using %i APOKASC objects ..." % (numpy.sum(indx))
        data= data[indx]
    bovy_plot.bovy_plot(data['J0']-data['K0'],
                        options.zsolar*10.**data['METALS'],
                        c=data['KASC_RG_LOGG_SCALE_2']-2.45,
                        s=20.,edgecolors='none',
                        scatter=True,colorbar=True,
                        overplot=True)
#                        mec='none',ms=3.)
    #Overplot cuts
    jks= numpy.linspace(0.5,0.8,201)
    bovy_plot.bovy_plot(jks,rcmodel.jkzcut(jks),
                        'k--',lw=2.,overplot=True)
    bovy_plot.bovy_plot(jks,rcmodel.jkzcut(jks,upper=True),
                        'k--',lw=2.,overplot=True)
    bovy_plot.bovy_end_print(options.outfilename)
    return None
    
if __name__ == '__main__':
    parser= get_options()
    plot_Z_jk_apokasc(parser)
    
