# Run with python plot_bird_psd.py ~/Desktop/test.png
import sys
import numpy
import bovy_psd
from plot_psd import _RCDX, _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, \
    _NNOISE, _SUBTRACTERRORS
from simple_spiral_simulation import simulate_vlos_spiral
from galpy.util import bovy_plot
from matplotlib import pyplot
_birdFile= '../pecvel/pecvel.npz'
_nSims= 8
_PLOTINDIV= True
_ADDDATALINE= True
def plot_bird_psd(plotfilename):
    #Read the Bird data
    birdData= numpy.load('../pecvel/pecvel.npz')
    #Get residuals for all simulations
    dx= _RCDX
    binsize= .8#.765
    scale= 4.*numpy.pi
    tmp= bovy_psd.psd1d(birdData['dVlos1'],dx,binsize=binsize) #just to get the size
    ks= tmp[0][1:-3]
    psds= numpy.zeros((len(tmp[1]),_nSims))
    if _SUBTRACTERRORS:
        for ii in range(_nSims):
            sim= ii+1
            tmpPsd= bovy_psd.psd1d(birdData['dVlos%i' % sim],
                                   dx,binsize=binsize)[1]
            #Simulations for the noise
            nnoise= _NNOISE
            noisepsd= numpy.empty((nnoise,len(tmpPsd)))
            for jj in range(nnoise):
                newresv= \
                    numpy.random.normal(size=birdData['dVlos%i' % sim].shape)\
                    *birdData['sig_dVlos%i' % sim].reshape((9,9))\
                    *(True-birdData['rc_mask'])
                noisepsd[jj,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1]
            psds[:,ii]= tmpPsd-numpy.median(noisepsd,axis=0)
    #Calculate median PSD and spread around this
    medPsd= scale*numpy.median(numpy.sqrt(psds),axis=1)[1:-3]
    flucPsd=\
        1.4826*scale*numpy.median(numpy.fabs(numpy.sqrt(psds)[1:-3]
                                             -numpy.tile(medPsd/scale,
                                                         (psds.shape[1],1)).T),axis=1)
    print medPsd, flucPsd
    #Now plot
    xrange=[.03,3.]
    yrange= [0.,11.9]
    bovy_plot.bovy_print(fig_width=5.5,fig_height=4.5)
    bovy_plot.bovy_plot(ks,medPsd,
                        'k-',lw=2.,
                        zorder=12,
                        xlabel=r'$k\,(\mathrm{kpc}^{-1})$',
                        ylabel=r'$\sqrt{P_k}\,(\mathrm{km\,s}^{-1})$',
                        semilogx=True,
                        xrange=xrange,yrange=yrange)
    goodIndx= True-numpy.isnan(flucPsd)
    pyplot.fill_between(ks[goodIndx],(medPsd-flucPsd)[goodIndx],
                        y2=(medPsd+flucPsd)[goodIndx],
                        color='0.65',zorder=1)
    if _PLOTINDIV:
        for ii in range(_nSims):
            bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psds[:,ii])[1:-3],
                                '-',color='0.8',overplot=True)
    if _ADDDATALINE:
        alpha= -12.5
        spvlos= simulate_vlos_spiral(alpha=alpha,
                                     gamma=1.2,
                                     xmin=_RCXMIN,xmax=_RCXMAX,
                                     ymin=_RCYMIN,ymax=_RCYMAX,
                                     dx=0.01)
        potscale= 1.35
        simpsd1d= bovy_psd.psd1d(spvlos*220.*potscale,0.01,binsize=binsize)
        tks= simpsd1d[0][1:-3]
        line1= bovy_plot.bovy_plot(tks,
                                   scale*numpy.sqrt(simpsd1d[1][1:-3]),
                                   'k--',lw=2.,overplot=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_bird_psd(sys.argv[1])
    
