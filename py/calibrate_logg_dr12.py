#Routine to calibrate logg for the DR12 ASPCAP analysis
import sys
import os
import numpy
import fitsio
import esutil
from galpy.util import bovy_plot
#Make sure that we're using the DR12 files
os.environ['APOGEE_REDUX']= 'v601'
os.environ['APOGEE_APOKASC_REDUX']= 'v7.3'
import apogee.tools.read as apread
from apogee.tools import paramIndx
def calibrate_logg_dr12(rgb=False):
    """Calibrate, using RC or RGV when rgb=True (the latter should reproduce the official calibration"""
    #Read the calibration file, APOKASC, and match them
    caldata= fitsio.read(os.path.join(os.getenv('APOGEE_DATA'),'cal_%s.fits' % os.getenv('APOGEE_REDUX')),1)
    apokasc= apread.apokasc()
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(caldata['RA'],caldata['DEC'],
                        apokasc['RA'],apokasc['DEC'],
                        2./3600.,maxmatch=1)
    caldata= caldata[m1]
    apokasc= apokasc[m2]
    #Select a decent calibration sample
    # Use stars that are definitely RGB or RC
    seismoState= numpy.char.strip(apokasc['SEISMO EVOL'])
    if rgb:
        indx= (seismoState == 'RGB') \
            + (seismoState == 'DWARF/SUBGIANT')
    else:
        indx= (seismoState == 'CLUMP')
    #Add low-logg giants of any kind
    indx+= (apokasc['KASC_RG_LOGG_SCALE_2'] < 2.)
    #rm bad data
    indx= (caldata['FPARAM'][:,paramIndx('logg')] > -1000.)
    indx*= (apokasc['KASC_RG_LOGG_SCALE_2'] > -1000.)
    #Apply limits
    indx*= (caldata['FPARAM'][:,paramIndx('logg')] > 1.)\
        *(caldata['FPARAM'][:,paramIndx('logg')] < 3.8)
    print "Using %i stars to calibrate logg ..." % numpy.sum(indx)
    #Now fit the difference
    fitOut= numpy.polyfit(caldata['FPARAM'][indx,paramIndx('logg')],
                          caldata['FPARAM'][indx,paramIndx('logg')]\
                              -apokasc['KASC_RG_LOGG_SCALE_2'][indx],
                          1)
    print fitOut
    if True:
        bovy_plot.bovy_print()
        bovy_plot.bovy_plot(caldata['FPARAM'][indx,paramIndx('logg')],
                            caldata['FPARAM'][indx,paramIndx('logg')]\
                                -apokasc['KASC_RG_LOGG_SCALE_2'][indx],
                            'k.',
                            xrange=[0.,5.],yrange=[-0.5,1.],
                            xlabel=r'$\log g_{\mathrm{ASPCAP}}$',
                            ylabel=r'$\log g_{\mathrm{ASPCAP}}-\log g_{\mathrm{seismo}}$')
        if not rgb:
            plotindx= (seismoState == 'RGB') \
                + (seismoState == 'DWARF/SUBGIANT')
        else:
            plotindx= (seismoState == 'CLUMP')
        plotindx*= (caldata['FPARAM'][:,paramIndx('logg')] > -1000.)
        plotindx*= (apokasc['KASC_RG_LOGG_SCALE_2'] > -1000.)
        bovy_plot.bovy_plot(caldata['FPARAM'][plotindx,paramIndx('logg')],
                            caldata['FPARAM'][plotindx,paramIndx('logg')]\
                                -apokasc['KASC_RG_LOGG_SCALE_2'][plotindx],
                            '.',color='0.65',overplot=True)
        xs= numpy.linspace(1.,3.8,1001)
        bovy_plot.bovy_plot(xs,fitOut[0]*xs+fitOut[1],'k-',overplot=True)
        bovy_plot.bovy_plot(xs,-0.14*xs+0.588,'k--',overplot=True)
        print numpy.amax(numpy.fabs(fitOut[0]*xs+fitOut[1]-(-0.14*xs+0.588)))
        print numpy.amax(numpy.fabs(fitOut[0]*xs+fitOut[1]-(-0.14*xs+0.588))[xs < 2.])
        bovy_plot.bovy_text(r'$\mathrm{diff} = %.3f \log g_{\mathrm{ASPCAP}} + %.3f$' % (fitOut[0],fitOut[1]),
                            bottom_left=True,fontsize=14.)
        bovy_plot.bovy_end_print('/Users/bovy/Desktop/test.png')
    return None

if __name__ == '__main__':
    calibrate_logg_dr12(rgb=(len(sys.argv) == 2))
