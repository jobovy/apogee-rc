import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
_EXT='png'
_ADDLLOGGCUT= True
_SCATTERPLOT= False
def plot_fehafe(basefilename):
    #First read the sample
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #R and z ranges to plot
    Rmins= [5.,7.,9.]
    Rmaxs= [7.,9.,11.]
    zmins= [0.,0.5,1.]
    zmaxs= [0.5,1.,2.]
    nRbins= len(Rmins)
    nzbins= len(zmins)
    for ii in range(nRbins):
        for jj in range(nzbins):
            #Cut data to relevant range
            tindx= (data['RC_GALR'] > Rmins[ii])\
                *(data['RC_GALR'] <= Rmaxs[ii])\
                *(numpy.fabs(data['RC_GALZ']) > zmins[jj])\
                *(numpy.fabs(data['RC_GALZ']) <= zmaxs[jj])
            tdata= data[tindx]
            bovy_plot.bovy_print()
            if _SCATTERPLOT or len(tdata) > 2000:
                bovy_plot.scatterplot(tdata['METALS'],
                                      tdata['ALPHAFE'],
                                      'k.',
                                      xrange=[-1.2,0.7],
                                      yrange=[-0.12,0.42],
                                      xlabel=r'$[\mathrm{Fe/H}]$',
                                      ylabel=r'$[\alpha/\mathrm{H}]$',
                                      onedhists=True,bins=31,onedhistsbins=31)
            else:
                bovy_plot.bovy_plot(tdata['METALS'],
                                    tdata['ALPHAFE'],
                                    'k.',
                                    xrange=[-1.2,0.7],
                                    yrange=[-0.12,0.42],
                                    xlabel=r'$[\mathrm{Fe/H}]$',
                                    ylabel=r'$[\alpha/\mathrm{H}]$',
                                    onedhists=True,bins=21)
            #Overplot ridge line
            def ridge1(feh):
                return -0.2*(feh+0.6)+0.23
            bovy_plot.bovy_plot([-0.8,0.3],[ridge1(-0.8),ridge1(0.3)],
                                '-',lw=2.,color='y',overplot=True)
            bovy_plot.bovy_text(r'$%i < R / \mathrm{kpc} \leq %i$' % (Rmins[ii],Rmaxs[ii]) +'\n'+r'$%.1f < |z| / \mathrm{kpc} \leq %.1f$' % (zmins[jj],zmaxs[jj]),
                                top_left=True,size=16.)
            bovy_plot.bovy_end_print(basefilename+'_%iR%i_%.1fz%.1f.' % (Rmins[ii],Rmaxs[ii],zmins[jj],zmaxs[jj])+_EXT)
    return None

if __name__ == '__main__':
    plot_fehafe(sys.argv[1])
