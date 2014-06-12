#Determine te solar motion by minimizing the power on large scales
import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
import pixelize_sample
from plot_2dkinematics import dvlosgal
from plot_psd import _ADDLLOGGCUT, \
    _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
import bovy_psd
def determine_vsolar(plotfilename):
    #Read the APOGEE-RC data and pixelate it
    #APOGEE-RC
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    data= data[indx]
    #Get velocity field
    pixrc= pixelize_sample.pixelXY(data,
                                   xmin=_RCXMIN,xmax=_RCXMAX,
                                   ymin=_RCYMIN,ymax=_RCYMAX,
                                   dx=_RCDX,dy=_RCDX)
    vsolars= numpy.linspace(0.,40.,51)
    lpower= large_scale_power(pixrc,vsolars,vc=218.,dx=_RCDX)
    bovy_plot.bovy_print()
    line1= bovy_plot.bovy_plot(vsolars,lpower,'k-',lw=2.,
                               xrange=[0.,40.],
                               yrange=[0.,30.],
                               xlabel=r'$V_\odot\,(\mathrm{km\,s}^{-1})$',
                               ylabel=r'$\sqrt\langle P_k(0.2 < k / (\mathrm{kpc}^{-1}) < 0.9)\rangle\,(\mathrm{km\,s}^{-1})$')
    #Find the minimum by fitting a second order polynomial
    p= numpy.polyfit(vsolars,lpower,2)
    minvsolar= -0.5*p[1]/p[0]
    bovy_plot.bovy_plot([minvsolar,minvsolar],[-10.,100.],'k-',lw=0.8,
                        overplot=True)
    bovy_plot.bovy_text(24.25,1.,
                        r'$V_\odot=%.1f\,\mathrm{km\,s}^{-1}$' % minvsolar,
                        size=15.)
    lpower240= large_scale_power(pixrc,vsolars,vc=240.,dx=_RCDX)
    line2= bovy_plot.bovy_plot(vsolars,lpower240,'k--',lw=2.,overplot=True)
    lpower200= large_scale_power(pixrc,vsolars,vc=200.,dx=_RCDX)
    line3= bovy_plot.bovy_plot(vsolars,lpower200,'k-.',lw=2.,overplot=True)
    #Add legend
    pyplot.legend((line3[0],line1[0],line2[0]),
                  (r'$V_c = 200\,\mathrm{km\,s}^{-1}$',
                   r'$V_c = 218\,\mathrm{km\,s}^{-1}$',
                   r'$V_c = 240\,\mathrm{km\,s}^{-1}$'),
                  loc='upper left',#bbox_to_anchor=(.91,.375),
                  numpoints=2,
                  prop={'size':14},
                  frameon=False)
    bovy_plot.bovy_end_print(plotfilename)
    return None

def large_scale_power(pix,vsolar,vc=218.,dx=None):
    """Determine the power on large scales in the residuals for different solarmotions"""
    out= numpy.empty_like(vsolar)
    binsize= 0.8
    scale= 0.522677552224
    for ii in range(len(vsolar)):
        resv= pix.plot(lambda x: dvlosgal(x,vc=vc,vtsun=vc+vsolar[ii]),
                       returnz=True,justcalc=True)
        psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
        indx= (psd1d[0] > 0.2)*(psd1d[0] < 0.9)
        out[ii]= scale*numpy.sqrt(numpy.sum(psd1d[1][indx]*(psd1d[1][indx]/psd1d[2][indx])**2.)/numpy.sum((psd1d[1][indx]/psd1d[2][indx])**2.))
    return out                      

if __name__ == '__main__':
    determine_vsolar(sys.argv[1])
