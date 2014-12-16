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
    lpower= large_scale_power(pixrc,vsolars,vc=220.,dx=_RCDX)
    bovy_plot.bovy_print()
    line1= bovy_plot.bovy_plot(vsolars,lpower,'k-',lw=2.,
                               xrange=[vsolars[0],vsolars[-1]],
                               yrange=[0.,22.9],
                               xlabel=r'$V_{\odot-c}\,(\mathrm{km\,s}^{-1})$',
                               ylabel=r'$\sqrt\langle P_k(0.2 < k / (\mathrm{kpc}^{-1}) < 0.9)\rangle\,(\mathrm{km\,s}^{-1})$')
    #Find the minimum by fitting a second order polynomial
    p= numpy.polyfit(vsolars,lpower,2)
    minvsolar= -0.5*p[1]/p[0]
    bovy_plot.bovy_plot([minvsolar,minvsolar],[-10.,100.],'k-',lw=0.8,
                        overplot=True)
    bovy_plot.bovy_text(22.65,1.,
                        r'$V_{\odot-c}=%.1f\,\mathrm{km\,s}^{-1}$' % minvsolar,
                        size=15.)
    lpower240= large_scale_power(pixrc,vsolars,vc=240.,dx=_RCDX)
    line2= bovy_plot.bovy_plot(vsolars,lpower240,'k--',lw=2.,overplot=True)
    lpower200= large_scale_power(pixrc,vsolars,vc=200.,dx=_RCDX)
    line3= bovy_plot.bovy_plot(vsolars,lpower200,'k-.',lw=2.,overplot=True)
    lpowerbetam0p2= large_scale_power(pixrc,vsolars,dx=_RCDX,beta=-0.1)
    line4= bovy_plot.bovy_plot(vsolars,lpowerbetam0p2,
                               '--',dashes=(20,10),
                               color='0.6',lw=3.,overplot=True)
    lpowerbetap0p2= large_scale_power(pixrc,vsolars,dx=_RCDX,beta=0.1)
    line5= bovy_plot.bovy_plot(vsolars,lpowerbetap0p2,
                               '--',dashes=(20,10),
                               color='0.4',lw=3.,overplot=True)
    lva= large_scale_power(pixrc,vsolars,vc=220.,dx=_RCDX,hs=6./8.,hR=3./8.)
    line6= bovy_plot.bovy_plot(vsolars,lva,
                               '-',color='0.8',lw=5.,overplot=True,zorder=-1)
    #Add legend
    legend1= pyplot.legend((line3[0],line1[0],line2[0]),
                           (r'$V_c = 200\,\mathrm{km\,s}^{-1}$',
                            r'$V_c = 220\,\mathrm{km\,s}^{-1}$',
                            r'$V_c = 240\,\mathrm{km\,s}^{-1}$'),
                           loc='upper left',#bbox_to_anchor=(.91,.375),
                           numpoints=2,
                           prop={'size':14},
                  frameon=False)
    legend2= pyplot.legend((line4[0],line5[0]),
                           (r'$\beta = -0.1$',
                            r'$\beta = \phantom{-}0.1$'),
                           loc='upper right',bbox_to_anchor=(1.,.96),
                           numpoints=2,
                           prop={'size':14},
                           frameon=False)
    legend3= pyplot.legend((line6[0],),
                           (r'$\Delta V_a(R_0) =5\,\mathrm{km\,s}^{-1}$',),
                           loc='lower left',#bbox_to_anchor=(1.,.96),
                           numpoints=2,
                           prop={'size':12},
                           frameon=False)
    pyplot.gca().add_artist(legend1)
    pyplot.gca().add_artist(legend2)
    pyplot.gca().add_artist(legend3)
    bovy_plot.bovy_end_print(plotfilename)
    return None

def large_scale_power(pix,vsolar,vc=220.,dx=None,beta=0.,hs=33.3,hR=3./8.):
    """Determine the power on large scales in the residuals for different solarmotions"""
    out= numpy.empty_like(vsolar)
    binsize= 0.8
    scale= 4.*numpy.pi #0.522677552224
    for ii in range(len(vsolar)):
        resv= pix.plot(lambda x: dvlosgal(x,vc=vc,vtsun=vc+vsolar[ii],
                                          beta=beta,hR=hR,hs=hs),
                       returnz=True,justcalc=True)
        psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
        indx= (psd1d[0] > 0.2)*(psd1d[0] < 0.9)
        out[ii]= scale*numpy.sqrt(numpy.sum(psd1d[1][indx]*(psd1d[1][indx]/psd1d[2][indx])**2.)/numpy.sum((psd1d[1][indx]/psd1d[2][indx])**2.))
    return out                      

if __name__ == '__main__':
    determine_vsolar(sys.argv[1])
