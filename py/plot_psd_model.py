import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import bovy_psd
from plot_psd import _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
import galpy_simulations
_HIVRESSTR= '_hivres'
def plot_psd_model(plotfilename,type):
    #Load fiducial
    spvlos= galpy_simulations.vlos('../sim/spiral_rect_omegas0.33_alpha-14%s.sav' % _HIVRESSTR)
    potscale= 0.85
    simpsd1d= bovy_psd.psd1d(spvlos*potscale,0.33333333,binsize=0.8)
    tks= simpsd1d[0][1:-3]
    xrange=[.08,3.]
    if type.lower() == 'elliptical':
        eres= 31
        p= 0.
        vloscp= galpy_simulations.vlos_elliptical(res=eres,cp=0.02,sp=0.,p=p)
        vlossp= galpy_simulations.vlos_elliptical(res=eres,sp=0.02,cp=0.,p=p)
        vloscpsp= galpy_simulations.vlos_elliptical(res=eres,p=p,
                                                    sp=0.02/numpy.sqrt(2.),
                                                    cp=0.02/numpy.sqrt(2.))
        p=2.
        vloscpp2= galpy_simulations.vlos_elliptical(res=eres,cp=0.01,sp=0.,p=p)
        vlosspp2= galpy_simulations.vlos_elliptical(res=eres,sp=0.01,cp=0.,p=p)
        vloscpspp2= galpy_simulations.vlos_elliptical(res=eres,p=p,
                                                      sp=0.01/numpy.sqrt(2.),
                                                      cp=0.01/numpy.sqrt(2.))
        p=-3.
        vloscppm3= galpy_simulations.vlos_elliptical(res=eres,cp=0.05,sp=0.,p=p)
        vlossppm3= galpy_simulations.vlos_elliptical(res=eres,sp=0.05,cp=0.,p=p)
        vloscpsppm3= galpy_simulations.vlos_elliptical(res=eres,p=p,
                                                       sp=0.05/numpy.sqrt(2.),
                                                       cp=0.05/numpy.sqrt(2.))
        xgrid= numpy.linspace((_RCXMIN-8.)/8.+_RCDX/8./2.,
                              (_RCXMAX-8.)/8.-_RCDX/8./2.,
                              eres)
        dx= (xgrid[1]-xgrid[0])*8.
        psdcp= bovy_psd.psd1d(vloscp,dx,binsize=0.8)
        psdsp= bovy_psd.psd1d(vlossp,dx,binsize=0.8)
        psdcpsp= bovy_psd.psd1d(vloscpsp,dx,binsize=0.8)
        psdcpp2= bovy_psd.psd1d(vloscpp2,dx,binsize=0.8)
        psdspp2= bovy_psd.psd1d(vlosspp2,dx,binsize=0.8)
        psdcpspp2= bovy_psd.psd1d(vloscpspp2,dx,binsize=0.8)
        psdcppm3= bovy_psd.psd1d(vloscppm3,dx,binsize=0.8)
        psdsppm3= bovy_psd.psd1d(vlossppm3,dx,binsize=0.8)
        psdcpsppm3= bovy_psd.psd1d(vloscpsppm3,dx,binsize=0.8)
        ks= psdcp[0][1:-3]
        scale= 4.*numpy.pi*220.
        bovy_plot.bovy_print(fig_width=7.,axes_labelsize=20)
        line1= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdcp[1][1:-3]),
                                   'k-',lw=2.,
                                   semilogx=True,
                                   xlabel=r'$k\,(\mathrm{kpc}^{-1})$',
                                   ylabel=r'$\sqrt{P_k}\,(\mathrm{km\,s}^{-1})$',
                                   xrange=xrange,
                                   yrange=[0.,11.9],zorder=1)
        line2= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdsp[1][1:-3]),
                                   'k--',lw=2.,
                                   overplot=True,zorder=1)
        line3= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdcpsp[1][1:-3]),
                                   'k-.',lw=2.,zorder=1,
                                   overplot=True)
        line4= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdcpp2[1][1:-3]),
                                   '-',lw=2.,color='y',zorder=1,
                                   overplot=True)
        line5= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdspp2[1][1:-3]),
                                   '--',lw=2.,color='y',zorder=1,
                                   overplot=True)
        line6= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdcpspp2[1][1:-3]),
                                   '-.',lw=2.,color='y',zorder=1,
                                   overplot=True)
        line7= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdcppm3[1][1:-3]),
                                   '-',lw=2.,color='r',zorder=1,
                                   overplot=True)
        line8= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdsppm3[1][1:-3]),
                                   '--',lw=2.,color='r',zorder=1,
                                   overplot=True)
        line9= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdcpsppm3[1][1:-3]),
                                   '-.',lw=2.,color='r',zorder=1,
                                   overplot=True)
        pyplot.annotate(r'$\mathrm{Elliptical\ perturbation}\ (m=2\ \mathrm{mode})$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=20.)
        l1= pyplot.legend((line1[0],line2[0],line3[0]),
                          (r'$\phi_b = 0^\circ$',
                           r'$\phi_b = 45^\circ$',
                           r'$\phi_b = 90^\circ$'),
                          loc='lower left',#bbox_to_anchor=(.91,.375),
                          numpoints=8,
                          prop={'size':16},
                          frameon=False)
        l2= pyplot.legend((line1[0],line4[0],line7[0]),
                          (r'$\epsilon(R) = 0.02$',
                           r'$\epsilon(R) = 0.01\,\left(\frac{R}{R_0}\right)^2$',
                           r'$\epsilon(R) = 0.05\,\left(\frac{R}{R_0}\right)^{-3}$'),
                          loc='upper right',#bbox_to_anchor=(.91,.375),
                          numpoints=8,
                          prop={'size':16},
                          frameon=False)
        pyplot.gca().add_artist(l1)
        pyplot.gca().add_artist(l2)
    elif type.lower() == 'bar':
        vlosbar= galpy_simulations.vlos('../sim/bar_rect%s.sav' % _HIVRESSTR)
        vlosslowbar= galpy_simulations.vlos('../sim/bar_rect_slow%s.sav' % _HIVRESSTR)
        eres= 19
        xgrid= numpy.linspace((_RCXMIN-8.)/8.+_RCDX/8./2.,
                              (_RCXMAX-8.)/8.-_RCDX/8./2.,
                              eres)
        dx= (xgrid[1]-xgrid[0])*8.
        psdbar= bovy_psd.psd1d(vlosbar,dx,binsize=0.8)
        psdslowbar= bovy_psd.psd1d(vlosslowbar,dx,binsize=0.8)
        ks= psdbar[0][1:-3]
        scale= 4.*numpy.pi*220.
        bovy_plot.bovy_print(fig_width=7.,axes_labelsize=20)
        line1= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdbar[1][1:-3]),
                                   'r-',lw=2.,
                                   semilogx=True,
                                   xlabel=r'$k\,(\mathrm{kpc}^{-1})$',
                                   ylabel=r'$\sqrt{P_k}\,(\mathrm{km\,s}^{-1})$',
                                   xrange=xrange,
                                   yrange=[0.,11.9],zorder=1)
        line2= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psdslowbar[1][1:-3]),
                                   'y-',lw=2.,
                                   overplot=True,zorder=1)
        pyplot.annotate(r'$\mathrm{Bar\ perturbation\ (rotating}\ m=2\ \mathrm{mode})$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=20.)
        l1= pyplot.legend((line1[0],line2[0]),
                          (r'$\mathrm{Fast\ bar\ growth}$',
                           r'$\mathrm{Slow\ bar\ growth}$'),
                          loc='upper right',#bbox_to_anchor=(.91,.375),
                          numpoints=8,
                          prop={'size':16},
                          frameon=False)    
    elif type.lower() == 'spiral':
        vlosfid= galpy_simulations.vlos('../sim/spiral_rect_omegas0.33_alpha-14%s.sav' % _HIVRESSTR)
        vloslpitch= galpy_simulations.vlos('../sim/spiral_rect_omegas0.33_alpha-7%s.sav' % _HIVRESSTR)
        vlosdiffgamma= galpy_simulations.vlos('../sim/spiral_rect_omegas0.33_gamma0.3%s.sav' % '')#_HIVRESSTR)
        vlosdiffomegas= galpy_simulations.vlos('../sim/spiral_rect_alpha-14%s.sav' % _HIVRESSTR)
        potscale= 0.85
        eres= 19
        xgrid= numpy.linspace((_RCXMIN-8.)/8.+_RCDX/8./2.,
                              (_RCXMAX-8.)/8.-_RCDX/8./2.,
                              eres)
        dx= (xgrid[1]-xgrid[0])*8.
        psdfid= bovy_psd.psd1d(vlosfid,dx,binsize=0.8)
        psdlpitch= bovy_psd.psd1d(vloslpitch,dx,binsize=0.8)
        psddiffgamma= bovy_psd.psd1d(vlosdiffgamma,dx,binsize=0.8)
        psddiffomegas= bovy_psd.psd1d(vlosdiffomegas,dx,binsize=0.8)
        ks= psdfid[0][1:-3]
        scale= 4.*numpy.pi*220.
        bovy_plot.bovy_print(fig_width=7.,axes_labelsize=20)
        line1= bovy_plot.bovy_plot(ks,potscale*scale*numpy.sqrt(psdfid[1][1:-3]),
                                   'k-',lw=2,
                                   semilogx=True,
                                   xlabel=r'$k\,(\mathrm{kpc}^{-1})$',
                                   ylabel=r'$\sqrt{P_k}\,(\mathrm{km\,s}^{-1})$',
                                   xrange=xrange,
                                   yrange=[0.,11.9],zorder=1)
        potscale= 0.25
        line2= bovy_plot.bovy_plot(ks,potscale*scale*numpy.sqrt(psdlpitch[1][1:-3]),
                                   'y-',lw=2.,zorder=1,
                                   overplot=True)
        potscale= 0.5
        line3= bovy_plot.bovy_plot(ks,potscale*scale*numpy.sqrt(psddiffgamma[1][1:-3]),
                                   'r-',lw=2.,zorder=1,
                                   overplot=True)
        line4= bovy_plot.bovy_plot(ks,4.*scale*numpy.sqrt(psddiffomegas[1][1:-3]),
                                   'b-',lw=2.,zorder=1,
                                   overplot=True)
        pyplot.annotate(r'$\mathrm{Spiral\ perturbation}$',
                        (0.5,1.08),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=20.)
        l1= pyplot.legend((line1[0],line2[0],line3[0],line4[0]),
                          (r'$\mathrm{Fiducial}$',
                           r'$\mathrm{Pitch\ angle} = 16^\circ$',
                           r'$\gamma = 17^\circ$',
                           r'$\Omega_s = 0.65\,\Omega_0$'),
                          loc='upper right',#bbox_to_anchor=(.91,.375),
                          numpoints=8,
                          prop={'size':16},
                          frameon=False)    
    #Also plot fiducial
    bovy_plot.bovy_plot(tks,
                        scale*numpy.sqrt(simpsd1d[1][1:-3]),
                        '-',color='0.65',lw=2.,overplot=True,zorder=0)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_psd_model(sys.argv[1],sys.argv[2])
