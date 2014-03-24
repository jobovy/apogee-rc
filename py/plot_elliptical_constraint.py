import sys
import math
import numpy
from galpy.util import bovy_plot
def plot_elliptical_constraint(plotfilename):
    cs= numpy.linspace(-0.25,0.25,201)
    ss= numpy.linspace(-0.25,0.25,201)
    cons= numpy.zeros((len(cs),len(ss),100))
    cnt= 0
    for ii in range(len(cs)):
        for jj in range(len(ss)):
            twophio= numpy.sqrt(cs[ii]**2.+ss[jj]**2.)
            phib= math.atan2(ss[jj],cs[ii])/2.
            q= 1-twophio
            for zz in range(7):
                for ww in range(zz,7):
                    cnt+= 1
                    cons[ii,jj,ww+zz*7]= q**2.*(1./numpy.sqrt(q**2.*numpy.cos((25.-zz*5.)*numpy.pi/180.-phib)**2.+numpy.sin((25.-zz*5.)*numpy.pi/180.-phib)**2.)-1./numpy.sqrt(q**2.*numpy.cos((20.-5.*ww)*numpy.pi/180.-phib)**2.+numpy.sin((20.-ww*5.)*numpy.pi/180.-phib)**2.))**2.
                    #cons[ii,jj,zz]= q**2.*(1./numpy.sqrt(q**2.*numpy.cos((25.*numpy.pi/180.-phib)**2.+numpy.sin(25.*numpy.pi/180.-phib)**2.)-1./numpy.sqrt(q**2.*numpy.cos((20.-5.*zz)*numpy.pi/180.-phib)**2.+numpy.sin((20.-zz*5.)*numpy.pi/180.-phib)**2.))**2.
            if ii == 0 and jj == 0: print cnt
    plotthis= numpy.all(cons < 1./40.**2.,axis=2)*0.65
    #Plot
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(plotthis.T,
                          origin='lower',
                          cmap='gist_yarg',
                          xrange=[-0.25,0.25],
                          yrange=[-0.25,0.25],
                          xlabel=r'$c_\Psi = \varepsilon_\Psi\,\cos 2 \phi_b$',
                          ylabel=r'$s_\Psi = \varepsilon_\Psi\,\sin 2 \phi_b$',
                          vmax=1.)
    bovy_plot.bovy_plot([-0.15,-0.12],[0.15,0.02],'k-',zorder=1,
                        overplot=True)
    bovy_plot.bovy_plot([-0.06363636363636363],[0.],'wx',overplot=True,ms=10.,
                        mew=2.,zorder=2.)
    bovy_plot.bovy_plot([-0.06363636363636363,-0.02],[0.,-0.14],'k-',zorder=1,
                        overplot=True)
    bovy_plot.bovy_text(-0.08,-0.2,r'$\mathrm{model\ with}$'+'\n'+r'$\mathrm{local}\ V_c - \mathrm{global}\ V_c = 14\,\mathrm{km\,s}^{-1}$',
                         size=13.)
    bovy_plot.bovy_text(r'$\mathrm{Elliptical\ disk\ models\ with}$'+'\n'+r'$\forall i \in [0,7] \cap \mathbb{Z}:  \forall j \in [i,7] \cap \mathbb{Z}:$'+'\n'+r'$|\Delta R(\phi=25^\circ - 5^\circ\,i,\phi=20^\circ-5^\circ\,j)| < 200\,\mathrm{pc}$',
                        size=13.,
                        top_left=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_elliptical_constraint(sys.argv[1])
