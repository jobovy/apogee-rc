import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import isodist
import rcmodel
from matplotlib.ticker import FuncFormatter, MultipleLocator
def G(g,gz):
    return g-0.1154-0.4175*gz-0.0497*gz**2.+0.0016*gz**3. #From Jordi et al. (2010)
def vifromgi(gi):
    if isinstance(gi,numpy.ndarray):
        out= 1.11*gi-0.52
        out[gi <= 2.1]= 0.675*gi[gi <= 2.1]+0.364
        return out
    else:
        if gi <= 2.1:
            return 0.675*gi+0.364
        else:
            return 1.11*gi-0.52
def plxerr(G,vi):
    z= numpy.amax([10.**(0.4*(12.-15.)),10.**(0.4*(G-15))])
    return 1./1000.*(9.3+658.1*z+4.568*z**2.)**0.5*(0.986+(1.-0.986)*vi)
def plot_gaiamag(plotfilename,plx=False):
    Zs= [0.0035,0.017,0.035]
    iso= isodist.PadovaIsochrone(type='sdss-2mass',parsec=True,
                                 Z=Zs)
    tages= [9.,9.3,9.5,9.7,9.85,10.,10.1]
    overplot= False
    bovy_plot.bovy_print()
    hs= [12.2,12.8,13.3,13.8]
    colors= ['r','y','g','b']
    for jj in range(len(hs)):
        h= hs[jj]
        for age in tages:
            for tZ in Zs:
                p= iso(age,tZ)
                jk= p['J']-p['Ks']
                indx= (jk < 0.8)*(jk > 0.5)\
                    *(tZ <= 0.06)\
                    *(tZ <= rcmodel.jkzcut(jk,upper=True))\
                    *(tZ >= rcmodel.jkzcut(jk))\
                    *(p['logg'] >= 1.8)\
                    *(p['logg'] <= 2.8)            
                tG= G(p['g'][indx],p['g'][indx]-p['z'][indx])-p['H'][indx]+h
                tdm= h+1.49
                tdist= 10.**(tdm/5.-2.) #kpc
                tplxerr= numpy.array([plxerr(tG[ii],
                                                   vifromgi(p['g'][indx][ii]-p['i'][indx][ii])) for ii in range(numpy.sum(indx))])
                if plx:
                    bovy_plot.bovy_plot(tG,
                                        tplxerr*tdist*100.,
                                        marker='o',color=colors[jj],mec='none',
                                        overplot=overplot,
                                        xlabel=r'$Gaia\ G$',
                                        ylabel=r'$\sigma_\pi/\pi\,(\%)$',
                                        xrange=[11.,21.],
                                        yrange=[0.,200.])
                else:
                    bovy_plot.bovy_plot(tG,
                                        4.74047*tdist*tplxerr*0.5,
                                        marker='o',color=colors[jj],mec='none',
                                        overplot=overplot,
                                        xlabel=r'$Gaia\ G$',
                                        ylabel=r'$\mathrm{transverse\ velocity\ uncertainty}\ (\mathrm{km\,s}^{-1})$',
                                        xrange=[11.,21.],
                                        yrange=[0.,8.])
                overplot= True
                #half a magnitude ah
                tvi= numpy.array([vifromgi(p['g'][indx][ii]-p['i'][indx][ii]) for ii in range(numpy.sum(indx))])
                ag= .5/0.18307*(0.8426-0.1187*tvi+0.0157*tvi**2.-0.0007*tvi**3.)
                tG= tG-0.5
                tplxerr= numpy.array([plxerr(tG[ii]+ag[ii],
                                                   vifromgi(p['g'][indx][ii]+.5/0.18307*1.20585-p['i'][indx][ii]-.5/0.18307*0.68319)) for ii in range(numpy.sum(indx))])
                tdm= h+1.49-0.5
                tdist= 10.**(tdm/5.-2.) #kpc
                #print numpy.median(ag)/.5, numpy.median(tvi), numpy.std(tvi)
                if plx:
                    bovy_plot.bovy_plot(tG+ag,
                                        tdist*tplxerr*100.,
                                        marker='o',color=colors[jj],mec='none',
                                        overplot=overplot)
                else:
                    bovy_plot.bovy_plot(tG+ag,
                                        4.74047*tdist*tplxerr*0.5,
                                        marker='o',color=colors[jj],mec='none',
                                        overplot=overplot)
                #Full magnitude of ah
                tvi= numpy.array([vifromgi(p['g'][indx][ii]-p['i'][indx][ii]) for ii in range(numpy.sum(indx))])
                ag= 1./0.18307*(0.8426-0.1187*tvi+0.0157*tvi**2.-0.0007*tvi**3.)
                tG= tG-0.5
                tplxerr= numpy.array([plxerr(tG[ii]+ag[ii],
                                                   vifromgi(p['g'][indx][ii]+1./0.18307*1.20585-p['i'][indx][ii]-1./0.18307*0.68319)) for ii in range(numpy.sum(indx))])
                tdm= h+1.49-1.
                tdist= 10.**(tdm/5.-2.) #kpc
                if plx:
                    bovy_plot.bovy_plot(tG+ag,
                                        tdist*tplxerr*100.,
                                        marker='o',color=colors[jj],mec='none',
                                        overplot=overplot)
                else:
                    bovy_plot.bovy_plot(tG+ag,
                                        4.74047*tdist*tplxerr*0.5,
                                        marker='o',color=colors[jj],mec='none',
                                        overplot=overplot)
                if True:
                    #1.5 magnitudes of ah
                    tvi= numpy.array([vifromgi(p['g'][indx][ii]-p['i'][indx][ii]) for ii in range(numpy.sum(indx))])
                    ag= 1.5/0.18307*(0.8426-0.1187*tvi+0.0157*tvi**2.-0.0007*tvi**3.)
                    tG= tG-.5
                    tplxerr= numpy.array([plxerr(tG[ii]+ag[ii],
                                                 vifromgi(p['g'][indx][ii]+1./0.18307*1.20585-p['i'][indx][ii]-1./0.18307*0.68319)) for ii in range(numpy.sum(indx))])
                    tdm= h+1.49-1.5
                    tdist= 10.**(tdm/5.-2.) #kpc
                    if plx:
                        bovy_plot.bovy_plot(tG+ag,
                                            tdist*tplxerr*100.,
                                            marker='o',color=colors[jj],mec='none',
                                            overplot=overplot)
                    else:
                        bovy_plot.bovy_plot(tG+ag,
                                            4.74047*tdist*tplxerr*0.5,
                                            marker='o',color=colors[jj],mec='none',
                                            overplot=overplot)
    bovy_plot.bovy_plot([20.,20.],[0.,300.],'k:',overplot=True)
    #Work on legend
    dot2= bovy_plot.bovy_plot([-10.],[-10.],'o',color=colors[0],overplot=True,mec='none')
    dot3= bovy_plot.bovy_plot([-10.],[-10.],'o',color=colors[1],overplot=True,mec='none')
    dot4= bovy_plot.bovy_plot([-10.],[-10.],'o',color=colors[2],overplot=True,mec='none')
    dot5= bovy_plot.bovy_plot([-10.],[-10.],'o',color=colors[3],overplot=True,mec='none')
    pyplot.legend((dot2[0],dot3[0],dot4[0],dot5[0]),
                  (r'$H = 12.2$',
                   r'$H = 12.8 \ \ \ \ \mathrm{each\ incl.}$',
                   r'$H = 13.3 \ \ \ A_H= 0, 0.5,$',
                   r'$H = 13.8 \ \ \ \ \ \ \ \ \ \ \ \ \, 1, 1.5$'),
                  loc='upper left',#bbox_to_anchor=(.91,.375),
                  numpoints=1,
                  prop={'size':16},
                  frameon=False)
    if plx:
        bovy_plot.bovy_plot([5.,25.],[10.,10.],'k--',lw=2.,overplot=True)
        bovy_plot.bovy_text(11.25,2.,r'$\mathrm{spectro-photometric\ RC\ precision}$',
                            size=16.)
        from matplotlib.patches import FancyArrowPatch
        ax=pyplot.gca()
        ax.add_patch(FancyArrowPatch((20.5,10.),(20.5,-1.),
                                     arrowstyle='->',mutation_scale=15,
                                     fill=True,
                                     lw=1.25,color='k'))
    #Add twin x axis w/ distsance
    ax= pyplot.gca()
    def my_formatter(x, pos):
        """distance in kpc for zero reddening"""
        xs= 10.**((x-0.71)/5.-2.)
        xs2= 10.**((x-0.71-4.)/5.-2.)
        if xs2 > 1:
            return r'$%.0f/%.0f$' % (xs,xs2)
        else:
            return r'$%.0f/%.1f$' % (xs,xs2)
    ax2= pyplot.twiny()
    major_formatter = FuncFormatter(my_formatter)
    ax2.xaxis.set_major_formatter(major_formatter)
    xstep= ax.xaxis.get_majorticklocs()
    xstep= xstep[1]-xstep[0]
    ax2.xaxis.set_minor_locator(MultipleLocator(xstep/5.))
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    xmin, xmax= ax.xaxis.get_view_interval()
    ax2.xaxis.set_view_interval(xmin,xmax,ignore=True)
    ax2.set_xlabel('$\mathrm{distance\ for}\ A_H = 0/1\,(\mathrm{kpc})$')
    ystep= ax.yaxis.get_majorticklocs()
    ystep= ystep[1]-ystep[0]
    ax2.yaxis.set_minor_locator(MultipleLocator(ystep/5.))
    #Save
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    plot_gaiamag(sys.argv[1],plx=len(sys.argv) > 2)
