#Plot the residual vlos velocity field of the RC and that of the model
import os, os.path
import sys
import copy
import pickle
import numpy
from galpy.util import bovy_plot, bovy_coords, save_pickles
from matplotlib import pyplot
from matplotlib import transforms
import apogee.tools.read as apread
import pixelize_sample
from plot_2dkinematics import dvlosgal
import galpy_simulations
from plot_psd import _ADDLLOGGCUT, \
    _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
def plot_rckinematics(plotfilename,subsun=False):
    #Set up 3 axes
    bovy_plot.bovy_print(fig_width=8.,axes_labelsize=14)
    axdx= 1./3.
    #APOGEE-RC observations
    tdy= (_RCYMAX-_RCYMIN+4.5)/(_RCXMAX-_RCXMIN+4.5)*axdx
    obsAxes= pyplot.axes([0.1,(1.-tdy)/2.,axdx,tdy])
    pyplot.sca(obsAxes)
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    data= data[indx]
    #Get velocity field
    pixrc= pixelize_sample.pixelXY(data,
                                   xmin=_RCXMIN-2.25,xmax=_RCXMAX+2.25,
                                   ymin=_RCYMIN-2.25,ymax=_RCYMAX+2.25,
                                   dx=_RCDX,dy=_RCDX)
    vmin, vmax= -16.,16.
    img= pixrc.plot(lambda x: dvlosgal(x,vtsun=220.+24.),
                    vmin=vmin,vmax=vmax,overplot=True,
                    colorbar=False)
    resv= pixrc.plot(lambda x: dvlosgal(x,vtsun=220.+24.),
                     justcalc=True,returnz=True) #for later
    pyplot.annotate(r'$\mathrm{APOGEE\!-\!RC\ data}$',
                    (0.5,1.09),xycoords='axes fraction',
                    horizontalalignment='center',
                    verticalalignment='top',size=10.)
    pyplot.axis([pixrc.xmin,pixrc.xmax,pixrc.ymin,pixrc.ymax])
    bovy_plot._add_ticks()
    bovy_plot._add_axislabels(r'$X_{\mathrm{GC}}\,(\mathrm{kpc})$',
                              r'$Y_{\mathrm{GC}}\,(\mathrm{kpc})$')
    #Colorbar
    cbaxes = pyplot.axes([0.1+axdx/2.,(1.-tdy)/2.+tdy+0.065,2.*axdx-0.195,0.02])
    CB1= pyplot.colorbar(img,orientation='horizontal',
                         cax=cbaxes)#,ticks=[-16.,-8.,0.,8.,16.])
    CB1.set_label(r'$\mathrm{median}\ \Delta V_{\mathrm{los,rot}}\,(\mathrm{km\,s}^{-1})$',labelpad=-35,fontsize=14.)
    #Now calculate the expected field
    expec_vlos= galpy_simulations.vlos_altrect('../sim/bar_altrect_alpha0.015_hivres.sav')*220.
    modelAxes= pyplot.axes([0.03+axdx,(1.-tdy)/2.,axdx,tdy])
    pyplot.sca(modelAxes)
    xlabel=r'$X_{\mathrm{GC}}\,(\mathrm{kpc})$'
    ylabel=r'$Y_{\mathrm{GC}}\,(\mathrm{kpc})$'
    indx= True-numpy.isnan(resv)
    plotthis= copy.copy(expec_vlos)
    plotthis[numpy.isnan(resv)]= numpy.nan #turn these off
    bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='jet',
                          interpolation='nearest',
                          xlabel=xlabel,ylabel=ylabel,
                          xrange=[_RCXMIN-2.25,_RCXMAX+2.25],
                          yrange=[_RCYMIN-2.25,_RCYMAX+2.25],
                          contours=False,
                          vmin=vmin,vmax=vmax,overplot=True,zorder=3)
    if True:
       #Now plot the pixels outside the APOGEE data set
        plotthis= copy.copy(expec_vlos)
        plotthis[True-numpy.isnan(resv)]= numpy.nan #turn these off
        bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='jet',
                              interpolation='nearest',
                              alpha=0.3,
                              xrange=[_RCXMIN-2.25,_RCXMAX+2.25],
                              yrange=[_RCYMIN-2.25,_RCYMAX+2.25],
                              contours=False,
                              vmin=vmin,vmax=vmax,overplot=True,
                              zorder=0)
    pyplot.annotate(r'$\mathrm{Favored\ bar\ model}$',
                    (1.02,1.09),xycoords='axes fraction',
                    horizontalalignment='center',
                    verticalalignment='top',size=10.,zorder=3)
    pyplot.axis([_RCXMIN-2.25,_RCXMAX+2.25,_RCYMIN-2.25,_RCYMAX+2.25])
    bovy_plot._add_ticks()
    bovy_plot._add_axislabels(xlabel,r'$ $')
    #Finally, add a polar plot of the whole disk
    res= 51
    rmin, rmax= 0.2, 2.4
    xgrid= numpy.linspace(0.,2.*numpy.pi*(1.-1./res/2.),
                          2.*res)
    ygrid= numpy.linspace(rmin,rmax,res)
    expec_vlos= galpy_simulations.vlos_polar('../sim/bar_polar_alpha0.015_hivres.sav')*220.
    plotxgrid= numpy.linspace(xgrid[0]-(xgrid[1]-xgrid[0])/2.,
                              xgrid[-1]+(xgrid[1]-xgrid[0])/2.,
                              len(xgrid)+1)
    plotygrid= numpy.linspace(ygrid[0]-(ygrid[1]-ygrid[0])/2.,
                           ygrid[-1]+(ygrid[1]-ygrid[0])/2.,
                           len(ygrid)+1)
    fullmodelAxes= pyplot.axes([-0.05+2.*axdx,(1.-tdy)/2.,axdx,tdy],polar=True)
    ax= fullmodelAxes
    pyplot.sca(fullmodelAxes)
    out= ax.pcolor(plotxgrid,plotygrid,expec_vlos.T,cmap='jet',
                   vmin=vmin,vmax=vmax,clip_on=False)
    from matplotlib.patches import FancyArrowPatch
    arr= FancyArrowPatch(posA=(numpy.pi+0.1,1.8),
                         posB=(3*numpy.pi/2.+0.1,1.8),
                         arrowstyle='->', 
                         connectionstyle='arc3,rad=%4.2f' % (numpy.pi/8.-0.05),
                         shrinkA=2.0, shrinkB=2.0, mutation_scale=20.0, 
                         mutation_aspect=None,fc='k')
    ax.add_patch(arr)
    bovy_plot.bovy_text(numpy.pi+0.17,1.7,r'$\mathrm{Galactic\ rotation}$',
                        rotation=-30.,size=9.)
    radii= numpy.array([0.5,1.,1.5,2.,2.5])
    labels= []
    for r in radii:
        ax.plot(numpy.linspace(0.,2.*numpy.pi,501,),
                numpy.zeros(501)+r,ls='-',color='0.65',zorder=1,lw=0.5)
        labels.append(r'$%i$' % int(r*8.))
    pyplot.rgrids(radii,labels=labels,angle=147.5)
    thetaticks = numpy.arange(0,360,45)
    # set ticklabels location at x times the axes' radius
    ax.set_thetagrids(thetaticks,frac=1.16)
    bovy_plot.bovy_text(3.*numpy.pi/4.+0.06,2.095,r'$\mathrm{kpc}$',size=10.)
    pyplot.ylim(0.,2.8)
    #Plot the box
    xs= numpy.linspace(_RCXMIN-2.25,_RCXMAX+2.25,101)
    ys= numpy.ones(101)*(_RCYMIN-2.25)
    rs= numpy.sqrt(xs**2.+ys**2.)/8.
    phis= numpy.arctan2(ys,xs)    
    ax.plot(phis,rs,'--',lw=1.25,color='k')
    #Plot the box
    xs= numpy.linspace(_RCXMIN-2.25,_RCXMAX+2.25,101)
    ys= numpy.ones(101)*(_RCYMAX+2.25)
    rs= numpy.sqrt(xs**2.+ys**2.)/8.
    phis= numpy.arctan2(ys,xs)    
    ax.plot(phis,rs,'--',lw=1.25,color='k')
    #Plot the box
    ys= numpy.linspace(_RCYMIN-2.25,_RCYMAX+2.25,101)
    xs= numpy.ones(101)*(_RCXMIN-2.25)
    rs= numpy.sqrt(xs**2.+ys**2.)/8.
    phis= numpy.arctan2(ys,xs)    
    ax.plot(phis,rs,'--',lw=1.25,color='k')
    #Plot the box
    ys= numpy.linspace(_RCYMIN-2.25,_RCYMAX+2.25,101)
    xs= numpy.ones(101)*(_RCXMAX+2.25)
    rs= numpy.sqrt(xs**2.+ys**2.)/8.
    phis= numpy.arctan2(ys,xs)    
    ax.plot(phis,rs,'--',lw=1.25,color='k')
    #Plot the connectors on the modelAxes
    xlow=-4.*8.
    ylow= 2.77*8.
    xs= numpy.linspace(xlow,(_RCXMAX+2.25),101)
    ys= (ylow-(_RCYMAX+2.25))/(xlow-(_RCXMAX+2.25))*(xs-xlow)+ylow
    rs= numpy.sqrt(xs**2.+ys**2.)/8.
    phis= numpy.arctan2(ys,xs)    
    line= ax.plot(phis,rs,':',lw=1.,color='k')
    line[0].set_clip_on(False)
    xlow=-4.*8.
    ylow= -2.77*8.
    xs= numpy.linspace(xlow,(_RCXMAX+2.25),101)
    ys= (ylow-(_RCYMIN-2.25))/(xlow-(_RCXMAX+2.25))*(xs-xlow)+ylow
    rs= numpy.sqrt(xs**2.+ys**2.)/8.
    phis= numpy.arctan2(ys,xs)    
    line= ax.plot(phis,rs,':',lw=1.,color='k')
    line[0].set_clip_on(False)
    bovy_plot.bovy_end_print(plotfilename,dpi=300)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 2:
        plot_rckinematics(sys.argv[1],subsun=True)
    else:
        plot_rckinematics(sys.argv[1],subsun=False)
