import sys
import os, os.path
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter
import fitsio
import apogee.tools.read as apread
import pixelize_sample
import hackGCS
from plot_2dkinematics import dvlosgal
from plot_psd import _ADDLLOGGCUT, \
    _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX, \
    _RAVEXMIN, _RAVEXMAX, _RAVEYMIN, _RAVEYMAX, _RAVEDX, \
    _GCSXMIN, _GCSXMAX, _GCSYMIN, _GCSYMAX, _GCSDX
def plot_velocityfield(plotfilename,err=False):
    #Read the 3 data files and pixelate them
    #APOGEE-RC
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    data= data[indx]
    #Get residuals
    pixrc= pixelize_sample.pixelXY(data,
                                   xmin=_RCXMIN,xmax=_RCXMAX,
                                   ymin=_RCYMIN,ymax=_RCYMAX,
                                   dx=_RCDX,dy=_RCDX)
    #RAVE
    data= fitsio.read(os.path.join(os.getenv('DATADIR'),'rave','ravedr4_rc.fits'))
    pixrave= pixelize_sample.pixelXY(data,
                                     xmin=_RAVEXMIN,xmax=_RAVEXMAX,
                                     ymin=_RAVEYMIN,ymax=_RAVEYMAX,
                                     dx=_RAVEDX,dy=_RAVEDX)
    #GCS
    data= hackGCS.hackGCS()
    pixgcs= pixelize_sample.pixelXY(data,
                                    xmin=_GCSXMIN,xmax=_GCSXMAX,
                                    ymin=_GCSYMIN,ymax=_GCSYMAX,
                                    dx=_GCSDX,dy=_GCSDX)
    #Set up 3 axes
    bovy_plot.bovy_print(fig_width=8.,axes_labelsize=14)
    axdx= 1./3.
    #APOGEE-RC
    tdy= (_RCYMAX-_RCYMIN)/(_RCXMAX-_RCXMIN)*axdx
    rcAxes= pyplot.axes([0.1,(1.-tdy)/2.,axdx,tdy])
    pyplot.sca(rcAxes)
    if err:
        img= pixrc.plot(lambda x: dvlosgal(x),vmin=0.,vmax=8.,overplot=True,
                        func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                        colorbar=False)
    else:
        img= pixrc.plot(lambda x: dvlosgal(x),vmin=-16.,vmax=16.,overplot=True,
                        colorbar=False)
    pyplot.axis([pixrc.xmin,pixrc.xmax,pixrc.ymin,pixrc.ymax])
    bovy_plot._add_ticks()
    bovy_plot._add_axislabels(r'$ $',
                              r'$Y_{\mathrm{GC}}\,(\mathrm{kpc})$')
    #Plot RAVE box
    plot_box(_RAVEXMIN,_RAVEXMAX,_RAVEYMIN,_RAVEYMAX,ls='--')
    line= pyplot.plot([_RAVEXMIN,14.9],[_RAVEYMAX,_RCYMAX],'k:')
    line[0].set_clip_on(False)
    line= pyplot.plot([_RAVEXMIN,14.9],[_RAVEYMIN,_RCYMIN],'k:')
    line[0].set_clip_on(False)
    #Colorbar
    cbaxes = pyplot.axes([0.1625,(1.-tdy)/2.+tdy+0.04,axdx-0.125,0.02])
    CB1= pyplot.colorbar(img,orientation='horizontal',
                         cax=cbaxes,ticks=[-16.,-8.,0.,8.,16.])
    CB1.set_label(r'$\mathrm{APOGEE-RC}$',labelpad=-35,fontsize=14.)
    #RAVE
    tdy= (_RAVEYMAX-_RAVEYMIN)/(_RAVEXMAX-_RAVEXMIN)*axdx
    raveAxes= pyplot.axes([0.05+axdx,(1.-tdy)/2.,axdx,tdy])
    pyplot.sca(raveAxes)
    if err:
        img= pixrave.plot(lambda x: dvlosgal(x,vtsun=230.),vmin=-0.,vmax=8.,
                          func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                          overplot=True,colorbar=False)
    else:
        img= pixrave.plot(lambda x: dvlosgal(x,vtsun=230.),vmin=-7.,vmax=7.,
                          overplot=True,colorbar=False)
    pyplot.axis([pixrave.xmin,pixrave.xmax,pixrave.ymin,pixrave.ymax])
    bovy_plot._add_ticks()
    bovy_plot._add_axislabels(r'$X_{\mathrm{GC}}\,(\mathrm{kpc})$',
                              r'$ $')
    #Plot GCS box
    plot_box(_GCSXMIN+8.,_GCSXMAX+8.,_GCSYMIN,_GCSYMAX,ls='--')
    line= pyplot.plot([_GCSXMIN+8.,9.48],[_GCSYMAX,_RAVEYMAX],'k:')
    line[0].set_clip_on(False)
    line= pyplot.plot([_GCSXMIN+8.,9.48],[_GCSYMIN,_RAVEYMIN],'k:')
    line[0].set_clip_on(False)
    #Colorbar
    cbaxes = pyplot.axes([0.1125+axdx,(1.-tdy)/2.+tdy+0.04,axdx-0.125,0.02])
    CB1= pyplot.colorbar(img,orientation='horizontal',
                    cax=cbaxes,ticks=[-6.,-3.,0.,3.,6.])
    CB1.set_label(r'$\mathrm{RAVE-RC}$',labelpad=-35,fontsize=14.)
    #GCS
    tdy= (_GCSYMAX-_GCSYMIN)/(_GCSXMAX-_GCSXMIN)*axdx
    gcsAxes= pyplot.axes([0.0+2.*axdx,(1.-tdy)/2.,axdx,tdy])
    pyplot.sca(gcsAxes)
    if err:
        img= pixgcs.plot(lambda x: x['VVel']-numpy.median(data['VVel']),
                         func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                         vmin=0.,vmax=8.,overplot=True,colorbar=False)
    else:
        img= pixgcs.plot(lambda x: x['VVel']-numpy.median(data['VVel']),
                         vmin=-4.,vmax=4.,overplot=True,colorbar=False)
    def my_formatter(x, pos):
        return r'$%g$' % (8.+x)
    major_formatter = FuncFormatter(my_formatter)
    gcsAxes.xaxis.set_major_formatter(major_formatter)
    pyplot.axis([pixgcs.xmin,pixgcs.xmax,pixgcs.ymin,pixgcs.ymax])
    gcsAxes.xaxis.set_ticks([-0.05,0.,0.05])
    gcsAxes.yaxis.set_ticks([-0.05,0.,0.05])
    bovy_plot._add_ticks()
    #Colorbar
    cbaxes = pyplot.axes([0.0625+2.*axdx,(1.-tdy)/2.+tdy+0.04,axdx-0.125,0.02])
    CB1= pyplot.colorbar(img,orientation='horizontal',
                    cax=cbaxes,ticks=[-4.,-2.,0.,2.,4.])
    CB1.set_label(r'$\mathrm{GCS}$',labelpad=-35,fontsize=14.)
    bovy_plot.bovy_end_print(plotfilename)
    return None

def plot_box(xmin,xmax,ymin,ymax,ls='-',color='k'):
    pyplot.plot([xmin,xmin],[ymin,ymax],ls=ls,color=color)
    pyplot.plot([xmax,xmax],[ymin,ymax],ls=ls,color=color)
    pyplot.plot([xmin,xmax],[ymin,ymin],ls=ls,color=color)
    pyplot.plot([xmin,xmax],[ymax,ymax],ls=ls,color=color)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 2:
        plot_velocityfield(sys.argv[1],err=True)
    else:
        plot_velocityfield(sys.argv[1],err=False)
