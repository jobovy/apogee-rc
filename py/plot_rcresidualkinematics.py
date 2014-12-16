#Plot the vlos velocity field of the RC, with and without subtracting the Solar motion
import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
import pixelize_sample
from plot_2dkinematics import dvlosgal
from plot_psd import _ADDLLOGGCUT, \
    _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
def plot_rcresidualkinematics(plotfilename,sbd10=False,vc=220.):
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
    bovy_plot.bovy_print()
    vmin, vmax= -16., 16.
    if sbd10:
        pixrc.plot(lambda x: dvlosgal(x,vc=vc,vtsun=vc+12.),
                   zlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los,rot}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax)
        pyplot.annotate(r'$V_c = %i\,\mathrm{km\,s}^{-1}, V_{\odot-c}= 12\,\mathrm{km\,s}^{-1}$' % vc,
                        (0.5,1.1),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=18.)
    else: 
        resv= pixrc.plot(lambda x: dvlosgal(x,vc=vc,vtsun=vc+24.),
                   zlabel=r'$\mathrm{median}\ \Delta V_{\mathrm{los,rot}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax,returnz=True)
        notNan= True-numpy.isnan(resv)
        print numpy.sum(notNan)
        print numpy.median(resv[notNan])
        print 1.4826*numpy.median(numpy.fabs(resv[notNan]-numpy.median(resv[notNan])))
        print numpy.mean(resv[notNan])
        print numpy.std(resv[notNan])
        pyplot.annotate(r'$V_c = %i\,\mathrm{km\,s}^{-1}, V_{\odot-c}= 24\,\mathrm{km\,s}^{-1}$' % vc,
                        (0.5,1.1),xycoords='axes fraction',
                        horizontalalignment='center',
                        verticalalignment='top',size=18.)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 3:
        plot_rcresidualkinematics(sys.argv[1],sbd10=True,vc=float(sys.argv[2]))
    elif len(sys.argv) > 2:
        plot_rcresidualkinematics(sys.argv[1],sbd10=False,vc=float(sys.argv[2]))
    else:
        plot_rcresidualkinematics(sys.argv[1],sbd10=False,vc=220.)
