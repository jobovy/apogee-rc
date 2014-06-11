#Plot the vlos velocity field of the RC, with and without subtracting the Solar motion
import sys
import numpy
from galpy.util import bovy_plot
import apogee.tools.read as apread
import pixelize_sample
from plot_2dkinematics import vlosgal
from plot_psd import _ADDLLOGGCUT, \
    _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
def plot_rckinematics(plotfilename,subsun=False):
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
                                   xmin=_RCXMIN-1.5,xmax=_RCXMAX+1.5,
                                   ymin=_RCYMIN-1.5,ymax=_RCYMAX+1.5,
                                   dx=_RCDX,dy=_RCDX)
    bovy_plot.bovy_print()
    if subsun:
        vmin, vmax= 0., 250.
        pixrc.plot(lambda x: vlosgal(x),
                   func=lambda x: numpy.fabs(numpy.median(x)),
                   zlabel=r'$|\mathrm{median}\ V^{\mathrm{GC}}_{\mathrm{los}}|\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax)
    else:
        vmin, vmax= -75., 75.
        pixrc.plot('VHELIO_AVG',
                   zlabel=r'$\mathrm{median}\ V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$',
                   vmin=vmin,vmax=vmax)
        bovy_plot.bovy_text(r'$\mathrm{typical\ uncertainty\!:}\ 3\,\mathrm{km\,s}^{-1}$',
                            bottom_left=True,size=18.)
        bovy_plot.bovy_text(r'$|Z| < 250\,\mathrm{pc}$',top_right=True,size=18.)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 2:
        plot_rckinematics(sys.argv[1],subsun=True)
    else:
        plot_rckinematics(sys.argv[1],subsun=False)
