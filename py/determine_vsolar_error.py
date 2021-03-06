#Determine the error on the solar motion inferred by minimizing the large-scale power
import sys
import numpy
from galpy.df import dehnendf
from galpy.potential import SteadyLogSpiralPotential
import apogee.tools.read as apread
import pixelize_sample
from plot_2dkinematics import dvlosgal
from plot_psd import _ADDLLOGGCUT, \
    _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
import bovy_psd
from determine_vsolar import large_scale_power
import galpy_simulations
def determine_vsolar_error(spiral=False):
    #Read the APOGEE-RC data and pixelate it
    #APOGEE-RC
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    data= data[indx]
    #Set up the Dehnen DF that we will sample from
    dfc= dehnendf(beta=0.,profileParams=(3./8.,33.3,31.4/220.))
    trueVsolar= 24./220.
    ntrials= 100
    mockvsolars= []
    for ii in range(ntrials):
        print "Working on %i / %i" % (ii+1,ntrials)
        mockvsolars.append(determine_vsolar_mock(data,dfc,trueVsolar,spiral))
    mockvsolars= numpy.array(mockvsolars)
    print mockvsolars
    print numpy.mean(mockvsolars-24.), numpy.std(mockvsolars)
    print numpy.sum(numpy.fabs(mockvsolars-24.) > 1.)/float(ntrials)
    return None

def determine_vsolar_mock(data,dfc,trueVsolar,spiral=False):
    #At the position of each real data point, generate a mock velocity
    data= create_mock_sample(data,dfc,trueVsolar,spiral=spiral)
    #Get velocity field
    pix= pixelize_sample.pixelXY(data,
                                 xmin=_RCXMIN,xmax=_RCXMAX,
                                 ymin=_RCYMIN,ymax=_RCYMAX,
                                 dx=_RCDX,dy=_RCDX)
    vsolars= numpy.linspace(0.,40.,51)
    lpower= large_scale_power(pix,vsolars,vc=220.,dx=_RCDX,beta=0.)
    p= numpy.polyfit(vsolars,lpower,2)
    minvsolar= -0.5*p[1]/p[0]
    return minvsolar


def create_mock_sample(data,dfc,trueVsolar,spiral=False):
    if spiral:
        #Read spiral vlos field
        spvlos= galpy_simulations.vlos('../sim/bar_rect_alpha0.015_hivres.sav')
        potscale= 0.85
        xgrid= numpy.linspace(_RCXMIN,_RCXMAX-_RCDX,19)
        ygrid= numpy.linspace(_RCYMIN,_RCYMAX-_RCDX,19)
    ndata= len(data)
    mockvel= numpy.empty(ndata)
    dphil= data['RC_GALPHI']+data['GLON']/180.*numpy.pi
    cosdphil= numpy.cos(dphil)
    sindphil= numpy.sin(dphil)
    cosl= numpy.cos(data['GLON']/180.*numpy.pi)
    sinl= numpy.sin(data['GLON']/180.*numpy.pi)
    for ii in range(ndata):
        vrvt= dfc.sampleVRVT(data['RC_GALR'][ii]/8.,n=1.)
        mockvel[ii]= -vrvt[0,0]*cosdphil[ii]\
            +vrvt[0,1]*sindphil[ii]\
            -10.5*cosl[ii]/220.\
            -(1.+trueVsolar)*sinl[ii]
        if spiral:
            #Find which pixel this sample is in
            pixIndxX= numpy.argmin(numpy.fabs(data['RC_GALR'][ii]\
                                                  *numpy.cos(data['RC_GALPHI'][ii])\
                                                  -xgrid))
            pixIndxY= numpy.argmin(numpy.fabs(data['RC_GALR'][ii]\
                                                  *numpy.sin(data['RC_GALPHI'][ii])\
                                                  -ygrid))
            spiraldev= spvlos[pixIndxX,pixIndxY]
            mockvel[ii]+= spiraldev*potscale
    data['VHELIO_AVG']= mockvel*220.
    return data

if __name__ == '__main__':
    if len(sys.argv) > 1:
        determine_vsolar_error(spiral=True)
    else:
        determine_vsolar_error()
