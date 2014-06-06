#Read the GCS data and hack it into a similar format as the RC data such that we can re-use RC code
import numpy
from astropy.coordinates import SkyCoord
import esutil
from galpy.util import bovy_coords
import readGCS
def hackGCS(cutbinaries=True,cutvalidfeh=True,cutfeh=True,
            casagrande=False,irfm=True):
    data= readGCS.readGCS(cutbinaries=cutbinaries,
                          cutvalidfeh=cutvalidfeh,
                          cutfeh=cutfeh,
                          casagrande=casagrande,irfm=irfm)
    #To allow for XY pixelization, we will hack these
    data= esutil.numpy_util.add_fields(data,[('RC_GALR', float),
                                             ('RC_GALPHI', float),
                                             ('RA',float),
                                             ('DEC',float),
                                             ('GLON',float),
                                             ('GLAT',float)])
    ras= numpy.empty(len(data))
    decs= numpy.empty(len(data))
    for ii in range(len(data)):
        c = SkyCoord(str(data['RAh'][ii])+'h'+str(data['RAm'][ii])+'m'
                     +str(data['RAs'][ii])+'s',
                     str(data['DE-'][ii])+str(data['DEd'][ii])+'d'
                     +str(data['DEm'][ii])+'m'+str(data['DEs'][ii])+'s',
                     'icrs')
        ras[ii]= c.ra.degree
        decs[ii]= c.dec.degree
    llbb= bovy_coords.radec_to_lb(ras,decs,degree=True)
    data['RA']= ras
    data['DEC']= decs
    data['GLON']= llbb[:,0]
    data['GLAT']= llbb[:,1]
    data['RC_GALR']= data['Dist']/1000.*numpy.cos(llbb[:,1]/180.*numpy.pi)
    data['RC_GALPHI']= llbb[:,0]/180.*numpy.pi
    return data
