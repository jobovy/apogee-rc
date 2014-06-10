#Module to read the RAVE DR4 data
import os, os.path
import numpy
import asciitable
import esutil
from galpy.util import bovy_coords
import isodist
import rcmodel
_DATAFILE= os.path.join(os.getenv('DATADIR'),'rave','ravedr4.dat')
_DATAREADME= os.path.join(os.getenv('DATADIR'),'rave','ReadMe')
def readRave():
    data= asciitable.read(_DATAFILE,
                          readme=_DATAREADME,
                          Reader=asciitable.cds.Cds,
                          guess=False,
                          fill_values=[('', '-999')])
    return data

def raveRC():
    data= readRave()
    jk= data['Jmag2']-data['Kmag2']-0.17*numpy.exp(data['Av'])
    z= isodist.FEH2Z(data['[M/H]K'],zsolar=0.017)
    logg= data['loggK']
    indx= (jk < 0.8)*(jk >= 0.5)\
        *(z <= 0.06)\
        *(z <= rcmodel.jkzcut(jk,upper=True))\
        *(z >= rcmodel.jkzcut(jk))\
        *(logg >= rcmodel.loggteffcut(data['TeffK'],z,upper=False))\
        *(logg <= rcmodel.loggteffcut(data['TeffK'],z,upper=True))
    data= data[indx]
    #To allow for XY pixelization
    data= esutil.numpy_util.add_fields(data,[('RC_GALR', float),
                                             ('RC_GALPHI', float),
                                             ('RC_GALZ', float),
                                             ('VHELIO_AVG', float)])
    XYZ= bovy_coords.lbd_to_XYZ(data['GLON'],
                                data['GLAT'],
                                data['Dist'],
                                degree=True)
    R,phi,Z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                          XYZ[:,1],
                                          XYZ[:,2],
                                          Xsun=8.,Zsun=0.025)
    data['RC_GALR']= R
    data['RC_GALPHI']= phi
    data['RC_GALZ']= Z   
    data['VHELIO_AVG']= data['HRV']
    return data
