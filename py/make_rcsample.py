import sys
import numpy
import fitsio
import esutil
from galpy.util import bovy_coords
import isodist
import apogee.tools.read as apread
import rcmodel
_ADDHAYDENDIST= True
def make_rcsample(savefilename):
    #Read the base-sample
    data= apread.allStar(adddist=_ADDHAYDENDIST)
    #Select red-clump stars
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'])
    logg= data['LOGG']
    indx= (jk < 0.75)*(jk > 0.5)\
        *(z <= rcmodel.jkzcut(jk,upper=True))\
        *(z >= rcmodel.jkzcut(jk))\
        *(logg >= 1.8)\
        *(logg <= 2.8)            
    data= data[indx]
    #Add distances
    data= esutil.numpy_util.add_fields(data,[('RC_DIST', float),
                                             ('RC_DM', float),
                                             ('RC_GALR', float),
                                             ('RC_GALPHI', float),
                                             ('RC_GALZ', float)])
    rcd= rcmodel.rcdist('../../rcdist-apogee/data/rcmodel_mode_jkz_ks.sav')
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'])
    data['RC_DIST']= rcd(jk,z,appmag=data['K0'])
    data['RC_DM']= 5.*numpy.log10(data['RC_DIST'])+10.
    XYZ= bovy_coords.lbd_to_XYZ(data['GLON'],
                                data['GLAT'],
                                data['RC_DIST'],
                                degree=True)
    R,phi,Z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                          XYZ[:,1],
                                          XYZ[:,2],
                                          Xsun=8.)
    data['RC_GALR']= R
    data['RC_GALPHI']= phi
    data['RC_GALZ']= Z
    #Save
    fitsio.write(savefilename,data,clobber=True)
    return None

if __name__ == '__main__':
    make_rcsample(sys.argv[1])
