import esutil
from galpy.util import bovy_coords
import apogee.tools.read as apread
def readAndHackHoltz():
    alldata= apread.allStar(adddist=True,distredux='v402')
    jk= alldata['J0']-alldata['K0']
    data= alldata[(jk >0.8)*(alldata['DISO_GAL'] > 0.)]
    #To allow for XY pixelization, we will hack these
    data= esutil.numpy_util.add_fields(data,[('RC_GALR', float),
                                             ('RC_GALPHI', float),
                                             ('RC_GALZ', float)])
    XYZ= bovy_coords.lbd_to_XYZ(data['GLON'],
                                data['GLAT'],
                                data['DISO_GAL'],
                                degree=True)
    R,phi,Z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                          XYZ[:,1],
                                          XYZ[:,2],
                                          Xsun=8.,Zsun=0.025)
    data['RC_GALR']= R
    data['RC_GALPHI']= phi
    data['RC_GALZ']= Z
    return data
