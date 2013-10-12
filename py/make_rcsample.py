import sys
import os, os.path
import numpy
import fitsio
import esutil
from astroquery.vizier import Vizier
from astropy import units as u
import astropy.coordinates as coord
from galpy.util import bovy_coords
import isodist
import apogee.tools.read as apread
import rcmodel
_ADDHAYDENDIST= False
_ERASESTR= "                                                                                "
def make_rcsample(savefilename):
    #Read the base-sample
    data= apread.allStar(adddist=_ADDHAYDENDIST)
    #Select red-clump stars
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)
    logg= data['LOGG']
    indx= (jk < 0.8)*(jk > 0.5)\
        *(z <= 0.06)\
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
    rcd= rcmodel.rcdist('../../rcdist-apogee/data/rcmodel_mode_jkz_ks_parsec.sav')
    jk= data['J0']-data['K0']
    z= isodist.FEH2Z(data['METALS'],zsolar=0.017)
    data['RC_DIST']= rcd(jk,z,appmag=data['K0'])
    data['RC_DM']= 5.*numpy.log10(data['RC_DIST'])+10.
    XYZ= bovy_coords.lbd_to_XYZ(data['GLON'],
                                data['GLAT'],
                                data['RC_DIST'],
                                degree=True)
    R,phi,Z= bovy_coords.XYZ_to_galcencyl(XYZ[:,0],
                                          XYZ[:,1],
                                          XYZ[:,2],
                                          Xsun=8.,Zsun=0.025)
    data['RC_GALR']= R
    data['RC_GALPHI']= phi
    data['RC_GALZ']= Z
    #Get proper motions
    pmfile= savefilename.split('.')[0]+'_pms.fits'
    if os.path.exists(pmfile):
        pmdata= fitsio.read(pmfile,1)
    else:
        pmdata= numpy.recarray(len(data),
                               formats=['f8','f8','f8','f8','f8','f8','i4'],
                               names=['RA','DEC','PMRA','PMDEC',
                                      'PMRA_ERR','PMDEC_ERR','PMMATCH'])
        rad= u.Quantity(4./3600.,u.degree)
        v= Vizier(columns=['RAJ2000','DEJ2000','pmRA','pmDE','e_pmRA','e_pmDE'])
        for ii in range(len(data)):
            #if ii > 100: break
            sys.stdout.write('\r'+"Getting pm data for point %i / %i" % (ii+1,len(data)))
            sys.stdout.flush()
            pmdata.RA[ii]= data['RA'][ii]
            pmdata.DEC[ii]= data['DEC'][ii]
            co= coord.ICRSCoordinates(ra=data['RA'][ii],
                                      dec=data['DEC'][ii],
                                      unit=(u.degree, u.degree))
            tab= v.query_region(co,rad,catalog='I/322') #UCAC-4 catalog
            if len(tab) == 0:
                pmdata.PMMATCH[ii]= 0
                print "Didn't find a match for %i ..." % ii
                continue
            else:
                pmdata.PMMATCH[ii]= len(tab)
                if len(tab[0]['pmRA']) > 1:
                    print "Found more than 1 match for %i ..." % ii
            try:
                pmdata.PMRA[ii]= float(tab[0]['pmRA'])
            except TypeError:
                jj= 1
                while len(tab[0]['pmRA']) > 1 and jj < 4: 
                    trad= u.Quantity((4.-jj)/3600.,u.degree)
                    tab= v.query_region(co,trad,catalog='I/322') #UCAC-4 catalog
                    jj+= 1
                if len(tab) == 0:
                    pmdata.PMMATCH[ii]= 0
                    print "Didn't find a unambiguous match for %i ..." % ii
                    continue               
                pmdata.PMRA[ii]= float(tab[0]['pmRA'])
            pmdata.PMDEC[ii]= float(tab[0]['pmDE'])
            pmdata.PMRA_ERR[ii]= float(tab[0]['e_pmRA'])
            pmdata.PMDEC_ERR[ii]= float(tab[0]['e_pmDE'])
            if numpy.isnan(float(tab[0]['pmRA'])): pmdata.PMMATCH[ii]= 0
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        fitsio.write(pmfile,pmdata,clobber=True)
        #To make sure we're using the same format below
        pmdata= fitsio.read(pmfile,1)
    #Match proper motions
    try: #These already exist currently, but may not always exist
        data= esutil.numpy_util.remove_fields(data,['PMRA','PMDEC'])
    except ValueError:
        pass
    data= esutil.numpy_util.add_fields(data,[('PMRA', numpy.float),
                                             ('PMDEC', numpy.float),
                                             ('PMRA_ERR', numpy.float),
                                             ('PMDEC_ERR', numpy.float),
                                             ('PMMATCH',int)])
    data['PMMATCH']= 0
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(pmdata['RA'],pmdata['DEC'],
                        data['RA'],data['DEC'],
                        2./3600.,maxmatch=1)
    data['PMRA'][m2]= pmdata['PMRA'][m1]
    data['PMDEC'][m2]= pmdata['PMDEC'][m1]
    data['PMRA_ERR'][m2]= pmdata['PMRA_ERR'][m1]
    data['PMDEC_ERR'][m2]= pmdata['PMDEC_ERR'][m1]
    data['PMMATCH'][m2]= pmdata['PMMATCH'][m1].astype(int)
    pmindx= data['PMMATCH'] == 1
    data['PMRA'][True-pmindx]= -9999.99
    data['PMDEC'][True-pmindx]= -9999.99
    data['PMRA_ERR'][True-pmindx]= -9999.99
    data['PMDEC_ERR'][True-pmindx]= -9999.99
    #Calculate Galactocentric velocities
    data= esutil.numpy_util.add_fields(data,[('GALVR', numpy.float),
                                             ('GALVT', numpy.float),
                                             ('GALVZ', numpy.float)])
    lb= bovy_coords.radec_to_lb(data['RA'],data['DEC'],degree=True)
    XYZ= bovy_coords.lbd_to_XYZ(lb[:,0],lb[:,1],data['RC_DIST'],degree=True)
    pmllpmbb= bovy_coords.pmrapmdec_to_pmllpmbb(data['PMRA'],data['PMDEC'],
                                                data['RA'],data['DEC'],
                                                degree=True)
    vxvyvz= bovy_coords.vrpmllpmbb_to_vxvyvz(data['VHELIO_AVG'],
                                             pmllpmbb[:,0],
                                             pmllpmbb[:,1],
                                             lb[:,0],lb[:,1],data['RC_DIST'],
                                             degree=True)
    vR, vT, vZ= bovy_coords.vxvyvz_to_galcencyl(vxvyvz[:,0],
                                                vxvyvz[:,1],
                                                vxvyvz[:,2],
                                                8.-XYZ[:,0],
                                                XYZ[:,1],
                                                XYZ[:,2]+0.025,
                                                vsun=[-11.1,30.24*8.,7.25])#Assumes proper motion of Sgr A* and R0=8 kpc, zo= 25 pc
    data['GALVR']= vR
    data['GALVT']= vT
    data['GALVZ']= vZ
    data['GALVR'][True-pmindx]= -9999.99
    data['GALVT'][True-pmindx]= -9999.99
    data['GALVZ'][True-pmindx]= -9999.99
    #Get proper motions
    pmfile= savefilename.split('.')[0]+'_pms_ppmxl.fits'
    if os.path.exists(pmfile):
        pmdata= fitsio.read(pmfile,1)
    else:
        pmdata= numpy.recarray(len(data),
                               formats=['f8','f8','f8','f8','f8','f8','i4'],
                               names=['RA','DEC','PMRA','PMDEC',
                                      'PMRA_ERR','PMDEC_ERR','PMMATCH'])
        rad= u.Quantity(4./3600.,u.degree)
        v= Vizier(columns=['RAJ2000','DEJ2000','pmRA','pmDE','e_pmRA','e_pmDE'])
        for ii in range(len(data)):
            #if ii > 100: break
            sys.stdout.write('\r'+"Getting pm data for point %i / %i" % (ii+1,len(data)))
            sys.stdout.flush()
            pmdata.RA[ii]= data['RA'][ii]
            pmdata.DEC[ii]= data['DEC'][ii]
            co= coord.ICRSCoordinates(ra=data['RA'][ii],
                                      dec=data['DEC'][ii],
                                      unit=(u.degree, u.degree))
            tab= v.query_region(co,rad,catalog='I/317') #PPMXL catalog
            if len(tab) == 0:
                pmdata.PMMATCH[ii]= 0
                print "Didn't find a match for %i ..." % ii
                continue
            else:
                pmdata.PMMATCH[ii]= len(tab)
                if len(tab[0]['pmRA']) > 1:
                    pass
                    #print "Found more than 1 match for %i ..." % ii
            try:
                pmdata.PMRA[ii]= float(tab[0]['pmRA'])
            except TypeError:
                #Find nearest
                cosdists= numpy.zeros(len(tab[0]['pmRA']))
                for jj in range(len(tab[0]['pmRA'])):
                    cosdists[jj]= cos_sphere_dist(tab[0]['RAJ2000'][jj],
                                                  tab[0]['DEJ2000'][jj],
                                                  data['RA'][ii],
                                                  data['DEC'][ii])
                closest= numpy.argmax(cosdists)
                pmdata.PMRA[ii]= float(tab[0]['pmRA'][closest])
                pmdata.PMDEC[ii]= float(tab[0]['pmDE'][closest])
                pmdata.PMRA_ERR[ii]= float(tab[0]['e_pmRA'][closest])
                pmdata.PMDEC_ERR[ii]= float(tab[0]['e_pmDE'][closest])
                if numpy.isnan(float(tab[0]['pmRA'][closest])): pmdata.PMMATCH[ii]= 0
            else:
                pmdata.PMDEC[ii]= float(tab[0]['pmDE'])
                pmdata.PMRA_ERR[ii]= float(tab[0]['e_pmRA'])
                pmdata.PMDEC_ERR[ii]= float(tab[0]['e_pmDE'])
                if numpy.isnan(float(tab[0]['pmRA'])): pmdata.PMMATCH[ii]= 0
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
        fitsio.write(pmfile,pmdata,clobber=True)
        #To make sure we're using the same format below
        pmdata= fitsio.read(pmfile,1)
    #Match proper motions to ppmxl
    data= esutil.numpy_util.add_fields(data,[('PMRA_PPMXL', numpy.float),
                                             ('PMDEC_PPMXL', numpy.float),
                                             ('PMRA_ERR_PPMXL', numpy.float),
                                             ('PMDEC_ERR_PPMXL', numpy.float),
                                             ('PMMATCH_PPMXL',int)])
    data['PMMATCH_PPMXL']= 0
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(pmdata['RA'],pmdata['DEC'],
                        data['RA'],data['DEC'],
                        2./3600.,maxmatch=1)
    data['PMRA_PPMXL'][m2]= pmdata['PMRA'][m1]
    data['PMDEC_PPMXL'][m2]= pmdata['PMDEC'][m1]
    data['PMRA_ERR_PPMXL'][m2]= pmdata['PMRA_ERR'][m1]
    data['PMDEC_ERR_PPMXL'][m2]= pmdata['PMDEC_ERR'][m1]
    data['PMMATCH_PPMXL'][m2]= pmdata['PMMATCH'][m1].astype(int)
    pmindx= data['PMMATCH_PPMXL'] == 1
    data['PMRA_PPMXL'][True-pmindx]= -9999.99
    data['PMDEC_PPMXL'][True-pmindx]= -9999.99
    data['PMRA_ERR_PPMXL'][True-pmindx]= -9999.99
    data['PMDEC_ERR_PPMXL'][True-pmindx]= -9999.99
    #Calculate Galactocentric velocities
    data= esutil.numpy_util.add_fields(data,[('GALVR_PPMXL', numpy.float),
                                             ('GALVT_PPMXL', numpy.float),
                                             ('GALVZ_PPMXL', numpy.float)])
    lb= bovy_coords.radec_to_lb(data['RA'],data['DEC'],degree=True)
    XYZ= bovy_coords.lbd_to_XYZ(lb[:,0],lb[:,1],data['RC_DIST'],degree=True)
    pmllpmbb= bovy_coords.pmrapmdec_to_pmllpmbb(data['PMRA_PPMXL'],
                                                data['PMDEC_PPMXL'],
                                                data['RA'],data['DEC'],
                                                degree=True)
    vxvyvz= bovy_coords.vrpmllpmbb_to_vxvyvz(data['VHELIO_AVG'],
                                             pmllpmbb[:,0],
                                             pmllpmbb[:,1],
                                             lb[:,0],lb[:,1],data['RC_DIST'],
                                             degree=True)
    vR, vT, vZ= bovy_coords.vxvyvz_to_galcencyl(vxvyvz[:,0],
                                                vxvyvz[:,1],
                                                vxvyvz[:,2],
                                                8.-XYZ[:,0],
                                                XYZ[:,1],
                                                XYZ[:,2]+0.025,
                                                vsun=[-11.1,30.24*8.,7.25])#Assumes proper motion of Sgr A* and R0=8 kpc, zo= 25 pc
    data['GALVR_PPMXL']= vR
    data['GALVT_PPMXL']= vT
    data['GALVZ_PPMXL']= vZ
    data['GALVR_PPMXL'][True-pmindx]= -9999.99
    data['GALVT_PPMXL'][True-pmindx]= -9999.99
    data['GALVZ_PPMXL'][True-pmindx]= -9999.99
    #Save
    fitsio.write(savefilename,data,clobber=True)
    return None

def cos_sphere_dist(theta,phi,theta_o,phi_o):
    """
    NAME:
       cos_sphere_dist
    PURPOSE:
       computes the cosine of the spherical distance between two
       points on the sphere
    INPUT:
       theta  - polar angle [0,pi]
       phi    - azimuth [0,2pi]
       theta  - polar angle of center of the disk
       phi_0  - azimuth of the center of the disk
    OUTPUT:
       spherical distance
    HISTORY:
       2010-04-29 -Written - Bovy (NYU)
    """
    return (numpy.sin(theta)*numpy.sin(theta_o)*(numpy.cos(phi_o)*numpy.cos(phi)+
                                                 numpy.sin(phi_o)*numpy.sin(phi))+
            numpy.cos(theta_o)*numpy.cos(theta))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        make_rcsample(sys.argv[1])
    else:
        #Create savefilename
        import apogee.tools.path as appath
        savefilename= os.path.join(appath._APOGEE_DATA,
                                   'rcsample_'+appath._APOGEE_REDUX+'.fits')
        print "Saving to %s ..." % savefilename
        make_rcsample(savefilename)
