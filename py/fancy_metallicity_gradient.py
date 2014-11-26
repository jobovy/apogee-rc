# fancy_metallicity_gradient: make a fancy plot a la Fan or Zhu & Menard of spectra as a function of R to illustrate the metallicity gradient
import sys
import numpy
from galpy.util import bovy_plot, save_pickles
import apogee.tools.read as apread
def fancy_metallicity_gradient(plotfilename,savefilename):
    # First read the RC catalog and cut it to stars near the plane
    data= apread.rcsample()
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    data= data[indx]
    # Now go through bins in R
    Rmin, Rmax, dR= 5.5, 13., 0.1
    Rs= numpy.arange(Rmin+dR/2.,Rmax+dR/2.,dR)
    nR= len(Rs)
    # Read one spectrum to learn the size of the array
    spec, hdr= apread.aspcapStar(4424,'2M00025587+5849278',ext=1)
    if os.path.exists(savefilename):
        # Reload previously calculated median spectrum
        savefile= open(savefilename,'rb')
        median_spec= pickle.load(savefile)
        savefile.close()
    else:
        median_spec= numpy.empty((len(spec),nR))
        # Now run through all the spectra and get the median
        for ii in range(nR):
            indx= (data['RC_GALR'] >= Rs[ii]-dR/2.)\
                *(data['RC_GALR'] < Rs[ii]+dR/2.)
            allspec= numpy.empty((len(spec),numpy.sum(indx)))
            for jj in range(numpy.sum(indx)):
                specdata= \
                    apread.aspcapStar(data['LOCATION_ID'][indx][jj],
                                      data['APOGEE_ID'][indx][jj],
                                      ext=1,header=False)
                specerr= \
                    apread.aspcapStar(data['LOCATION_ID'][indx][jj],
                                      data['APOGEE_ID'][indx][jj],
                                      ext=2,header=False)
                allspec[:,jj]= specdata
                allspec[specerr > 1.,jj]= numpy.nan
            median_spec[:,ii]= numpy.median(allspec[True-numpy.isnan(allspec)],
                                            axis=1)
        save_pickles(savefilename,median_spec)
    # Now plot
    wave= 10.**(numpy.arange(hdr['CRVAL1'],
                             hdr['CRVAL1']+len(spec)*hdr['CDELT1'],
                             hdr['CDELT1']))
    bovy_plot.bovy_print()
    bovy_plot.bovy_dens2d(median_spec.T,origin='lower',cmap='afmhot',
                          vmin=0.5,vmax=1.,
                          xrange=[wave[0],wave[-1]],
                          yrange=[Rs[0]-dR/2.,Rs[-1]+dR/2.],
                          xlabel=r'$\lambda\,(\AA)$',
                          ylabel=r'$R\,(\mathrm{kpc})$')
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    fancy_metallicity_gradient(sys.argv[1],sys.argv[2]))
