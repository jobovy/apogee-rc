# bovy_metallicity_gradient: make a fancy plot a la Fan or Zhu & Menard of spectra as a function of R to illustrate the metallicity gradient
#
# Use as:
#
# python bovy_metallicity_gradient.py bovy_metal.png bovy_metal.sav
#
# Requires:
#
#   apogee package and associated data (https://github.com/jobovy/apogee; spectra will be downloaded automatically, but slow)
#   galpy: for plotting (https://github.com/jobovy/galpy)
#
import os
import sys
import pickle
import numpy
from galpy.util import bovy_plot, save_pickles
import apogee.tools.read as apread
def bovy_metallicity_gradient(plotfilename,savefilename):
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
        median_spec= numpy.zeros((len(spec),nR))
        # Now run through all the spectra and get the median
        tot= 0
        for ii in range(nR):
            indx= (data['RC_GALR'] >= Rs[ii]-dR/2.)\
                *(data['RC_GALR'] < Rs[ii]+dR/2.)
            tot+= numpy.sum(indx)
            print numpy.sum(indx), tot
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
            for jj in range(len(spec)):
                median_spec[jj,ii]= \
                    numpy.median(allspec[jj,True-numpy.isnan(allspec[jj,:])])
        save_pickles(savefilename,median_spec)
    # Normalization
    median_spec[median_spec > .98]= .98 #focus on real absorption lines
    roindx= numpy.argmin(numpy.fabs(Rs-8.))
    for jj in range(len(spec)):
        #Normalize by the solar radius
        median_spec[jj,:]/= numpy.median(median_spec[jj,roindx-3:roindx+4])
    # Now plot
    startindx, endindx= 3652, 4100#3915
    bovy_plot.bovy_print(fig_width=7.,fig_height=4.)
    bovy_plot.bovy_dens2d((1.-(median_spec[startindx:endindx,:]-1.)).T,
                          origin='lower',cmap='coolwarm',#cmap='afmhot',
                          vmin=0.965,vmax=1.035,
                          interpolation='bicubic',
                          aspect=(hdr['CRVAL1']+(endindx-0.5)*hdr['CDELT1']-hdr['CRVAL1']-(startindx-0.5)*hdr['CDELT1'])/(Rs[-1]+dR/2.-Rs[0]+dR/2.)/2.,
                          xrange=[hdr['CRVAL1']+(startindx-0.5)*hdr['CDELT1']-numpy.log10(15000.),
                                  hdr['CRVAL1']+(endindx-0.5)*hdr['CDELT1']-numpy.log10(15000)],
                          yrange=[Rs[0]-dR/2.,Rs[-1]+dR/2.],
                          xlabel=r'$\log \lambda / 15,000 \AA$',
                          ylabel=r'$R\,(\mathrm{kpc})$')
    # Draw solar line
    bovy_plot.bovy_plot([-1000000000,100000000],[8.,8.],'k--',lw=1.5,
                        overplot=True)
    bovy_plot.bovy_end_print(plotfilename)
    return None

if __name__ == '__main__':
    bovy_metallicity_gradient(sys.argv[1],sys.argv[2])
