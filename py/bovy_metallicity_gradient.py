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
# Plotting options
_PLOTMAD= False
_NORMFE= True
_HIZ= True
# Lines
line_labels= {}
line_labels['fe']= r'$\mathrm{Fe\ I}$'
line_labels['mg']= r'$\mathrm{Mg\ I}$'
line_labels['al']= r'$\mathrm{Al\ I}$'
line_labels['si']= r'$\mathrm{Si\ I}$'
line_labels['k']= r'$\mathrm{K\ I}$'
line_labels['ca']= r'$\mathrm{Ca\ I}$'
line_labels['ti']= r'$\mathrm{Ti\ I}$'
line_labels['cr']= r'$\mathrm{Cr\ I}$'
line_labels['ni']= r'$\mathrm{Ni\ I}$'
_FEI_lines= [15198.644,15211.682,15399.925,15494.572,15652.786,15969.229,
             16045.040,16157.660,16169.448]
_MGI_lines= [15745.017,15753.203,15770.108,15883.839,15890.541,15893.826,
             15958.836]
_ALI_lines= [16723.524,16767.938]
_SII_lines= [15365.359,15381.033,15837.928,15964.424,16064.397,16099.184,
             16220.100,16685.327,16832.756]
_KI_lines= [15167.211,15172.521]
_CAI_lines= [16141.232,16155.176,16159.650,16161.778]
_TII_lines= [15548.003,15607.106,15703.269,15719.867,16639.705]
_CRI_lines= [15684.348,15864.548]
_NII_lines= [15609.944,15636.926,16588.970,16593.827,16678.266,16820.064,
             16823.354]
def bovy_metallicity_gradient(plotfilename,savefilename,largewave=False):
    # First read the RC catalog and cut it to stars near the plane
    data= apread.rcsample()
    if _HIZ:
        indx= (numpy.fabs(data['RC_GALZ']) > 0.6)*(data['METALS'] > -1000.)
    else:
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
                if _PLOTMAD:
                    median_spec[jj,ii]= \
                        numpy.median(numpy.fabs(allspec[jj,True-numpy.isnan(allspec[jj,:])]-numpy.median(allspec[jj,True-numpy.isnan(allspec[jj,:])])))
                else:
                    median_spec[jj,ii]= \
                        numpy.median(allspec[jj,
                                             True-numpy.isnan(allspec[jj,:])])
        save_pickles(savefilename,median_spec)
    # Wavelengths
    wave= 10.**(numpy.arange(hdr['CRVAL1'],
                             hdr['CRVAL1']+len(spec)*hdr['CDELT1'],
                             hdr['CDELT1']))
    # Normalization, first calculate the spectrum near Ro
    if _NORMFE:
        absmax= 0.965
    else:
        absmax= 0.98
    absindx= median_spec > absmax
    median_spec[absindx]= absmax #focus on real absorption lines
    rospec= numpy.zeros(len(spec))
    roindx= numpy.argmin(numpy.fabs(Rs-8.))
    for jj in range(len(spec)):
        rospec[jj]= numpy.median(median_spec[jj,roindx-3:roindx+4])
    if _NORMFE and not _PLOTMAD:
        # Normalize by the difference in Fe
        feindx= numpy.zeros(len(spec),dtype='bool')
        for line in _FEI_lines:
            wavindx= numpy.argmin(numpy.fabs(wave-line))
            feindx[wavindx-2:wavindx+3]= True
            feindx[wavindx-4:wavindx+5]= True
        median_spec= numpy.log(median_spec)
        rospec= numpy.log(rospec)
        for jj in range(nR):
            fehdiff= numpy.median(median_spec[feindx,jj]-rospec[feindx])
#            print fehdiff, numpy.median(numpy.fabs(median_spec[feindx,jj]-rospec[feindx]-fehdiff))
            median_spec[:,jj]= (median_spec[:,jj]-rospec)-fehdiff
        median_spec[absindx]= 0.
        vmin, vmax= -0.02, 0.02
    elif not _PLOTMAD:
        # Normalization by spectrum at Ro
        roindx= numpy.argmin(numpy.fabs(Rs-8.))
        for jj in range(nR):
            # Normalize by the solar radius
            median_spec[:,jj]/= rospec
        median_spec-= 1. 
        vmin=-0.035
        vmax=0.035
    # Now plot
    if False:
        startindx, endindx= 3652, 4100#3915
    if largewave:
        startindx, endindx= 7000,7600#7375, 7889
    else:
        startindx, endindx= 2500, 3100
    bovy_plot.bovy_print(fig_width=7.,fig_height=4.)
    bovy_plot.bovy_dens2d(-median_spec[startindx:endindx,:].T,
                           origin='lower',cmap='coolwarm',#cmap='afmhot',
                           vmin=vmin,vmax=vmax,
#                          colorbar=True,shrink=0.78,
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
    # Label the lines
    _label_lines('fe',wave[startindx],wave[endindx])
    _label_lines('mg',wave[startindx],wave[endindx])
    _label_lines('si',wave[startindx],wave[endindx])
    _label_lines('al',wave[startindx],wave[endindx])
    _label_lines('k',wave[startindx],wave[endindx])
    _label_lines('cr',wave[startindx],wave[endindx])
    _label_lines('ca',wave[startindx],wave[endindx])
    _label_lines('ti',wave[startindx],wave[endindx])
    _label_lines('ni',wave[startindx],wave[endindx])
    bovy_plot.bovy_end_print(plotfilename)
    return None

def _label_lines(elem,wavemin,wavemax):
    if elem.lower() == 'fe':
        lines= _FEI_lines
    elif elem.lower() == 'mg':
        lines= _MGI_lines
    elif elem.lower() == 'al':
        lines= _ALI_lines
    elif elem.lower() == 'si':
        lines= _SII_lines
    elif elem.lower() == 'k':
        lines= _KI_lines
    elif elem.lower() == 'cr':
        lines= _CRI_lines
    elif elem.lower() == 'ca':
        lines= _CAI_lines
    elif elem.lower() == 'ti':
        lines= _TII_lines
    elif elem.lower() == 'ni':
        lines= _NII_lines
    for line in lines:
        if line > wavemin and line < wavemax:
            bovy_plot.bovy_plot([numpy.log10(line)-numpy.log10(15000.),
                                 numpy.log10(line)-numpy.log10(15000.)],
                                [12.5,13.],'w-',overplot=True)
            bovy_plot.bovy_text(numpy.log10(line)-numpy.log10(15000.),
                                13.1,line_labels[elem.lower()],
                                size=10.)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 3:
        bovy_metallicity_gradient(sys.argv[1],sys.argv[2],largewave=True)
    else:
        bovy_metallicity_gradient(sys.argv[1],sys.argv[2])
