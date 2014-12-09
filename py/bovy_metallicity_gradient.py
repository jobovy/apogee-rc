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
from matplotlib import pyplot
import matplotlib.ticker as ticker
from matplotlib.ticker import NullFormatter
import apogee.tools.read as apread
# Plotting options
_PLOTMAD= False
_HIZ= False
# Lines
line_labels= {}
line_labels['fe']= r'$\mathrm{Fe\kern 0.1em I}$'
line_labels['mg']= r'$\mathrm{Mg\kern 0.1em I}$'
line_labels['al']= r'$\mathrm{Al\kern 0.1em I}$'
line_labels['si']= r'$\mathrm{Si\kern 0.1em I}$'
line_labels['k']= r'$\mathrm{K\kern 0.1em I}$'
line_labels['ca']= r'$\mathrm{Ca\kern 0.1em I}$'
line_labels['ti']= r'$\mathrm{Ti\kern 0.1em I}$'
line_labels['cr']= r'$\mathrm{Cr\kern 0.1em I}$'
line_labels['ni']= r'$\mathrm{Ni\kern 0.1em I}$'
line_labels['na']= r'$\mathrm{Na\kern 0.1em I}$'
line_labels['mn']= r'$\mathrm{Mn\kern 0.1em I}$'
line_labels['s']= r'$\mathrm{S\kern 0.1em I}$'
line_labels['oh']= r'$\mathrm{OH}$'
line_labels['dib']= r'$\mathrm{DIB}$'
_FEI_lines= [15198.644,15211.682,15399.925,15494.572,15652.786,15969.229,
             16045.040,16157.660,16169.448,16697.635]
_MGI_lines= [15745.017,15753.203,15770.108,15883.839,15890.541,15893.826,
             15958.836]
_ALI_lines= [16723.524,16767.938]
_SII_lines= [15365.359,15381.033,15837.928,15964.424,16064.397,16099.184,
             16220.100,16685.327,16832.756]
_KI_lines= [15167.211,15172.521]
_CAI_lines= [16141.232,16155.176,16159.650,16161.778]
_TII_lines= [15548.003,15607.106,15703.269,15719.867,16639.705]
_CRI_lines= [15684.348,15864.548,15470.129]
_NII_lines= [15609.944,15636.926,16588.970,16593.827,16678.266,16820.064,
             16823.354]
_NAI_lines= [16378.346633274852,16393.340725803333]
_MNI_lines= [15677.437,16712.565]
_SI_lines= [15406.540,15426.490,15474.043,15482.712]
_OH_lines= [15505.5]
_DIB_lines= [15272.42] #from Zasowski et al. (2014)
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
                allspec[specerr == 0.,jj]= numpy.nan
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
    absmax= 0.98
    absindx= median_spec > absmax
    median_spec[absindx]= absmax #focus on real absorption lines
    rospec= numpy.zeros(len(spec))
    roindx= numpy.argmin(numpy.fabs(Rs-8.))
    for jj in range(len(spec)):
        rospec[jj]= numpy.median(median_spec[jj,roindx-3:roindx+4])
    if not _PLOTMAD:
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
    # Make N separate plots showing different wavelength regions
    bovy_plot.bovy_print(fig_width=8.,fig_height=3.,
                         axes_labelsize=10,text_fontsize=9,legend_fontsize=9,
                         xtick_labelsize=8,ytick_labelsize=8)
    startindxs= [322,1784,2707,3665,4880,5870,7178] #DIB is at 818
    endindxs= [355,1930,2857,3718,4925,5955,7400]
    nregions= len(startindxs)
    # Calculate the width of the plot
    dx= numpy.array([endindxs[ii]-startindxs[ii] for ii in range(nregions)],
                    dtype='float')
    specdx= numpy.sum(dx) # for later
    dx/= numpy.sum(dx)
    totdx= 0.85
    skipdx= 0.015
    dx*= (totdx-(nregions-1)*skipdx)
    for ii in range(nregions):
        # Setup the axes
        if ii == 0:
            left, bottom, width, height= 0.1, 0.1, dx[ii],0.8
        else:
            left, bottom, width, height= 0.1+numpy.cumsum(dx)[ii-1]+skipdx*ii,\
                0.1, dx[ii], 0.8
        thisax= pyplot.axes([left,bottom,width,height])
        fig= pyplot.gcf()
        fig.sca(thisax)
        startindx, endindx= startindxs[ii], endindxs[ii]
        aspect= (hdr['CRVAL1']+(endindx-0.5)*hdr['CDELT1']-hdr['CRVAL1']-(startindx-0.5)*hdr['CDELT1'])\
            /(Rs[-1]+dR/2.-Rs[0]+dR/2.)
        aspect/= dx[ii]*4.5
        yrange= [Rs[0]-dR/2.,Rs[-1]+dR/2.]
        xrange=[hdr['CRVAL1']+(startindx-0.5)*hdr['CDELT1']-numpy.log10(15000.),\
                    hdr['CRVAL1']+(endindx-0.5)*hdr['CDELT1']-numpy.log10(15000)]
        extent= [xrange[0],xrange[1],yrange[0],yrange[1]]
        pyplot.imshow(-median_spec[startindx:endindx,:].T,
                       origin='lower',cmap='coolwarm',
                       vmin=vmin,vmax=vmax,
                       extent=extent,
                       interpolation='nearest',
                       #interpolation='bicubic',
                       aspect=aspect)
        pyplot.axis(extent)
        pyplot.xlim(xrange[0],xrange[1])
        pyplot.ylim(yrange[0],yrange[1])
        thisax.xaxis.set_major_locator(ticker.MultipleLocator(0.0005))
        bovy_plot._add_ticks()
        if ii > 0:
            nullfmt   = NullFormatter()         # no labels
            thisax.yaxis.set_major_formatter(nullfmt)
        else:
            pyplot.ylabel(r'$R\,(\mathrm{kpc})$')
        # Remove spines between different wavelength regions
        if ii == 0:
            thisax.spines['right'].set_visible(False)
            thisax.tick_params(right=False,which='both')
        elif ii == (nregions-1):
            thisax.spines['left'].set_visible(False)
            thisax.tick_params(labelleft='off')
            thisax.tick_params(left=False,which='both')
        else:
            thisax.spines['left'].set_visible(False)
            thisax.spines['right'].set_visible(False)
            thisax.tick_params(labelleft='off')
            thisax.tick_params(left=False,which='both')
            thisax.tick_params(right=False,which='both')
        # Plot cut-out markers
        d = .015 # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=thisax.transAxes, color='k', clip_on=False)
        slope= 1./(dx[ii]+0.2*skipdx)/3.
        if ii == 0:
            thisax.plot((1-slope*d,1+slope*d),(-d,+d), **kwargs)
            thisax.plot((1-slope*d,1+slope*d),(1-d,1+d), **kwargs)
        elif ii == (nregions-1):
            thisax.plot((-slope*d,+slope*d),(-d,+d), **kwargs)
            thisax.plot((-slope*d,+slope*d),(1-d,1+d), **kwargs)
        else:
            thisax.plot((1-slope*d,1+slope*d),(-d,+d), **kwargs)
            thisax.plot((1-slope*d,1+slope*d),(1-d,1+d), **kwargs)
            thisax.plot((-slope*d,+slope*d),(-d,+d), **kwargs)
            thisax.plot((-slope*d,+slope*d),(1-d,1+d), **kwargs)
        # Draw solar line
        if ii == (nregions-1):
            xend= hdr['CRVAL1']+(endindx-0.5)*hdr['CDELT1']-numpy.log10(15000)
            # Total wavelength span, incl. skipdx parts
            totaldx= hdr['CDELT1']*specdx*\
                (1.+(nregions-1)*skipdx/(totdx-(nregions-1)*skipdx))
            thisax.plot((xend-totaldx,xend),(8.,8.),
                        color='k',ls='--',lw=1.5,marker='None',clip_on=False)
        # Label the lines
        _label_all_lines(wave[startindx],wave[endindx])
    # Add the x-axis label
    thisax= pyplot.axes([0.1,0.15,0.85,0.7])
    pyplot.gcf().sca(thisax)
    thisax.spines['left'].set_visible(False)
    thisax.spines['right'].set_visible(False)
    thisax.spines['bottom'].set_visible(False)
    thisax.spines['top'].set_visible(False)
    thisax.tick_params(labelleft='off')
    thisax.tick_params(left=False,which='both')
    thisax.tick_params(right=False,which='both')
    thisax.tick_params(labelbottom='off')
    thisax.tick_params(bottom=False,which='both')
    thisax.tick_params(top=False,which='both')
    pyplot.xlabel(r'$\log \lambda / 15,000 \AA$')
    thisax.set_zorder(-1)
    bovy_plot.bovy_end_print(plotfilename,dpi=300)
    return None

def _label_all_lines(wavemin,wavemax):
    _label_lines('fe',wavemin,wavemax)
    _label_lines('mg',wavemin,wavemax)
    _label_lines('si',wavemin,wavemax)
    _label_lines('al',wavemin,wavemax)
    _label_lines('k',wavemin,wavemax)
    _label_lines('cr',wavemin,wavemax)
    _label_lines('ca',wavemin,wavemax)
    _label_lines('ti',wavemin,wavemax)
    _label_lines('ni',wavemin,wavemax)
    _label_lines('na',wavemin,wavemax)
    _label_lines('mn',wavemin,wavemax)
    _label_lines('s',wavemin,wavemax)
    _label_lines('oh',wavemin,wavemax)
    _label_lines('dib',wavemin,wavemax)
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
    elif elem.lower() == 'na':
        lines= _NAI_lines
    elif elem.lower() == 'mn':
        lines= _MNI_lines
    elif elem.lower() == 's':
        lines= _SI_lines
    elif elem.lower() == 'oh':
        lines= _OH_lines
    elif elem.lower() == 'dib':
        lines= _DIB_lines
    fontsize= 5.5
    for line in lines:
        if line > wavemin and line < wavemax:
            bovy_plot.bovy_plot([numpy.log10(line)-numpy.log10(15000.),
                                 numpy.log10(line)-numpy.log10(15000.)],
                                [12.5,13.],'w-',overplot=True)
            if elem == 'ca' and line > 16154. and line < 16160.:
                bovy_plot.bovy_text(numpy.log10(line)-numpy.log10(15000.),
                                    13.4,line_labels[elem.lower()],
                                    size=fontsize)
            elif elem == 'mg' and line > 15887. and line < 15891.826:
                bovy_plot.bovy_text(numpy.log10(line)-numpy.log10(15000.),
                                    13.4,line_labels[elem.lower()],
                                    size=fontsize)
            else:
                bovy_plot.bovy_text(numpy.log10(line)-numpy.log10(15000.),
                                    13.1,line_labels[elem.lower()],
                                    size=fontsize)
    return None

if __name__ == '__main__':
    if len(sys.argv) > 3:
        bovy_metallicity_gradient(sys.argv[1],sys.argv[2],largewave=True)
    else:
        bovy_metallicity_gradient(sys.argv[1],sys.argv[2])
