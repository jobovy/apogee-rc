import numpy
from radialprofile import azimuthalAverage
def psd2d(image):
    """
    NAME:
       psd2d
    PURPOSE:
       Calculate the 2D power spectrum of an image
    INPUT:
       image - the image [NxM]
    OUTPUT:
       the 2D power spectrum using definitions of NR 13.4 eqn. (13.4.5)
    HISTORY:
       2014-06-06 - Written - Bovy (IAS)
    """
    #First rm NaN
    image[numpy.isnan(image)]= 0.
    #First take the 2D FFT
    image_fft= numpy.fft.fft2(image,s=(2*image.shape[0],2*image.shape[1]))
    #Then calculate the periodogram estimate of the power spectrum
    ret= numpy.abs(image_fft)**2.\
        +numpy.abs(numpy.roll(image_fft[::-1,::-1],1))**2.
    ret[0,0]*= 0.5 #correct zero order term
    if numpy.all([numpy.mod(image.shape[ii],2) for ii in range(2)] == 0):
        ret[image.shape[0]/2,image.shape[1]/2]*= 0.5
    return numpy.fft.fftshift(ret)/numpy.prod(image.shape)*2.

def psd1d(image,dx,binsize=1.):
    """
    NAME:
       psd1d
    PURPOSE:
       Calculate the 1D, azimuthally averaged power spectrum of an image
    INPUT:
       image - the image [NxM]
       dx- spacing in X and Y directions
       binsize= (1) radial binsize in terms of Nyquist frequency
    OUTPUT:
       the 1D power spectrum using definitions of NR 13.4 eqn. (13.4.5)
    HISTORY:
       2014-06-06 - Written - Bovy (IAS)
    """
    #First subtract DC component
    image-= numpy.mean(image[True-numpy.isnan(image)])
    #Calculate the 2D PSD
    ret2d= psd2d(image)
    nr, radii,ret= azimuthalAverage(ret2d,returnradii=False,binsize=binsize,
                                    dx=1./2./dx/image.shape[0],
                                    dy=1./2./dx/image.shape[1],
                                    interpnan=False,
                                    return_nr=True)
    return (radii/2./dx/(image.shape[0]/2.-1.),ret,ret/numpy.sqrt(nr))
