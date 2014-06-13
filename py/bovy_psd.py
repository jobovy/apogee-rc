import numpy
from radialprofile import azimuthalAverage
def psd(arr,dx,pad=True):
    """
    NAME:
       psd
    PURPOSE:
       Calculate the power spectrum of an array
    INPUT:
       arr - array
       dx - sampling interval
    OUTPUT:
       (frequencies,
       the power spectrum using definitions of NR 13.4 eqn. (13.4.5))
    HISTORY:
       2014-06-12 - Written - Bovy (IAS)
    """
    #First rm NaN
    arr[numpy.isnan(arr)]= 0.
    #Take the 1D FFT
    if pad:
        n= 2*len(arr)+1
    else:
        n= len(arr)
    arr_fft= numpy.fft.fftshift(numpy.fft.fft(arr,n=n))
    #Then calculate the periodogram estimate of the power spectrum
    ret= numpy.abs(arr_fft)**2.\
        +numpy.abs(arr_fft[::-1])**2.
    if pad:
        ret[len(arr)]*= 0.5 #correct zero order term
    else:
        ret[len(arr)/2]*= 0.5 #correct zero order term
    return (numpy.fft.fftshift(numpy.fft.fftfreq(n,dx)),
            ret/n**2.*2.)

def psd2d(image,pad=True):
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
    if pad:
        image_fft= numpy.fft.fftshift(numpy.fft.fft2(image,
                                                     s=(2*image.shape[0]+1,
                                                        2*image.shape[1]+1)))
    else:
        image_fft= numpy.fft.fftshift(numpy.fft.fft2(image,
                                                     s=(image.shape[0],
                                                        image.shape[1])))
    #Then calculate the periodogram estimate of the power spectrum
    ret= numpy.abs(image_fft)**2.\
        +numpy.abs(image_fft[::-1,::-1])**2.
    if pad:
        ret[image.shape[0],image.shape[1]]*= 0.5 #correct zero order term
        return ret/(2.*image.shape[0]+1)**2./(2.*image.shape[1]+1)**2.*4.
    else:
        ret[image.shape[0]/2,image.shape[1]/2]*= 0.5 #correct zero order term
        return ret/image.shape[0]**2./image.shape[1]**2.*16.

def psd1d(image,dx,binsize=1.,pad=True):
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
    ret2d= psd2d(image,pad=pad)
    nr, radii,ret= azimuthalAverage(ret2d,returnradii=False,binsize=binsize,
                                    dx=1./2./dx/image.shape[0],
                                    dy=1./2./dx/image.shape[1],
                                    interpnan=False,
                                    return_nr=True)
    return (radii/2.**pad/dx/(image.shape[0]/2.-0.),ret,ret/numpy.sqrt(nr))
