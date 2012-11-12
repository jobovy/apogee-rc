# KDE density estimation
import numpy
from scipy import maxentropy
class densKDE:
    """Class for KDE density estimation"""
    def __init__(self,data,kernel='biweight',w=None,
                 scale=True,fit=False,h=1.):
        """
        NAME:
           __init__
        PURPOSE:
           initialize
        INPUT:
           data - [N,dim] array of data of dimensionality dim
           kernel - function kernel(x,y,log=False) 
                    that returns the (properly normalized)
                    kernel for input x[ndata,dim], or one of ['biweight']
           scale= (True), if False, don't scale the data
           fit= (False) if True, determine best bandwith h
           h= (1.) bandwidth
        OUTPUT:
        HISTORY:
           2012-11-12 - Written - Bovy (IAS)
        """
        self._setup_kernel(kernel)
        if scale:
            self._data= self._scale_data(data)
        else:
            self._data= data
        self._ndata= self._data.shape[0]
        self._dim= self._data.shape[1]
        if not w is None:
            self._w= w/numpy.sum(w)
        else:
            self._w= numpy.ones(self._ndata)/float(self._ndata)
        if fit:
            raise NotImplementedError("fit=True not implemented yet")
        else:
            self._h= h

    def __call__(self,x,h=None,log=False):
        """
        NAME:
           __call__
        PURPOSE:
            return the density
        INPUT:
           log= (False) if True, return the log
           h= (None) if set, use this bandwidth
        OUTPUT:
           density (or log)
        HISTORY:
           2012-11-12 - Written - Bovy (IAS)
        """
        if h is None:
            thish= self._h
        else:
            thish= h
        x= self._prepare_x(x)
        thiskernel= self._kernel(numpy.tile(x,(self._ndata,1))/thish,
                                 self._data/thish,
                                 log=log)
        if log:
            return -self._dim*numpy.log(thish)\
                +maxentropy.logsumexp(thiskernel+numpy.log(self._w))
        else:
            return 1./thish**self._dim\
                *numpy.sum(self._w*thiskernel)

    def _prepare_x(self,x):
        if isinstance(x,list):
            x= numpy.reshape(numpy.array(x),(1,len(x)))
        elif isinstance(x,(float,numpy.float32,numpy.float64)):
            x= numpy.reshape(numpy.array([x]),(1,1))
        if self._scaled:
            x/= numpy.tile(self._scales,(x.shape[0],1))
        return x

    def _scale_data(self,data):
        """Scale the data"""
        self._scaled= True
        s= numpy.std(data,axis=0)
        self._scales= s
        return data/numpy.tile(s,(data.shape[0],1)) #simple scaling by std dev

    def _setup_kernel(self,kernel):
        """Parse the kernel input"""
        if isinstance(kernel,str):
            if kernel.lower() == 'biweight':
                self._kernel= kernel_biweight
            elif kernel.lower() == 'gauss' or kernel.lower() == 'gaussian':
                self._kernel= kernel_gauss
        else:
            self._kernel= kernel

def kernel_biweight(x,y,log=False):
    x, y= preparexy(x,y)
    r2= numpy.sum((x-y)**2.,axis=1)
    indx= (r2 < 1.)
    out= numpy.empty(x.shape[0])
    if log:
        out[indx]= numpy.log(3./numpy.pi*(1.-r2[indx])**2.)
    else:
        out[indx]= 3./numpy.pi*(1.-r2[indx])**2.
    if log:
        out[True-indx]= -numpy.finfo(numpy.dtype(numpy.float64)).max
    else:
        out[True-indx]= 0.
    return out

def kernel_gauss(x,y,log=False):
    x, y= preparexy(x,y)
    dim= x.shape[1]
    r2= numpy.sum((x-y)**2.,axis=1)
    if log:
        return -dim/2.*numpy.log(2.*numpy.pi)-r2/2.
    else:
        return 1./(2.*numpy.pi)**dim/2.*numpy.exp(-r2/2.)

def preparexy(x,y):
    return (preparex(x),preparex(y))
def preparex(x):
    if isinstance(x,list):
        x= numpy.array(x)
    elif isinstance(x,(float,numpy.float32,numpy.float64)):
        x= numpy.array([x])
    return x
