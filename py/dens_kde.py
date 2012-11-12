# KDE density estimation
import copy
import numpy
from scipy import maxentropy
class densKDE:
    """Class for KDE density estimation"""
    def __init__(self,data,kernel='biweight',w=None,
                 scale=True,fit=False,h=1.,
                 variable=False,variablenitt=2,variableexp=0.5):
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
           variable= (False) if True, do variable width
           variablenitt= number of iterations to settle on variable bandwidth
           variableexp= (0.5) exponent to use in variable lambda
        OUTPUT:
        HISTORY:
           2012-11-12 - Written - Bovy (IAS)
        """
        self._setup_kernel(kernel)
        self._ndata= data.shape[0]
        self._dim= data.shape[1]
        if not w is None:
            self._w= w/numpy.sum(w)
        else:
            self._w= numpy.ones(self._ndata)/float(self._ndata)
        if scale:
            self._data= self._scale_data(data,self._w)
        else:
            self._data= data
        #For variable
        self._lambda= numpy.ones((self._ndata,self._dim))
        if fit:
            raise NotImplementedError("fit=True not implemented yet")
        else:
            if isinstance(h,str) and 'scott' in h.lower():
                self._h= self._ndata**(-1./(self._dim+4.))
            elif isinstance(h,str) and 'silver' in h.lower():
                self._h= (self._ndata*(self._dim+2.)/4.)**(-1./(self._dim+4.))
            else:
                self._h= h
            if isinstance(h,str) and not 'gauss' in kernel.lower():
                self._h*= 2. #Larger for finite size kernel
        if variable:
            self._setup_variable(variablenitt,variableexp)

    def __call__(self,x,h=None,log=False,scale=True):
        """
        NAME:
           __call__
        PURPOSE:
            return the density
        INPUT:
           log= (False) if True, return the log
           h= (None) if set, use this bandwidth
           scale= (True), if False, don't rescale first
        OUTPUT:
           density (or log)
        HISTORY:
           2012-11-12 - Written - Bovy (IAS)
        """
        if h is None:
            thish= self._h
        else:
            thish= h
        x= self._prepare_x(x,scale)
        thiskernel= self._kernel(numpy.tile(x,(self._ndata,1))/thish/self._lambda,
                                 self._data/thish/self._lambda,
                                 log=log)
        if log:
            return -self._dim*numpy.log(thish)\
                +maxentropy.logsumexp(thiskernel+numpy.log(self._w)\
                                          -self._dim*numpy.log(self._lambda[:,0])) #latter assumes that lambda are spherical
        else:
            return 1./thish**self._dim\
                *numpy.sum(self._w*thiskernel/self._lambda[:,0]**self._dim)

    def _setup_variable(self,nitt,alpha):
        for ii in range(nitt):
            logdens= numpy.array([self(self._data[ii,:],log=True,scale=False) for ii in range(self._ndata)])
            logg= numpy.mean(logdens)
            self._lambda= numpy.tile(numpy.exp(-alpha*(logdens-logg)),(self._dim,1)).T
        return None

    def _prepare_x(self,x,scale):
        x= copy.copy(x) #To avoid changing input
        if isinstance(x,list):
            x= numpy.reshape(numpy.array(x),(1,len(x)))
        elif isinstance(x,(float,numpy.float32,numpy.float64)):
            x= numpy.reshape(numpy.array([x]),(1,1))
        else:
            x= numpy.reshape(x,(1,numpy.prod(x.shape))) #FOR NOW
        if scale and self._scaled:
            x-= numpy.tile(self._scalem,(x.shape[0],1))
            x/= numpy.tile(self._scales,(x.shape[0],1))
        return x

    def _scale_data(self,data,w):
        """Scale the data"""
        self._scaled= True
        m= numpy.sum(numpy.tile(w,(self._dim,1)).T*data,axis=0)/numpy.sum(w)
        s= numpy.sqrt(numpy.sum(numpy.tile(w,(self._dim,1)).T*data**2.,axis=0)/numpy.sum(w)-m**2.)
        self._scalem= m
        self._scales= s
        return (data-numpy.tile(m,(self._ndata,1)))\
            /numpy.tile(s,(self._ndata,1)) #simple scaling by std dev

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
