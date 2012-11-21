###############################################################################
# rcmodel.py: isochrone model for the distribution in (J-Ks,M_H)
#             along the RC
#
# MAY NEED:
# export PYTHONPATH=~/Repos/extreme-deconvolution-ngerrors/py:$PYTHONPATH
###############################################################################
import os, os.path
import math
import cPickle as pickle
import numpy
from scipy import maxentropy, integrate, special
import scipy.interpolate
#from extreme_deconvolution import extreme_deconvolution
#import xdtarget
import isodist, isodist.imf
try:
    from galpy.util import bovy_plot
    _BOVY_PLOT_LOADED= True
except ImportError:
    _BOVY_PLOT_LOADED= False   
import dens_kde
def jkzcut(jk,upper=False):
    """Return the cut in jk-Z"""
    if upper:
        alpha= 4.
        x= 0.4
        A= 0.022/((0.6-x)**alpha-(0.5-x)**alpha)
        B= 0.03-A*(0.6-x)**alpha
        return A*(jk-x)**alpha+B
    else:
        alpha= 5.
        x= 0.225
        A= 0.028/((0.73-x)**alpha-(0.5-x)**alpha)
        B= 0.03-A*(0.73-x)**alpha
        return A*(jk-x)**alpha+B
class rcdist:
    """Class that holds the RC mean mag"""
    def __init__(self,*args,**kwargs):
        """
        NAME:
           __init__
        PURPOSE:
           initialize rcdist
        INPUT:
           Either:
              - file that holds a pickle
              - 2D-array [jk,Z], jks, Zs
        OUTPUT:
           object
        HISTORY:
           2012-11-15 - Written - Bovy (IAS)
        """
        if isinstance(args[0],str):
            if os.path.exists(args[0]):
                savefile= open(args[0],'rb')
                self._meanmag= pickle.load(savefile)
                self._jks= pickle.load(savefile)
                self._zs= pickle.load(savefile)
                savefile.close()
            else:
                raise IOError(args[0]+' file does not exist')
        else:
            self._meanmag= args[0]
            self._jks= args[1]
            self._zs= args[2]
        #Interpolate
        self._interpMag= scipy.interpolate.RectBivariateSpline(self._jks,
                                                               self._zs,
                                                               self._meanmag,
                                                               kx=3,ky=3,s=0.)
        return None

    def __call__(self,jk,Z,appmag=None,dk=0.036):
        """
        NAME:
           __call__
        PURPOSE:
           calls 
        INPUT:
           jk - color
           Z - metal-content
           appmag - apparent magnitude
           dk= calibration offset (dm= m-M-dk)
        OUTPUT:
           Either:
              - absmag (if appmag is None)
              - distance in kpc (if appmag given)
        HISTORY:
           2012-11-15 - Written - Bovy (IAS)
        """
        if appmag is None:
            return self._interpMag.ev(jk,Z)+dk
        else:
            absmag= self._interpMag.ev(jk,Z)
            return 10.**((appmag-absmag-dk)/5-2.)
    
class rcmodel:
    """rcmodel: isochrone model for the distribution in (J-Ks,M_H) along the RC"""
    def __init__(self,imfmodel='lognormalChabrier2001',Z=None,
                 interpolate=False,expsfh=False,band='H',
                 dontgather=False,loggmin=None,loggmax=None,
                 basti=False,ngauss=None,fitlogg=False,
                 parsec=False,stage=None,
                 loofrac=0.1):
        """
        NAME:
           __init__
        PURPOSE:
           initialize rcmodel
        INPUT:
           Z= metallicity (if not set, use flat prior in Z over all Z; can be list)
           imfmodel= (default: 'lognormalChabrier2001') IMF model to use (see isodist.imf code for options)
           band= band to use for M_X (JHK)
           interpolate= (default: False) if True, interpolate the binned representation
           expsfh= if True, use an exponentially-declining star-formation history
           dontgather= if True, don't gather surrounding Zs
           loggmin= if set, cut logg at this minimum
           loggmax= if set, cut logg at this maximum
           basti= if True, use Basti isochrones (if False, use Padova)
           parsec= if True, use PARSEC isochrones
           stage= if True, only use this evolutionary stage
           ngauss= number of Gaussians for XD, if -1, determine best
           loofrac= (0.1) fraction to loo if determining best
           fitlogg= if True, also fit logg in XD
        OUTPUT:
           object
        HISTORY:
           2012-11-07 - Written - Bovy (IAS)
        """
        self._band= band
        self._fitlogg= fitlogg
        self._ngauss= ngauss
        self._loggmin= loggmin
        self._loggmax= loggmax
        self._expsfh= expsfh
        self._Z= Z
        self._imfmodel= imfmodel
        self._basti= basti
        #Read isochrones
        if basti:
            zs= numpy.array([0.0001,0.0003,0.0006,0.001,0.002,0.004,0.008,
                             0.01,0.0198,0.03,0.04])
        else:
            zs= numpy.arange(0.0005,0.03005,0.0005)
        if Z is None:
            Zs= zs
        elif isinstance(Z,float):
            if basti or dontgather:
                Zs= [Z]
            elif Z < 0.001 or Z > 0.0295:
                Zs= [Z] 
            elif Z < 0.0015 or Z > 0.029:
                Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
            elif Z < 0.01:
                Zs= [Z-0.001,Z-0.0005,Z,Z+0.0005,Z+0.001] #build up statistics
            else:
                Zs= [Z-0.0005,Z,Z+0.0005] #build up statistics
        if basti:
            p= isodist.BastiIsochrone(Z=Zs)
        else:
            p= isodist.PadovaIsochrone(Z=Zs,parsec=parsec)
        maxage= 9.+numpy.log10(10.) #BaSTI goes too old, so does Padova
        if basti:
            #Setup KDE for logage distribution
            #BaSTI is not sampled uniformly in logage, so we empirically 
            #determine the logage distribution
            this_logages= p.logages()[(p.logages() < maxage)]
            this_logages= numpy.reshape(this_logages,(len(this_logages),1))
            self._bastikde= dens_kde.densKDE(this_logages,
                                             h=0.35,kernel='biweight')
        #Get relevant data
        sample= []
        weights= []
        loggs= []
        for logage in p.logages():
            if logage > maxage: continue
            for zz in range(len(Zs)):
                thisiso= p(logage,Zs[zz],asrecarray=True,stage=stage)
                if len(thisiso.M_ini) == 0: continue
                #Calculate int_IMF for this IMF model
                if not imfmodel == 'lognormalChabrier2001': #That would be the default
                    if imfmodel == 'exponentialChabrier2001':
                        int_IMF= isodist.imf.exponentialChabrier2001(thisiso.M_ini,int=True)
                    elif imfmodel == 'kroupa2003':
                        int_IMF= isodist.imf.kroupa2003(thisiso.M_ini,int=True)
                    elif imfmodel == 'chabrier2003':
                        int_IMF= isodist.imf.chabrier2003(thisiso.M_ini,int=True)
                    else:
                        raise IOError("imfmodel option not understood (non-existing model)")
                elif basti:
                    int_IMF= isodist.imf.lognormalChabrier2001(thisiso.M_ini,int=True)
                else:
                    int_IMF= thisiso.int_IMF
                dN= numpy.roll(int_IMF,-1)-int_IMF
                for ii in range(1,len(int_IMF)-1):
                    if basti:
                        JK= 0.996*(thisiso.J[ii]-thisiso.K[ii])+0.00923
                        #JK= 0.972*(thisiso.J[ii]-thisiso.K[ii])-0.011
                    else:
                        JK= thisiso.J[ii]-thisiso.Ks[ii]
                    if band.lower() == 'h':
                        if basti:
                            raise NotImplementedError("'H' not implemented for BaSTI yet")
                            J= JK+thisiso.K[ii]-0.044
                            H= J-(0.980*(thisiso.J[ii]-thisiso.H[ii])-0.045)
                        else:
                            H= thisiso.H[ii]
                    elif band.lower() == 'j':
                        if basti:
                            raise NotImplementedError("'J' not implemented for BaSTI yet")
                            J= JK+thisiso.K[ii]-0.044
                        else:
                            H= thisiso.J[ii]
                    elif band.lower() == 'k' or band.lower() == 'ks':
                        if basti:
                            H= thisiso.K[ii]-0.046
                        else:
                            H= thisiso.Ks[ii]
                    if JK < 0.3 \
                            or (not loggmax is None and thisiso['logg'][ii] > loggmax) \
                            or (not loggmin is None and thisiso['logg'][ii] < loggmin):
                        continue
                    if dN[ii] > 0.: 
                        sample.append([JK,H])
                        loggs.append([thisiso.logg[ii]])
                        if basti:
                            if expsfh:
                                weights.append(numpy.exp(self._bastikde(logage,log=True))*dN[ii]*10**(logage-7.)*numpy.exp((10.**(logage-7.))/800.)) #e.g., Binney (2010)
                            else:
                                weights.append(numpy.exp(self._bastikde(logage,log=True))*dN[ii]*10**(logage-7.))
                        else:
                            if expsfh:
                                weights.append(dN[ii]*10**(logage-7.)*numpy.exp((10.**(logage-7.))/800.)) #e.g., Binney (2010)
                            else:
                                weights.append(dN[ii]*10**(logage-7.))
                    else: 
                        continue #no use in continuing here   
        #Form array
        sample= numpy.array(sample)
        loggs= numpy.array(loggs)
        weights= numpy.array(weights)
        #Cut out low weights
        if False:
            indx= (weights > 10.**-5.*numpy.sum(weights))
        else:
            indx= numpy.ones(len(weights),dtype='bool')
        self._sample= sample[indx,:]
        self._weights= weights[indx]
        self._loggs= loggs[indx]
        #Setup KDE
        self._kde= dens_kde.densKDE(self._sample,w=self._weights,
                                    h='scott',kernel='biweight',
                                    variable=True,variablenitt=3,
                                    variableexp=0.35)
        self._jkmin, self._jkmax= 0.5,0.75
        self._hmin, self._hmax= -3.,0.
        return None
        #Run XD
        if ngauss > 0:
            l, xamp, xmean, xcovar= self._run_xd(ngauss=ngauss,fitlogg=fitlogg)
        else:
            ngauss, xamp, xmean, xcovar= self._determine_ngauss(fitlogg=fitlogg,
                                                                loofrac=loofrac)
            self._ngauss= ngauss
            l, xamp, xmean, xcovar= self._run_xd(ngauss=ngauss,fitlogg=fitlogg,
                                                 _xamp=xamp,
                                                 _xmean=xmean,
                                                 _xcovar=xcovar,
                                                 likeonly=True)
        self._loglike= l
        self._amp= xamp
        self._mean= xmean
        self._covar= xcovar
        self._xdtarg= xdtarget.xdtarget(amp=xamp,mean=xmean,covar=xcovar)
#        return None
        #Histogram
        self._nbins= 101#26#49
        self._djk= (self._jkmax-self._jkmin)/float(self._nbins)
        self._dh= (self._hmax-self._hmin)/float(self._nbins)
        self._hist, self._edges= numpy.histogramdd(sample,weights=weights,
                                                   bins=self._nbins,
                                                   range=[[self._jkmin,self._jkmax],
                                                          [self._hmin,self._hmax]])
        #Save
        self._Zs= Zs
        self._interpolate= interpolate
        self._loghist= numpy.log(self._hist)
        self._loghist[(self._hist == 0.)]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        if interpolate:
            #Form histogram grid
            jks= numpy.linspace(self._jkmin+self._djk/2.,
                                self._jkmax-self._djk/2.,
                                self._nbins)
            hs= numpy.linspace(self._hmin+self._dh/2.,
                               self._hmax-self._dh/2.,
                               self._nbins)
            self._interpolatedhist= scipy.interpolate.RectBivariateSpline(jks,hs,self._hist,
                                                                          bbox=[self._jkmin,self._jkmax,
                                            self._hmin,self._hmax],
                                                                          s=3.)
        return None

    def __call__(self,jk,h,sjk=None):
        """
        NAME:
           __call__
        PURPOSE:
           calls 
        INPUT:
           jk - color
           h - magnitude (set by 'band' in __init__)
           sjk= if not None, error
        OUTPUT:
           -
        HISTORY:
           2012-11-07 - Skeleton - Bovy (IAS)
        """
        return self.logpjkh(jk,h,sjk=sjk)

    def logpjkh(self,jk,h,sjk=None):
        """
        NAME:
           logpjkh
        PURPOSE:
           return the probability of the (J-Ks,M_H) pair
        INPUT:
           jk - J-Ks
           h - M_H (absolute magnitude)
           sjk= if not None, error
        OUTPUT:
           log of the probability
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if isinstance(jk,(list,numpy.ndarray)):
            testxs= numpy.empty((len(jk),2))
            testxs[:,0]= jk
            testxs[:,1]= h
            if not sjk is None:
                sjk2= numpy.zeros((len(jk),2))
                sjk2[:,0]= sjk**2.
            else:
                sjk2= None
            return self._kde(testxs,log=True,sx2=sjk2)
        else:
            if not sjk is None:
                sjk2= numpy.reshape(numpy.array([sjk**2.,0.]),
                                    (1,2))
            else:
                sjk2= None
            return self._kde(numpy.reshape(numpy.array([jk,h]),
                                           (1,2)),log=True,sx2=sjk2)
        if self._interpolate:
            return numpy.log(self._interpolatedhist(jk,h))
        else:
            jkbin= numpy.floor((jk-self._jkmin)/self._djk)
            hbin= numpy.floor((h-self._hmin)/self._dh)
            if isinstance(jk,numpy.ndarray):
                out= numpy.zeros(len(jk))-numpy.finfo(numpy.dtype(numpy.float64)).max
                jkbin= jkbin.astype('int')
                hbin= hbin.astype('int')
                indx= (jkbin >= 0.)*(hbin >= 0.)*(jkbin < self._nbins)\
                    *(hbin < self._nbins)
                out[indx]= self._loghist[jkbin[indx],hbin[indx]]
                return out
            else:
                jkbin= int(jkbin)
                hbin= int(hbin)
                if jkbin < 0 or jkbin >= self._nbins:
                    return -numpy.finfo(numpy.dtype(numpy.float64)).max
                if hbin < 0 or hbin >= self._nbins:
                    return -numpy.finfo(numpy.dtype(numpy.float64)).max
                return self._loghist[jkbin,hbin]

    def plot_pdf(self,jk,sjk=None,**kwargs):
        """
        NAME:
           plot_pdf
        PURPOSE:
           plot the conditioned PDF
        INPUT:
           jk - J-Ks
           sjk - error in J-K
           +bovy_plot.bovy_plot kwargs
        OUTPUT:
           plot to output
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("galpy.util.bovy_plot could not be imported")
        xs, lnpdf= self.calc_pdf(jk,sjk=sjk)
        if self._band == 'J':
            xlabel= r'$M_J$'
            xlim=[0.,-2.]
        elif self._band == 'H':
            xlabel= r'$M_H$'
            xlim=[0.,-2.]
        elif self._band == 'K':
            xlabel= r'$M_K$'
            xlim=[0.,-3.]
        elif self._band == 'Ks':
            xlabel= r'$M_{K_s}$'
            xlim=[0.,-2.]
        return bovy_plot.bovy_plot(xs,numpy.exp(lnpdf),'k-',
                                   xrange=xlim,
                                   yrange=[0.,
                                           1.1*numpy.amax(numpy.exp(lnpdf))],
                                   xlabel=xlabel,**kwargs)       

    def calc_pdf(self,jk,sjk=None,nxs=1001):
        """
        NAME:
           calc_pdf
        PURPOSE:
           calculate the conditioned PDF
        INPUT:
           jk - J-Ks
           sjk - error in J-K
           nxs= number of M_X
        OUTPUT:
           (xs,lnpdf)
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First collapse
        #coamp, comean, cocovar= self.collapse(jk,sjk=sjk)
        #Calculate pdf
        xs= numpy.linspace(-3.,0.,nxs)
        lnpdf= self(jk*numpy.ones(nxs),xs)
#        for ii, x in enumerate(xs):
#            this_lnpdf= numpy.log(coamp)-numpy.log(cocovar)-0.5*(x-comean)**2./cocovar
#            lnpdf[ii]=maxentropy.logsumexp(this_lnpdf)
#            lnpdf[ii]= self(jk,x)
        lnpdf[numpy.isnan(lnpdf)]= -numpy.finfo(numpy.dtype(numpy.float64)).max
        lnpdf-= maxentropy.logsumexp(lnpdf)+numpy.log(xs[1]-xs[0])
        return (xs,lnpdf)
    
    def calc_invcumul(self,jk,sjk=None):
        """
        NAME:
           calc_invcumul
        PURPOSE:
           calculate the inverse of the cumulative distribution
        INPUT:
           jk - J-Ks
           sjk - error in J-K
        OUTPUT:
           interpolated inverse cumulative distribution
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,sjk=sjk,nxs=1001)
        pdf= numpy.exp(lnpdf)
        pdf= numpy.cumsum(pdf)
        pdf/= pdf[-1]
        return scipy.interpolate.InterpolatedUnivariateSpline(pdf,xs,k=3)

    def median(self,jk,sjk=None):
        """
        NAME:
           median
        PURPOSE:
           return the median of the M_x distribution at this J-K
        INPUT:
           jk - J-Ks
           sjk - error in J-K
        OUTPUT:
           median
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the inverse cumulative distribution
        interpInvCumul= self.calc_invcumul(jk,sjk=sjk)
        return interpInvCumul(0.5)

    def quant(self,q,jk,sjk=None,sigma=True):
        """
        NAME:
           quant
        PURPOSE:
           return the quantile of the M_x distribution at this J-K
        INPUT:
           q - desired quantile in terms of 'sigma'
           jk - J-Ks
           sjk - error in J-K
           sigma= if False, the quantile is the actual quantile
        OUTPUT:
           quantile
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the inverse cumulative distribution
        interpInvCumul= self.calc_invcumul(jk,sjk=sjk)
        if not sigma:
            return interpInvCumul(q)
        else:
            if q > 0.:
                arg= 1.-(1.-special.erf(q/numpy.sqrt(2.)))/2.
            else:
                arg= (1.-special.erf(-q/numpy.sqrt(2.)))/2.
            return interpInvCumul(arg)

    def mode(self,jk,sjk=None):
        """
        NAME:
           mode
        PURPOSE:
           return the moden of the M_x distribution at this J-K
        INPUT:
           jk - J-Ks
           sjk - error in J-K
        OUTPUT:
           mode
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,sjk=sjk,nxs=1001)
        return xs[numpy.argmax(lnpdf)]
    
    def sigmafwhm(self,jk,sjk=None,straight=False):
        """
        NAME:
           sigmafwhm
        PURPOSE:
           return the sigma of the M_X distribution based on the FWHM
        INPUT:
           jk - J-Ks
           sjk - error in J-K
           straight= (False) if True, return actual hm points
        OUTPUT:
           FWHM/2.35...
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,sjk=None,nxs=1001)
        tmode= xs[numpy.argmax(lnpdf)]
        lnpdf_mode= numpy.amax(lnpdf)
        lnpdf_hm= lnpdf_mode-numpy.log(2.)
        minxs= xs[(xs < tmode)]
        minlnpdf= lnpdf[(xs < tmode)]
        minhm= minxs[numpy.argmin((minlnpdf-lnpdf_hm)**2.)]
        maxxs= xs[(xs > tmode)]
        maxlnpdf= lnpdf[(xs > tmode)]
        maxhm= maxxs[numpy.argmin((maxlnpdf-lnpdf_hm)**2.)]
        if straight:
            return (minhm,maxhm)
        else:
            return (maxhm-minhm)/2.*numpy.sqrt(2.*numpy.log(2.))
    
    def sigma2sigma(self,jk,sjk=None):
        """
        NAME:
           sigma2sigma
        PURPOSE:
           return the sigma obtained by integrating out to 2 sigma
        INPUT:
           jk - J-Ks
           sjk - error in J-K
        OUTPUT:
           2 sigma
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First calculate the PDF
        xs, lnpdf= self.calc_pdf(jk,sjk=sjk,nxs=1001)
        pdf= numpy.exp(lnpdf)
        cpdf= numpy.cumsum(pdf)
        cpdf/= cpdf[-1]
        q= 2.
        indx= (cpdf >= (1.-special.erf(q/numpy.sqrt(2.)))/2.)\
            *(cpdf <= (1.-(1.-special.erf(q/numpy.sqrt(2.)))/2.))
        m= numpy.sum(pdf[indx]*xs[indx])/numpy.sum(pdf[indx])
        return numpy.sqrt(numpy.sum(pdf[indx]*xs[indx]**2.)/numpy.sum(pdf[indx])-m**2.)/0.773741 #this factor to get a 'Gaussian' sigma
    
    def mean(self,jk,sjk=None):
        """
        NAME:
           mean
        PURPOSE:
           return the mean of the M_x distribution at this J-K
        INPUT:
           jk - J-Ks
           sjk - error in J-K
        OUTPUT:
           mean
        HISTORY:
           2012-11-07 - Written - Bovy (IAS)
        """
        #First collapse
        coamp, comean, cocovar= self.collapse(jk,sjk=sjk)
        #Then calculate mean
        return numpy.sum(coamp*comean)
    
    def sigma(self,jk,sjk=None):
        """
        NAME:
           sigma
        PURPOSE:
           return the standard dev of the M_x distribution at this J-K
           sjk - error in J-K
        INPUT:
           jk - J-Ks
        OUTPUT:
           standard deviation
        HISTORY:
           2012-11-07 - Written - Bovy (IAS)
        """
        #First collapse
        coamp, comean, cocovar= self.collapse(jk,sjk=sjk)
        #Then calculate mean
        tmean= numpy.sum(coamp*comean)
        #Then calculate squared
        tsq= numpy.sum(coamp*(cocovar+comean**2.))
        return numpy.sqrt(tsq-tmean**2.)

    def skew(self,jk,sjk=None):
        """
        NAME:
           skew
        PURPOSE:
           return the skew dev of the M_x distribution at this J-K
           sjk - error in J-K
        INPUT:
           jk - J-Ks
        OUTPUT:
           skew
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First collapse
        coamp, comean, cocovar= self.collapse(jk,sjk=sjk)
        #Then calculate mean
        tmean= numpy.sum(coamp*comean)
        #Then calculate squared
        tsq= numpy.sum(coamp*(cocovar+comean**2.))
        tsigma= numpy.sqrt(tsq-tmean**2.)
        #then calculate 3rd moment
        tt= numpy.sum(coamp*(3.*comean*cocovar+comean**3.))
        return (tt-3.*tmean*tsigma**2.-tmean**3.)/tsigma**3.

    def kurtosis(self,jk,sjk=None):
        """
        NAME:
           kurtosis
        PURPOSE:
           return the kurtosis of the M_x distribution at this J-K
           sjk - error in J-K
        INPUT:
           jk - J-Ks
        OUTPUT:
           kurtosis
        HISTORY:
           2012-11-09 - Written - Bovy (IAS)
        """
        #First collapse
        coamp, comean, cocovar= self.collapse(jk,sjk=sjk)
        #Then calculate mean
        tmean= numpy.sum(coamp*comean)
        #Then calculate squared
        tsq= numpy.sum(coamp*(cocovar+comean**2.))
        tsigma= numpy.sqrt(tsq-tmean**2.)
        #then calculate 3rd moment
        tt= numpy.sum(coamp*(3.*comean*cocovar+comean**3.))
        #then calculate 4th moment
        tth= numpy.sum(coamp*(3.*cocovar**2.+6.*comean**2.*cocovar+comean**4.))
        return (tth-4.*tmean*tt+6.*tmean**2.*tsigma**2.+3.*tmean**4.)/tsigma**4.-3.

    def collapse(self,jk,sjk=None):
        """
        NAME:
           collapse
        PURPOSE:
           collapse the distribution onto this jk (condition!)
        INPUT:
           jk - J-Ks
           sjk - error in J-K
        OUTPUT:
           amp,mean,var
        HISTORY:
           2012-11-08 - Written - Bovy (IAS)
        """
        if self._fitlogg:
            raise NotImplementedError("Collapse for fitlogg not implemented yet")
        #First calculate amplitudes
        logrelamp= numpy.log(self._amp)-0.5*numpy.log(self._covar[:,0,0])\
            -0.5*(jk-self._mean[:,0])**2./self._covar[:,0,0]
        logrelamp-= maxentropy.logsumexp(logrelamp)
        #Then calculate conditioned means and variances
        comean= self._mean[:,1]+self._covar[:,0,1]/(self._covar[:,0,0]+sjk**2.)*(jk-self._mean[:,0])
        cocovar= self._covar[:,1,1]-self._covar[:,0,1]**2./(self._covar[:,0,0]+sjk**2.)
        return (numpy.exp(logrelamp),comean,cocovar)

    def _determine_ngauss(self,fitlogg=False,loofrac=0.1):
        """Determine the optimal number of Gaussians"""
        ngauss_s= range(1,15)
        #First do loo
        perm= numpy.random.permutation(len(self._weights))
        tperm= perm[0:int(numpy.ceil(loofrac*len(self._weights)))]
        tsample= self._sample[tperm,:]
        tweights= self._weights[tperm]
        tloggs= self._loggs[tperm]
        #loo sample
        lperm= perm[int(numpy.ceil(loofrac*len(self._weights))):len(self._weights)]
        lsample= self._sample[lperm,:]
        lweights= self._weights[lperm]
        lloggs= self._loggs[lperm]
        #!!LOO!!
        loos= numpy.zeros(len(ngauss_s))
        xamps= []
        xmeans= []
        xcovars= []
        for ii, ngauss in enumerate(ngauss_s):
            print ngauss
            #First run XD
            l, xamp, xmean, xcovar= self._run_xd(ngauss=ngauss,fitlogg=fitlogg,
                                                 _sample=tsample,
                                                 _weights=tweights,
                                                 _loggs=tloggs)
            tl, xamp, xmean, xcovar= self._run_xd(ngauss=ngauss,
                                                  fitlogg=fitlogg,
                                                  likeonly=True,
                                                 _sample=lsample,
                                                 _weights=lweights,
                                                 _loggs=lloggs,
                                                  _xamp=xamp,
                                                  _xmean=xmean,
                                                  _xcovar=xcovar)
            loos[ii]= tl*len(lweights)
            xamps.append(xamp)
            xmeans.append(xmean)
            xcovars.append(xcovar)
        bovy_plot.bovy_plot(ngauss_s,loos,'ko',
                            yrange=[0.9*numpy.amin(loos),1.1*numpy.amax(loos)],
                            xlabel=r'$\#\ \mathrm{Gaussian\ components}$',
                            ylabel=r'$\mathrm{loo\ log\ likelihood}$')
        indx= numpy.argmax(loos)
        print ngauss_s[indx]
        return (ngauss_s[numpy.argmax(loos)],
                xamps[indx],
                xmeans[indx],
                xcovars[indx])

    def _run_xd(self,ngauss=None,fitlogg=False,
                _weights=None,_sample=None,_loggs=None,
                likeonly=False,
                _xamp=None,_xmean=None,_xcovar=None):
        """Run XD"""
        if _weights == None:
            _weights= self._weights
        if _sample == None:
            _sample= self._sample
        if _loggs == None:
            _loggs= self._loggs
        if fitlogg:
            dim= 3
        else:
            dim= 2
        ydata= numpy.reshape(_sample,(len(_weights),2))
        if fitlogg:
            tmp_ydata= numpy.zeros((len(_weights),dim))
            tmp_ydata[:,0]= ydata[:,0]
            tmp_ydata[:,1]= ydata[:,1]
            tmp_ydata[:,2]= logit(_loggs,
                                  self._loggmin-0.2,self._loggmax+0.2)[:,0]
            ydata= tmp_ydata
        ycovar= numpy.zeros((len(_weights),dim))
        if _xamp is None:
            xamp= numpy.ones(ngauss)/float(ngauss)
        else:
            xamp= _xamp
        datamean= numpy.mean(ydata,axis=0)
        datacov= numpy.cov(ydata,rowvar=0)
        datastd= numpy.sqrt(numpy.diag(datacov))
        if _xmean is None:
            xmean= datamean+numpy.random.normal(size=(ngauss,dim))\
                *numpy.tile(datastd,(ngauss,1))
        else:
            xmean= _xmean
        if _xcovar is None:
            xcovar= numpy.tile(datacov,(ngauss,1,1))*4.
        else:
            xcovar= _xcovar
        l= extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,
                                 weight=_weights,likeonly=likeonly,
                                 splitnmerge=0,w=10.**-6.)
        return (l, xamp, xmean, xcovar)

    def plot(self,log=False,conditional=False,nbins=None,sjk=None):
        """
        NAME:
           plot
        PURPOSE:
           plot the resulting (J-Ks,H) distribution
        INPUT:
           log= (default: False) if True, plot log
           conditional= (default: False) if True, plot conditional distribution
                        of H given J-Ks
           nbins= if set, set the number of bins
           sjk= error
        OUTPUT:
           plot to output device
        HISTORY:
           2012-02-17 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("'galpy.util.bovy_plot' plotting package not found")
        #Form histogram grid
        if nbins is None:
            nbins= self._nbins
        djk= (self._jkmax-self._jkmin)/float(nbins)
        dh= (self._hmax-self._hmin)/float(nbins)
        jks= numpy.linspace(self._jkmin+djk/2.,
                            self._jkmax-djk/2.,
                            nbins)
        hs= numpy.linspace(self._hmax-dh/2.,#we reverse
                           self._hmin+dh/2.,
                           nbins)
        plotthis= numpy.zeros((nbins,nbins))
        for ii in range(nbins):
            for jj in range(nbins):
                plotthis[ii,jj]= self(jks[ii],hs[jj],sjk=sjk)
        if not log:
            plotthis= numpy.exp(plotthis)
        if conditional: #normalize further
            for ii in range(nbins):
                plotthis[ii,:]/= numpy.nanmax(plotthis[ii,:])/numpy.nanmax(plotthis)
        if self._band == 'J':
            ylabel= r'$M_J$'
            ylim=[0.,-3.]
        elif self._band == 'H':
            ylabel= r'$M_H$'
            ylim=[0.,-3.]
        elif self._band == 'K':
            ylabel= r'$M_K$'
            ylim=[0.,-3.]
        elif self._band == 'Ks':
            ylabel= r'$M_{K_s}$'
            ylim=[0.,-3.]
        return bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='gist_yarg',
                                     xrange=[self._jkmin,self._jkmax],
                                     yrange=ylim,
                                     aspect=(self._jkmax-self._jkmin)/(self._hmax-self._hmin),
                                     xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                                     ylabel=ylabel,
                                     interpolation='nearest')
    
    def plot_samples(self):
        """
        NAME:
           plot_samples
        PURPOSE:
           plot the samples that the histogramming is based on 
        INPUT:
        OUTPUT:
           plot to output device
        HISTORY:
           2012-06-15 - Written - Bovy (IAS)
        """
        if not _BOVY_PLOT_LOADED:
            raise ImportError("'galpy.util.bovy_plot' plotting package not found")
        if self._band == 'J':
            ylabel= r'$M_J$'
            ylim=[0.,-3.]
        elif self._band == 'H':
            ylabel= r'$M_H$'
            ylim=[0.,-3.]
        elif self._band == 'K':
            ylabel= r'$M_K$'
            ylim=[0.,-3.]
        elif self._band == 'Ks':
            ylabel= r'$M_{K_s}$'
            ylim=[0.,-3.]
        return bovy_plot.bovy_plot(self._sample[:,0],self._sample[:,1],
                                   xrange=[0.5,0.75],
                                   yrange=ylim,
                                   xlabel=r'$(J-K_s)_0\ [\mathrm{mag}]$',
                                   ylabel=ylabel,
                                   scatter=True,
                                   c=self._weights,
                                   edgecolors='none',
                                   s=numpy.ones(len(self._weights))*10.,
                                   alpha=0.5,
                                   marker='o',
                                   colorbar=True)
    
    def plot_model(self,nsamples=1000):
        """
        NAME:
           plot_model
        PURPOSE:
           plot the model distribution as a a set of samples and Gaussians
        INPUT:
        OUTPUT:
           plot to output device
        HISTORY:
           2012-07-08 - Written - Bovy (IAS)
        """
        self._xdtarg.sample(nsample=nsamples)
        if self._band == 'J':
            ylabel= r'$M_J$'
            ylim=[0.,-3.]
        elif self._band == 'H':
            ylabel= r'$M_H$'
            ylim=[0.,-3.]
        elif self._band == 'K':
            ylabel= r'$M_K$'
            ylim=[0.,-3.]
        elif self._band == 'Ks':
            ylabel= r'$M_{K_s}$'
            ylim=[0.,-3.]
        if self._fitlogg:
            scatter= True
            plotc= inv_logit(numpy.array(self._xdtarg.samples[:,2]),
                             self._loggmin-0.2,self._loggmax+0.2)
            colorbar= True
            clabel= r'$\log g$'
            crange=[self._loggmin,self._loggmax]
            marker= 'o'
            color= None
            alpha=0.5
            s= numpy.ones_like(plotc)*10.,
        else:
            scatter= True
            plotc= 'k'
            colorbar= False
            clabel= None
            crange= [0.,0.]
            marker= ','
            color= 'k'
            alpha=1.
            s=1.
        return self._xdtarg.scatterplot(0,1,
                                        xlabel=r'$J-K_s$',
                                        ylabel=ylabel,
                                        xrange=[0.5,0.75],
                                        yrange=ylim,
                                        scatter=scatter,
                                        cmap='jet',
                                        color=color,
                                        c=plotc,
                                        marker=marker,
                                        s=s,
                                        colorbar=colorbar,
                                        clabel=clabel,
                                        edgecolors='none',
                                        vmin=crange[0],vmax=crange[1],
                                        crange=crange,
                                        alpha=alpha,
                                        hoggscatter=False)
def logit(x,xmin,xmax):
    p= (x-xmin)/(xmax-xmin)
    return numpy.log(p/(1.-p))

def inv_logit(p,xmin,xmax):
    x= numpy.exp(p)/(1.+numpy.exp(p))
    return x*(xmax-xmin)+xmin
