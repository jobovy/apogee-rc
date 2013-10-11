import sys
import numpy
from scipy import optimize
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
def plot_sigmar(plotfilename,plotfilename2):
    #First read the sample
    data= apread.rcsample()
    data= data[(data['PMMATCH'] == 1)] #Good velocities
    galy= data['RC_GALR']*numpy.sin(data['RC_GALPHI'])
    indx= (numpy.fabs(galy) < .75)*(numpy.fabs(data['RC_GALZ']) < 0.05) #(numpy.fabs(data['GLAT']) < 1.51)
    data= data[indx]
    print len(data)
    rs= numpy.arange(3.5,16.501,1.)
    sigrs= numpy.zeros_like(rs)+numpy.nan
    mvrs= numpy.zeros_like(rs)+numpy.nan
    emvrs= numpy.zeros_like(rs)+numpy.nan
    esigrs= numpy.zeros_like(rs)+numpy.nan
    pmesigrs= numpy.zeros_like(rs)+numpy.nan
    for ii in range(len(rs)-1):
        tindx= (data['RC_GALR'] > rs[ii])\
            *(data['RC_GALR'] <= rs[ii+1])
        print rs[ii], numpy.sum(tindx)
        sigrs[ii]= numpy.std(data['GALVR'][tindx])
        mvrs[ii]= numpy.median(data['GALVR'][tindx])
        emvrs[ii]= sigrs[ii]/numpy.sqrt(numpy.sum(tindx))
        esigrs[ii]= sigrs[ii]/numpy.sqrt(2.*numpy.sum(tindx))
        pmesigrs[ii]= 4.74047*numpy.sqrt(numpy.mean(data['RC_DIST'][tindx]**2.\
                                                   *data['PMRA_ERR'][tindx]**2.\
                                                   *numpy.sin(data['GLON'][tindx]/180.*numpy.pi+data['RC_GALPHI'][tindx])**2.))
    #Now plot
    bovy_plot.bovy_print()
    print sigrs
    print esigrs
    bovy_plot.bovy_plot(rs+0.5*(rs[1]-rs[0]),sigrs,'ko',semilogy=True,
                        xlabel=r'$R\,(\mathrm{kpc})$',
                        ylabel=r'$\sigma_R\,(\mathrm{km\,s}^{-1})$',
                        xrange=[0.,18.],
                        yrange=[10.,100.])
    pyplot.errorbar(rs+0.5*(rs[1]-rs[0]),sigrs,yerr=esigrs,marker='None',
                    ls='none',color='k')
    #bovy_plot.bovy_plot(rs+0.5*(rs[1]-rs[0]),pmesigrs,'ks',overplot=True)
    #Fit
    indx= True-numpy.isnan(sigrs)
    #esigrs[3]= 100.
    bestfit= optimize.curve_fit(expfit,(rs+0.5*(rs[1]-rs[0]))[indx],
                                sigrs[indx],sigma=esigrs[indx],
                                p0=(numpy.log(32.),numpy.log(12.)))
    print numpy.exp(bestfit[0]), bestfit[1]
    bovy_plot.bovy_plot(rs+0.5*(rs[1]-rs[0]),
                        numpy.exp(bestfit[0][0])\
                            *numpy.exp(-(rs+0.5*(rs[1]-rs[0])-8.)/numpy.exp(bestfit[0][1])),
                        'k--',overplot=True)
    pyplot.yticks([10,20,30,40,50,60,70,80,90,100],[r'$10$',r'$20$',r'$30$',r'$40$',r'$50$',r'$60$',r'$70$',r'$80$',r'$90$',r'$100$'])
    bovy_plot.bovy_text(r'$|Z| < 50\,\mathrm{pc}$',top_right=True,size=18.)
    bovy_plot.bovy_text(r'$\sigma_R = %.0f\,\mathrm{km\,s}^{-1}\,\exp\left(-\left[R-8\,\mathrm{kpc}\right] / %.0f\,\mathrm{kpc}\right)$' % (numpy.exp(bestfit[0][0]),numpy.exp(bestfit[0][1])),
                        bottom_left=True,size=14.)
    bovy_plot.bovy_end_print(plotfilename)
    #Now plot vr
    #First read the sample
    data= apread.rcsample()
    data= data[(data['PMMATCH'] == 1)] #Good velocities
    galy= data['RC_GALR']*numpy.sin(data['RC_GALPHI'])
    indx= (numpy.fabs(galy) < .75)*(numpy.fabs(data['RC_GALZ']) < 0.05) #(numpy.fabs(data['GLAT']) < 1.51)
    data= data[indx]
    print len(data)
    rs= numpy.arange(3.5,16.501,.75)
    sigrs= numpy.zeros_like(rs)+numpy.nan
    mvrs= numpy.zeros_like(rs)+numpy.nan
    emvrs= numpy.zeros_like(rs)+numpy.nan
    esigrs= numpy.zeros_like(rs)+numpy.nan
    pmesigrs= numpy.zeros_like(rs)+numpy.nan
    for ii in range(len(rs)-1):
        tindx= (data['RC_GALR'] > rs[ii])\
            *(data['RC_GALR'] <= rs[ii+1])
        print rs[ii], numpy.sum(tindx)
        sigrs[ii]= numpy.std(data['GALVR'][tindx])
        mvrs[ii]= numpy.median(data['GALVR'][tindx])
        emvrs[ii]= sigrs[ii]/numpy.sqrt(numpy.sum(tindx))
        esigrs[ii]= sigrs[ii]/numpy.sqrt(2.*numpy.sum(tindx))
        pmesigrs[ii]= 4.74047*numpy.sqrt(numpy.mean(data['RC_DIST'][tindx]**2.\
                                                   *data['PMRA_ERR'][tindx]**2.\
                                                   *numpy.sin(data['GLON'][tindx]/180.*numpy.pi+data['RC_GALPHI'][tindx])**2.))
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(rs+0.5*(rs[1]-rs[0]),mvrs,'ko',
                        xlabel=r'$R\,(\mathrm{kpc})$',
                        ylabel=r'$\mathrm{median}\ v_R\,(\mathrm{km\,s}^{-1})$',
                        xrange=[0.,18.],
                        yrange=[-10.,10.])
    pyplot.errorbar(rs+0.5*(rs[1]-rs[0]),mvrs,yerr=emvrs,marker='None',
                    ls='none',color='k')
    bovy_plot.bovy_text(r'$|Z| < 50\,\mathrm{pc}$',top_right=True,size=18.)
    bovy_plot.bovy_end_print(plotfilename2)
    

def expfit(x,lnsro,lnhsr):
    return numpy.exp(lnsro)*numpy.exp(-(x-8.)/numpy.exp(lnhsr))

if __name__ == '__main__':
    plot_sigmar(sys.argv[1],sys.argv[2])
