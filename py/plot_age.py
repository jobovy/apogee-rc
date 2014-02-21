import sys
import os, os.path
import pickle
import numpy
from scipy import interpolate, optimize
from galpy.util import bovy_plot
from calc_avg_rcmks import localzdist
_FIT= False
def plot_age(plotfilename,massfile,omegafile):
    if os.path.exists(massfile):
        savefile= open(massfile,'rb')
        mass= pickle.load(savefile)
        zs= pickle.load(savefile)
        lages= pickle.load(savefile)
        savefile.close()
    else:
        raise IOError("file %s has to exist ..." % massfile)
    if os.path.exists(omegafile):
        savefile= open(omegafile,'rb')
        omega= pickle.load(savefile)
        savefile.close()
    else:
        raise IOError("file %s has to exist ..." % omegafile)
    agezdist= omega/mass
    pz= numpy.exp(numpy.array([localzdist(z,zsolar=0.017) for z in zs]))
    pz/= numpy.sum(pz)
    page= 1.
    #Tile
    pz= numpy.tile(pz,(len(lages),1)).T
    #Build age pdfs
    postage= page*numpy.nansum(pz*agezdist,axis=0)/numpy.nansum(pz,axis=0)/10.**lages
    postage= postage[lages > numpy.log10(0.8)]
    postage/= numpy.nanmax(postage)
    #Interpolate
    lages= lages[lages > numpy.log10(0.8)]
    postage_spline= interpolate.InterpolatedUnivariateSpline(lages,
                                                             numpy.log(postage),
                                                             k=3)
    plages= numpy.linspace(0.8,10.,1001)
    plpostage= numpy.exp(postage_spline(numpy.log10(plages)))
    plpostage/= numpy.nansum(plpostage)*(plages[1]-plages[0])
    bovy_plot.bovy_print(fig_width=7.)
    bovy_plot.bovy_plot(plages,plpostage,
                        'k-',lw=2.,
                        xlabel=r'$\mathrm{Age}\,(\mathrm{Gyr})$',
                        ylabel=r'$p(\mathrm{RC | population\ Age})$',
                        xrange=[0.,10.],
                        yrange=[0.,0.4],
                        zorder=10,loglog=False)
    #Baseline
    bovy_plot.bovy_plot(plages,1./9.2*numpy.ones(len(plages)),
                        '-',color='0.4',overplot=True,zorder=3,lw=2.)
    #Fit the dependence
    if _FIT:
        opt= optimize.fmin_powell(chi2,[-1.,0.,0.,0.,-1.,0.,0.,0.],#[-4.,0.,0.25,0.,0.,1.],
                                  args=(numpy.log10(plages),
                                        numpy.log(plpostage)))
        print opt
        bovy_plot.bovy_plot(plages,numpy.exp(_fit_func(opt,numpy.log10(plages))),
                            'r-',overplot=True,zorder=11)
        a= numpy.log10(plages)
        bovy_plot.bovy_plot(plages[a <= 0.23],numpy.exp(polyfit1(a[a <= 0.23])),'b-',overplot=True,zorder=12)
        bovy_plot.bovy_plot(plages[a > 0.23],numpy.exp(polyfit2(a[a > 0.23])),'b-',overplot=True,zorder=12)
    bovy_plot.bovy_end_print(plotfilename)
    return None

def chi2(p,a,pa):
    #relamp= p[2]
    #if relamp < 0.: return 100000000000.
    #elif relamp > 1.: return 100000000000.
    fpa= _fit_func(p,a)
    return 0.5*numpy.sum((pa-fpa)**2.)

def _fit_func(p,a):
    #pl= p[0]
    #pl2= p[1]
    #relamp= p[2]
    #cut= numpy.exp(p[3])
    #typa2= numpy.exp(p[4])
    np1= 4
    np2= 4
    pol= 0.
    pol2= 0.
    for ii in range(np1):
        pol+= p[ii]*a**ii
    for ii in range(np2):
        pol2+= p[ii+np1]*a**ii
    out= pol
    out[a > 0.23]= pol2[a > 0.23]
    return out
    return pol#*((1.-relamp)*a**pl+relamp*pl2*numpy.exp(-(a-typa2)/cut))

def polyfit1(a):
    return -1.6314+3.8230*a+2.2133*a**2.-35.7414*a**3.
def polyfit2(a):
    return -1.0666+1.8664*a-9.0617*a**2.+4.5860*a**3.

if __name__ == '__main__':
    plot_age(sys.argv[1],sys.argv[2],sys.argv[3])
