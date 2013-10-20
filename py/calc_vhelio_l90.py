import sys
import os, os.path
import pickle
import math
import numpy
from scipy import interpolate
import apogee.tools.read as apread
from galpy.orbit import Orbit
from galpy.df import dehnendf
from galpy.util import save_pickles, bovy_coords, bovy_plot
#Quick def
def safe_dl_to_rphi(d,l):
    R= math.sqrt(1.+d**2.-2.*d*math.cos(l))
    if R == 0.:
        R= 0.0001
        d+= 0.0001
    if 1./math.cos(l) < d and math.cos(l) > 0.:
        theta= math.pi-math.asin(d/R*math.sin(l))
    else:
        theta= math.asin(d/R*math.sin(l))
    return (R,theta,d,l)

if __name__ == '__main__':
    nds= 101
    ds= numpy.linspace(0.001,8./8.,nds)
    dfc= dehnendf(profileParams=(3./8.,2.,36./218.),beta=0.,correct=False)
    if False:
        #First calculate mean
        mean_savefilename= 'vhelio_l90.sav'
        if os.path.exists(mean_savefilename):
            savefile= open(mean_savefilename,'rb')
            plotds= pickle.load(savefile)
            plotvlos= pickle.load(savefile)
            savefile.close()
        else:
            vlos= numpy.zeros(nds)
            for jj in range(nds):
                print jj
                d,l= ds[jj], 90./180.*numpy.pi
                R,theta,d,l= safe_dl_to_rphi(d,l)
                vlos[jj]= dfc.meanvT(R)*math.sin(theta+l)
            plotds= ds*8.
            plotvlos= vlos*218.-218
            save_pickles(mean_savefilename,
                         plotds,plotvlos)
    #Now calculate the whole vlos distribution at each d
    dist_savefilename= 'vhelio_l90_dist.sav'
    if os.path.exists(dist_savefilename):
        savefile= open(dist_savefilename,'rb')
        plotds= pickle.load(savefile)
        plotvloss= pickle.load(savefile)
        vlos_dist= pickle.load(savefile)
        savefile.close()
    else:
        nvloss= 101
        vloss= numpy.linspace(-0.7,0.25,nvloss)
        vlos_dist= numpy.zeros((nds,nvloss))
        ra, dec= bovy_coords.lb_to_radec(90.,0.,degree=True)
        for ii in range(nds):
            print ii
            for jj in range(nvloss):
                #Setup this Orbit
                o= Orbit([ra,dec,ds[ii],0.,0.,vloss[jj]],
                         radec=True,vo=1.,ro=1.,zo=0.,
                         solarmotion=[0.,25./218.,0.])
                vlos_dist[ii,jj]= dfc(o,marginalizeVperp=True)
        plotds= ds*8.
        plotvloss= vloss*218.
        save_pickles(dist_savefilename,plotds,plotvloss,vlos_dist)
    #Now plot
    bovy_plot.bovy_print()
    for ii in range(nds):
        vlos_dist[ii,:]/= numpy.sum(vlos_dist[ii,:])
    bovy_plot.bovy_dens2d(vlos_dist.T,origin='lower',cmap='gist_yarg',
                          xrange=[0.001-(ds[1]-ds[0])/2.,
                                  8.+(ds[1]-ds[0])/2.],
                          yrange=[-0.7*218.-(plotvloss[1]-plotvloss[0])/2.,
                                   0.25*218.+(plotvloss[1]-plotvloss[0])/2.],
                          xlabel=r'$\mathrm{distance}\,(\mathrm{kpc})$',
                          ylabel=r'$\mathrm{line\!-\!of\!-\!sight\ velocity}\,(\mathrm{km\,s}^{-1})$')                         
    mid50= numpy.zeros(len(ds))
    up95= numpy.zeros(len(ds))
    down95= numpy.zeros(len(ds))
    for ii in range(len(ds)):
        tdist= numpy.cumsum(vlos_dist[ii,:])/numpy.sum(vlos_dist[ii,:])
        mid50[ii]= plotvloss[numpy.argmin(numpy.fabs(tdist-0.5))]
        down95[ii]= plotvloss[numpy.argmin(numpy.fabs(tdist-0.025))]
        up95[ii]= plotvloss[numpy.argmin(numpy.fabs(tdist-0.975))]
    interpmid50= interpolate.UnivariateSpline(ds,mid50,k=3)
    interpup95= interpolate.UnivariateSpline(ds,up95,k=3)
    interpdown95= interpolate.UnivariateSpline(ds,down95,k=3)
    #Plot median
    bovy_plot.bovy_plot(plotds,interpmid50(ds)+15.,
                        'w--',lw=2.,overplot=True)
    bovy_plot.bovy_plot(plotds,interpmid50(ds),'w-',lw=2.,overplot=True)
    #Plot 95%
    bovy_plot.bovy_plot(plotds,interpup95(ds),'-',color='0.3',overplot=True)
    bovy_plot.bovy_plot(plotds,interpdown95(ds),'-',color='0.3',overplot=True)
    data= apread.rcsample()
    data= data[(numpy.fabs(data['GLON']-90.) < 5.)*(numpy.fabs(data['RC_GALZ']) < 0.05)]
    print numpy.amin(data['GLON'])
    print numpy.amax(data['GLON'])
    print len(data)
    bovy_plot.bovy_plot(data['RC_DIST'],data['VHELIO_AVG'],
                        'o',mfc='.45',mec='0.45',
                        ms=5.,overplot=True,zorder=10)
    bovy_plot.bovy_text(r'$l=90^\circ,\, |Z| < 50\,\mathrm{pc}$',
                        bottom_left=True,size=18.)
    bovy_plot.bovy_text(r'$\mathrm{solar\ motion}\ =\ 25\,\mathrm{km\,s}^{-1}$'
                        +'\n'+r'$\mathrm{dashed}\ \rightarrow\ 10\,\mathrm{km\,s}^{-1}$',
                        top_right=True,size=18.)
    bovy_plot.bovy_end_print(sys.argv[1])
