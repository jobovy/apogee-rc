import sys
import os, os.path
import csv
import socket
import numpy
from scipy import interpolate
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter
import fitsio
import apogee.tools.read as apread
import pixelize_sample
import bovy_psd
from plot_2dkinematics import dvlosgal
import hackGCS
import readAndHackHoltz
from simple_spiral_simulation import simulate_vlos_spiral
import galpy_simulations
_EXT='png'
_ADDLLOGGCUT= True
_ADDGCS= True
_ADDRED= False
_ADDRAVE= True
_ADDRIX= True
_ADDSIMPLESPIRAL= True
_INTERP= False
_NNOISE= 1000
_PLOTBAND= False
_SIGNIF= 0.95
_SUBTRACTERRORS= 1.
_PROPOSAL= False
#Parameters of the pixelizations
#APOGEE-RC
_RCXMIN= 5.5
_RCXMAX= 12.25
_RCYMIN= -3.
_RCYMAX= 3.75
_RCDX= 0.75
#RAVE
_RAVEXMIN= 6.75
_RAVEXMAX= 8.75
_RAVEYMIN= -1.75
_RAVEYMAX= 0.25
_RAVEDX= 0.25
#GCS
_GCSXMIN= -0.0625
_GCSXMAX= 0.0625
_GCSYMIN= -0.0625
_GCSYMAX= 0.0625
_GCSDX= 0.025
# output to file
_DUMP2FILE= True
def plot_psd(plotfilename):
    data= apread.rcsample()
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #Cut
    indx= (numpy.fabs(data['RC_GALZ']) < 0.25)*(data['METALS'] > -1000.)
    print "Using %i stars for low-Z 2D kinematics analysis" % numpy.sum(indx)
    data= data[indx]
    #Get residuals
    dx= _RCDX
    binsize= .8#.765
    pix= pixelize_sample.pixelXY(data,
                                 xmin=_RCXMIN,xmax=_RCXMAX,
                                 ymin=_RCYMIN,ymax=_RCYMAX,
                                 dx=dx,dy=dx)
    resv= pix.plot(lambda x: dvlosgal(x,vtsun=220.+22.5),
                   returnz=True,justcalc=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-4))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][1:-3]
    scale= 4.*numpy.pi#3.25/numpy.median(numpy.sqrt(noisepsd))
    print scale
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-4))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][1:-3]
    ks= psd1d[0][1:-3]
    if _ADDGCS:
        xrange=[.03,110.]
    else:
        xrange= [0.,1.]
    yrange= [0.,11.9]
    if _PROPOSAL:
        bovy_plot.bovy_print(fig_width=7.5,fig_height=3.)
    else:
        bovy_plot.bovy_print(fig_width=7.5,fig_height=4.5)
    apop= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][1:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                        'ko',lw=2.,
                        zorder=12,
                        xlabel=r'$k\,(\mathrm{kpc}^{-1})$',
                        ylabel=r'$\sqrt{P_k}\,(\mathrm{km\,s}^{-1})$',
                        semilogx=_ADDGCS,
                        xrange=xrange,yrange=yrange)
    if _DUMP2FILE:
        with open('bovy-apogee-psd.dat','w') as csvfile:
            writer= csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            csvfile.write('#APOGEE\n')
            for ii in range(len(ks)):
                writer.writerow([ks[ii],
                                 (scale*numpy.sqrt(psd1d[1][1:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)))[ii],
                                 (scale*0.5*(psd1d[2][1:-3]**2.+_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)**2.)**0.5/numpy.sqrt(psd1d[1][1:-3]))[ii]])
    if _PROPOSAL:
        pyplot.gcf().subplots_adjust(bottom=0.15)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][1:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*(psd1d[2][1:-3]**2.
                                    +_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)**2.)**0.5/numpy.sqrt(psd1d[1][1:-3]),
                    marker='None',ls='none',color='k')
    if _PLOTBAND:
#        bovy_plot.bovy_plot(ks,
#                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
#                            '--',lw=2.,zorder=10,
#                            color='0.85',overplot=True)
#        bovy_plot.bovy_plot(ks,
#                            scale*numpy.sqrt(numpy.sort(noisepsd
#                                                        -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0),axis=0)[int(numpy.floor(0.99*_NNOISE)),:]),
#                            zorder=1,ls='--',overplot=True,
#                            color='0.65')
        pyplot.fill_between(ks,
                            scale*numpy.sqrt(numpy.sort(noisepsd
                                                        -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0),axis=0)[int(numpy.floor(_SIGNIF*_NNOISE)),:]),
                            zorder=0,
                            color='0.65')
    #Add a simple model of a spiral potential
    if _ADDSIMPLESPIRAL:
        if False:
            alpha= -12.5
            spvlos= simulate_vlos_spiral(alpha=alpha,
                                         gamma=1.2,
                                         xmin=_RCXMIN,xmax=_RCXMAX,
                                         ymin=_RCYMIN,ymax=_RCYMAX,
                                         dx=0.01)
            potscale= 1.35*1.2
            print numpy.arctan(2./alpha)/numpy.pi*180., numpy.sqrt(0.035/numpy.fabs(alpha)/2.)*potscale*220., numpy.sqrt(0.035/numpy.fabs(alpha))*potscale*220.
            simpsd1d= bovy_psd.psd1d(spvlos*220.*potscale,0.01,binsize=binsize)
            tks= simpsd1d[0][1:-3]
            bovy_plot.bovy_plot(tks,
                                scale*numpy.sqrt(simpsd1d[1][1:-3]),
                                'k--',lw=2.,overplot=True)
        #A better simulation
#        spvlos= galpy_simulations.vlos('../sim/spiral_rect_omegas0.33_alpha-14.sav')
        spvlos= galpy_simulations.vlos('../sim/bar_rect_alpha0.015_hivres.sav')
        potscale= 1.
        simpsd1d= bovy_psd.psd1d(spvlos*220.*potscale,0.33333333,binsize=binsize)
        tks= simpsd1d[0][1:-3]
#        alpha=-14.
#        print numpy.arctan(2./-14.)/numpy.pi*180., numpy.sqrt(0.075/numpy.fabs(alpha)/2.)*potscale*220., numpy.sqrt(0.075/numpy.fabs(alpha))*potscale*220.
        line1= bovy_plot.bovy_plot(tks,
                                   scale*numpy.sqrt(simpsd1d[1][1:-3]),
                                   'k--',lw=2.,overplot=True)
        if _DUMP2FILE:
            with open('bovy-bar-psd.dat','w') as csvfile:
                writer= csv.writer(csvfile, delimiter=',',
                                   quotechar='|', quoting=csv.QUOTE_MINIMAL)
                for ii in range(len(ks)):
                    writer.writerow([tks[ii],(scale*numpy.sqrt(simpsd1d[1][1:-3]))[ii]])
        #bovy_plot.bovy_plot(tks[tks > 0.7],
        #                    scale*numpy.sqrt(simpsd1d[1][1:-3][tks > 0.7]+4./scale**2.*(1.-numpy.tanh(-(tks[tks > 0.7]-0.9)/0.1))/2.),
        #                    'k-.',lw=2.,overplot=True)
        #line2= bovy_plot.bovy_plot(tks,
        #                           scale*numpy.sqrt(simpsd1d[1][1:-3]+4./scale**2.),
        #                           'k-.',lw=2.,overplot=True,dashes=(10,5,3,5))
        l1= pyplot.legend((line1[0],),
#                      (r'$\mathrm{Spiral}:\ \delta \phi_{\mathrm{rms}} = (10\,\mathrm{km\,s}^{-1})^2,$'+'\n'+r'$\mathrm{pitch\ angle} = 8^\circ$'+'\n'+r'$\mathrm{Sun\ near\ 2\!:\!1\ Lindblad\ resonance}$',),
                      (r'$\mathrm{Bar}:\ F_{R,\mathrm{bar}} / F_{R,\mathrm{axi}} = 1.5\%,$'+'\n'+r'$\mathrm{angle} = 25^\circ,$'+'\n'+r'$\mathrm{Sun\ near\ 2\!:\!1\ Lindblad\ resonance}$',),
                      loc='upper right',#bbox_to_anchor=(.91,.375),
                      numpoints=8,
                      prop={'size':14},
                      frameon=False)
    #Add the lopsided and ellipticity constraints from Rix/Zaritsky
    if _ADDRIX:
        pyplot.errorbar([1./16.],[5.6],
                        yerr=[5.6/2.],
                        marker='d',color='0.6',
                        mew=1.5,mfc='none',mec='0.6')
        pyplot.errorbar([1./8./numpy.sqrt(2.)],[6.4],
                        yerr=numpy.reshape(numpy.array([6.4/0.045*0.02,6.4/0.045*0.03]),(2,1)),
                        marker='d',color='0.6',mec='0.6',
                        mew=1.5,mfc='none')
    if _ADDGCS:
        ks_gcs, psd_gcs, e_psd_gcs, gcsp= plot_psd_gcs()
    if _ADDRAVE:
        ks_rave, psd_rave, e_psd_rave, ravep= plot_psd_rave()
    if _ADDRED:
        plot_psd_red()
    l2= pyplot.legend((apop[0],ravep[0],gcsp[0]),
                      (r'$\mathrm{APOGEE}$',
                       r'$\mathrm{RAVE}$',
                       r'$\mathrm{GCS}$'),
                      loc='upper right',bbox_to_anchor=(.95,.750),
                      numpoints=1,
                      prop={'size':14},
                      frameon=False)
    pyplot.gca().add_artist(l1)
    pyplot.gca().add_artist(l2)
    #Plot an estimate of the noise, based on looking at the bands
    nks= numpy.linspace(2.,120.,2)
    if not 'mba23' in socket.gethostname():
        pyplot.fill_between(nks,[2.,2.],hatch='/',color='k',facecolor=(0,0,0,0),
                            lw=0.)
    #                        edgecolor=(0,0,0,0))
    pyplot.plot(nks,[2.,2.],color='k',lw=1.5)
    nks= numpy.linspace(0.17,2.,2)
    def linsp(x):
        return .8/(numpy.log(0.17)-numpy.log(2.))*(numpy.log(x)-numpy.log(0.17))+2.8
    if not 'mba23' in socket.gethostname():
        pyplot.fill_between(nks,linsp(nks),hatch='/',color='k',facecolor=(0,0,0,0),
                            lw=0.)
    #                        edgecolor=(0,0,0,0))   
    #Plot lines
    pyplot.plot([0.17,0.17],[0.,2.8],color='k',lw=1.5)
    pyplot.plot(nks,linsp(nks)+0.01,color='k',lw=1.5)
    bovy_plot.bovy_text(0.19,.5,r'$95\,\%\,\mathrm{noise\ range}$',
                        bbox=dict(facecolor='w',edgecolor='w'),fontsize=14.)
    if _INTERP:
        interpks= list(ks[:-5])
        interppsd= list((scale*numpy.sqrt(psd1d[1][1:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)))[:-5])
        interppsd_w= (scale*0.5*psd1d[2][1:-3]/numpy.sqrt(psd1d[1][1:-3]))[:-5]
        interppsd_w[:8]*= 0.025 #fiddling to get a decent fit
        interppsd_w[3:5]*= 0.001
        interppsd_w= list(interppsd_w)
        interpks.append(0.025)
        interppsd.append(10.**-5.)
        interppsd_w.append(0.00001)
        interpks.append(110.)
        interppsd.append(1.)
        interppsd_w.append(0.00001)
        if _ADDGCS:
            interpks.extend(ks_gcs)
            interppsd.extend(psd_gcs)
            interppsd_w.extend(e_psd_gcs)
        if _ADDRAVE:
            interpks.extend(ks_rave[5:])
            interppsd.extend(psd_rave[5:])
            interppsd_w.extend(e_psd_rave[5:])
        interpks= numpy.array(interpks)
        sortindx= numpy.argsort(interpks)
        interpks= interpks[sortindx]
        interppsd= numpy.array(interppsd)[sortindx]
        interppsd_w= interppsd/numpy.array(interppsd_w)[sortindx]
        interpindx= True-numpy.isnan(interppsd)
    #interpspec= interpolate.InterpolatedUnivariateSpline(interpks[interpindx],
        interpspec= interpolate.UnivariateSpline(numpy.log(interpks[interpindx]),
                                                 numpy.log(interppsd[interpindx]/3.),
                                                 w=interppsd_w,
                                                 k=3,s=len(interppsd_w)*0.8)
        pks= numpy.linspace(interpks[0],interpks[-1],201)
        #bovy_plot.bovy_plot(pks,
        #                    3.*numpy.exp(interpspec(numpy.log(pks))),
#                    'k-',overplot=True)
    def my_formatter(x, pos):
        return r'$%g$' % x
    def my_formatter2(x, pos):
        return r'$%g$' % (1./x)
    major_formatter = FuncFormatter(my_formatter)
    major_formatter2 = FuncFormatter(my_formatter2)
    ax= pyplot.gca()
    ax.xaxis.set_major_formatter(major_formatter)
    if not _PROPOSAL:
        ax2= pyplot.twiny()
        xmin, xmax= ax.xaxis.get_view_interval()
        ax2.set_xscale('log')
        ax2.xaxis.set_view_interval(1./xmin,1./xmax,ignore=True)
        ax2.set_xlabel('$\mathrm{Approximate\ scale}\,(\mathrm{kpc})$',
                       fontsize=12.,ha='center',x=0.5)
        ax2.xaxis.set_major_formatter(major_formatter)
    bovy_plot.bovy_end_print(plotfilename,dpi=300)
    return None

def plot_psd_gcs():
    data= hackGCS.hackGCS()
    dx= _GCSDX
    binsize= .8
    pix= pixelize_sample.pixelXY(data,
                                 xmin=_GCSXMIN,xmax=_GCSXMAX,
                                 ymin=_GCSYMIN,ymax=_GCSYMAX,
                                 dx=dx,dy=dx)
    resv= pix.plot('VVel',returnz=True,justcalc=True)
    resvunc= pix.plot('VVel',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-4))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][1:-3]
    scale= 4.*numpy.pi#3.25/numpy.median(numpy.sqrt(noisepsd))
    print scale
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-4))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][1:-3]
    ks= psd1d[0][1:-3]
    gcsp= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][1:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),'kx',mew=2.,
                        overplot=True)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][1:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*(psd1d[2][1:-3]**2.
                                    +_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)**2.)**0.5/numpy.sqrt(psd1d[1][1:-3]),
                    marker='None',ls='none',color='k')#,lolims=True)
    if _DUMP2FILE:
        with open('bovy-apogee-psd.dat','a') as csvfile:
            writer= csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            csvfile.write('#GCS\n')
            for ii in range(len(ks)):
                writer.writerow([ks[ii],
                                 (scale*numpy.sqrt(psd1d[1][1:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)))[ii],
                                 (scale*0.5*(psd1d[2][1:-3]**2.+_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)**2.)**0.5/numpy.sqrt(psd1d[1][1:-3]))[ii]])
    if _PLOTBAND:
#        bovy_plot.bovy_plot(ks,
#                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
#                            '-',lw=2.,zorder=9,
#                            color='0.65',overplot=True)
#        bovy_plot.bovy_plot(ks,
#                            scale*numpy.sqrt(numpy.sort(noisepsd
#                                                        -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0),axis=0)[int(numpy.floor(0.99*_NNOISE)),:]),
#                            zorder=1,ls='--',overplot=True,
#                            color='0.65')
        pyplot.fill_between(ks,
                            scale*numpy.sqrt(numpy.sort(noisepsd
                                                        -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0),axis=0)[int(numpy.floor(_SIGNIF*_NNOISE)),:]),
                            zorder=0,
                            color='0.65')
    return (ks,
            scale*numpy.sqrt(psd1d[1][1:-3]
                             -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
            scale*0.5*psd1d[2][1:-3]/numpy.sqrt(psd1d[1][1:-3]),
            gcsp)

def plot_psd_rave():
    data= fitsio.read(os.path.join(os.getenv('DATADIR'),'rave','ravedr4_rc.fits'))
    dx= _RAVEDX
    binsize= 0.8#.735
    pix= pix= pixelize_sample.pixelXY(data,
                                      xmin=_RAVEXMIN,xmax=_RAVEXMAX,
                                      ymin=_RAVEYMIN,ymax=_RAVEYMAX,
                                      dx=dx,dy=dx)
    resv= pix.plot(lambda x: dvlosgal(x,vtsun=230.),returnz=True,justcalc=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-4))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][1:-3]
    scale= 4.*numpy.pi#3.25/numpy.median(numpy.sqrt(noisepsd))
    print scale
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-4))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][1:-3]
    ks= psd1d[0][1:-3]
    ravep= bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][1:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),'k+',mew=2.,
                        overplot=True)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][1:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*(psd1d[2][1:-3]**2.
                                    +_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)**2.)**0.5/numpy.sqrt(psd1d[1][1:-3]),
                    marker='None',ls='none',color='k')
    if _DUMP2FILE:
        with open('bovy-apogee-psd.dat','a') as csvfile:
            writer= csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
            csvfile.write('#RAVE\n')
            for ii in range(len(ks)):
                writer.writerow([ks[ii],
                                 (scale*numpy.sqrt(psd1d[1][1:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)))[ii],
                                 (scale*0.5*(psd1d[2][1:-3]**2.+_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)**2.)**0.5/numpy.sqrt(psd1d[1][1:-3]))[ii]])
    if False:
        interpindx= True-numpy.isnan(psd1d[1][1:-3])
        interpspec= interpolate.InterpolatedUnivariateSpline(ks[interpindx],
                                                             scale*numpy.sqrt(psd1d[1][1:-3])[interpindx],
                                                             k=3)
        pks= numpy.linspace(ks[0],ks[-1],201)
        bovy_plot.bovy_plot(pks,interpspec(pks),'k-',overplot=True)
    if _PLOTBAND:
        #bovy_plot.bovy_plot(ks,
        #                    scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
        #                    '-',lw=8.,zorder=9,
        #                    color='0.65',overplot=True)
#        bovy_plot.bovy_plot(ks,
#                            scale*numpy.sqrt(numpy.sort(noisepsd
#                                                        -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0),axis=0)[int(numpy.floor(0.99*_NNOISE)),:]),
#                            zorder=1,ls='--',overplot=True,
#                            color='0.65')
        pyplot.fill_between(ks,
                            scale*numpy.sqrt(numpy.sort(noisepsd
                                                        -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0),axis=0)[int(numpy.floor(_SIGNIF*_NNOISE)),:]),
                            zorder=0,
                            color='0.65')
        #bovy_plot.bovy_plot(numpy.tile(ks,(_NNOISE,1)).T,
        #                    scale*numpy.sqrt(noisepsd).T,
        #                    '-',zorder=0,alpha=0.5,
        #                    color='0.45',overplot=True)
    return (ks,
            scale*numpy.sqrt(psd1d[1][1:-3]
                             -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
            scale*0.5*psd1d[2][1:-3]/numpy.sqrt(psd1d[1][1:-3]),
            ravep)

def plot_psd_red():
    data= readAndHackHoltz.readAndHackHoltz()
    dx= 1.
    binsize= .735
    pix= pixelize_sample.pixelXY(data,xmin=5.,xmax=13,
                                 ymin=-3.,ymax=7.,
                                 dx=dx,dy=dx)
    resv= pix.plot(lambda x: dvlosgal(x),returnz=True,justcalc=True)
    resvunc= pix.plot('VHELIO_AVG',
                      func=lambda x: 1.4826*numpy.median(numpy.fabs(x-numpy.median(x)))/numpy.sqrt(len(x)),
                      returnz=True,justcalc=True)
    psd1d= bovy_psd.psd1d(resv,dx,binsize=binsize)
    print psd1d
    #Simulations for constant 3.25 km/s
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*3.25
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    scale= 4.*numpy.pi#3.25/numpy.median(numpy.sqrt(noisepsd))
    print scale
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    ks= psd1d[0][0:-3]
    bovy_plot.bovy_plot(ks,scale*numpy.sqrt(psd1d[1][0:-3]
                                            -_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),'ko',mew=2.,
                        mfc='none',overplot=True)
    pyplot.errorbar(ks,scale*numpy.sqrt(psd1d[1][0:-3]-_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)),
                    yerr=scale*0.5*(psd1d[2][0:-3]**2.
                                    +_SUBTRACTERRORS*numpy.median(noisepsd,axis=0)**2.)**0.5/numpy.sqrt(psd1d[1][0:-3]),
                    marker='None',ls='none',color='k')
    if False:
        interpindx= True-numpy.isnan(psd1d[1][0:-3])
        interpspec= interpolate.InterpolatedUnivariateSpline(ks[interpindx],
                                                             scale*numpy.sqrt(psd1d[1][0:-3])[interpindx],
                                                             k=3)
        pks= numpy.linspace(ks[0],ks[-1],201)
        bovy_plot.bovy_plot(pks,interpspec(pks),'k-',overplot=True)
    #Simulations for the actual noise
    nnoise= _NNOISE
    noisepsd= numpy.empty((nnoise,len(psd1d[0])-3))
    for ii in range(nnoise):
        newresv= numpy.random.normal(size=resv.shape)*resvunc
        noisepsd[ii,:]= bovy_psd.psd1d(newresv,dx,binsize=binsize)[1][0:-3]
    if _PLOTBAND:
        bovy_plot.bovy_plot(ks,
                            scale*numpy.median(numpy.sqrt(noisepsd),axis=0),
                            '-',lw=8.,zorder=9,
                            color='0.65',overplot=True)
    return None

if __name__ == '__main__':
    plot_psd(sys.argv[1])
