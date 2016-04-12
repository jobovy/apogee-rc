import sys
import numpy
from scipy import special
from galpy.util import bovy_plot
from matplotlib import pyplot
from extreme_deconvolution import extreme_deconvolution
import monoAbundanceMW as maps
import apogee.tools.read as apread
_EXT='png'
_ADDLLOGGCUT= False
_SCATTERPLOT= False
_ADDBOVYRIX= False
_NOLINE= False
def plot_fehafe(basefilename):
    #First read the sample
    data= apread.rcsample()
    data= data[data['SNR'] >= 70.]
    if _ADDLLOGGCUT:
        data= data[data['ADDL_LOGG_CUT'] == 1]
    #R and z ranges to plot
    Rmins= [5.,7.,9.]
    Rmaxs= [7.,9.,11.]
    zmins= [0.,0.5,1.]
    zmaxs= [0.5,1.,2.]
    nRbins= len(Rmins)
    nzbins= len(zmins)
    for ii in range(nRbins):
        for jj in range(nzbins):
            #Cut data to relevant range
            tindx= (data['RC_GALR'] > Rmins[ii])\
                *(data['RC_GALR'] <= Rmaxs[ii])\
                *(numpy.fabs(data['RC_GALZ']) > zmins[jj])\
                *(numpy.fabs(data['RC_GALZ']) <= zmaxs[jj])
            tdata= data[tindx]
            bovy_plot.bovy_print()
            if (_SCATTERPLOT or len(tdata) > 1500) and not _ADDBOVYRIX:
                bovy_plot.scatterplot(tdata['METALS'],
                                      tdata['ALPHAFE'],
                                      'k.',
                                      levels=[0.68,.90],
                                      #special.erf(numpy.arange(1,3)/numpy.sqrt(2.)),
                                      xrange=[-1.2,0.7],
                                      yrange=[-0.12,0.42],
                                      xlabel=r'$[\mathrm{Fe/H}]$',
                                      ylabel=r'$[\alpha/\mathrm{H}]$',
                                      onedhists=True,bins=31,onedhistsbins=21)
            else:
                axScatter,axHistx, axHisty=\
                    bovy_plot.bovy_plot(tdata['METALS'],
                                        tdata['ALPHAFE'],
                                        'k.',
                                        xrange=[-1.2,0.7],
                                        yrange=[-0.12,0.42],
                                        xlabel=r'$[\mathrm{Fe/H}]$',
                                        ylabel=r'$[\alpha/\mathrm{H}]$',
                                        onedhists=True,bins=21)
            if _ADDBOVYRIX:
                add_bovy_rix(axScatter,Rmins[ii],Rmaxs[ii],zmins[jj],zmaxs[jj])
            if not _NOLINE:
                #Overplot ridge line
                bovy_plot.bovy_plot([-0.8,0.3],[ridge1(-0.8),ridge1(0.3)],
                                    '-',lw=2.,color='y',overplot=True)
            bovy_plot.bovy_text(r'$%i < R / \mathrm{kpc} \leq %i$' % (Rmins[ii],Rmaxs[ii]) +'\n'+r'$%.1f < |z| / \mathrm{kpc} \leq %.1f$' % (zmins[jj],zmaxs[jj]),
                                    top_left=True,size=16.)
            bovy_plot.bovy_end_print(basefilename+'_%iR%i_%.1fz%.1f.' % (Rmins[ii],Rmaxs[ii],zmins[jj],zmaxs[jj])+_EXT)
    #Plot z bins only
    zmins= [0.,0.5,1.,2.]
    zmaxs= [0.5,1.,2.,4.]
    nzbins= len(zmins)
    for jj in range(nzbins):
        #Cut data to relevant range
        tindx= (data['RC_GALR'] > 4.)\
            *(data['RC_GALR'] <= 9.)\
                *(numpy.fabs(data['RC_GALZ']) > zmins[jj])\
                *(numpy.fabs(data['RC_GALZ']) <= zmaxs[jj])
        tdata= data[tindx]
        bovy_plot.bovy_print()
        if (_SCATTERPLOT or len(tdata) > 1500) and not _ADDBOVYRIX:
            bovy_plot.scatterplot(tdata['METALS'],
                                  tdata['ALPHAFE'],
                                  'k.',
                                  levels=[0.68,.90],
                                  xrange=[-1.2,0.7],
                                  yrange=[-0.12,0.42],
                                  xlabel=r'$[\mathrm{Fe/H}]$',
                                  ylabel=r'$[\alpha/\mathrm{H}]$',
                                  onedhists=True,bins=31,onedhistsbins=31)
        else:
            axScatter,axHistx, axHisty=\
                bovy_plot.bovy_plot(tdata['METALS'],
                                    tdata['ALPHAFE'],
                                    'k.',
                                    xrange=[-1.2,0.7],
                                    yrange=[-0.12,0.42],
                                    xlabel=r'$[\mathrm{Fe/H}]$',
                                    ylabel=r'$[\alpha/\mathrm{H}]$',
                                    onedhists=True,bins=21)
        if _ADDBOVYRIX:
            add_bovy_rix(axScatter,4.,9.,zmins[jj],zmaxs[jj])
        #Overplot ridge line
        bovy_plot.bovy_plot([-0.8,0.3],[ridge1(-0.8),ridge1(0.3)],
                            '-',lw=2.,color='y',overplot=True)
        bovy_plot.bovy_text(r'$%.1f < |z| / \mathrm{kpc} \leq %.1f$' % (zmins[jj],zmaxs[jj]),
                                top_left=True,size=16.)
        bovy_plot.bovy_end_print(basefilename+'_%.1fz%.1f.' % (zmins[jj],zmaxs[jj])+_EXT)
    #Plot scatter in the high-alpha sequence
    #Histograms in 3 R bins
    bovy_plot.bovy_print(fig_width=7.)
    overplot= False
    colors= ['y','k','r']
    labelposx= 0.235
    labelposy= [13.,11.,9.]
    for ii in range(nRbins):
        tindx= (data['RC_GALR'] > Rmins[ii])\
            *(data['RC_GALR'] <= Rmaxs[ii])\
            *(numpy.fabs(data['RC_GALZ']) <= 3.)
        tdata= data[tindx]
        hindx= determine_highalpha(tdata)
        tdata= tdata[hindx]
        bovy_plot.bovy_hist(tdata['ALPHAFE']-ridge1(tdata['METALS']),
                            color=colors[ii],histtype='step',
                            bins=28,range=[-0.275,0.25],yrange=[0.,15.],
                            weights=numpy.ones(len(tdata)),
                            xlabel=r'$\delta[\alpha\mathrm{/Fe}]$',
                            overplot=overplot,normed=True)
        dafe= tdata['ALPHAFE']-ridge1(tdata['METALS'])
        dafe= dafe[(numpy.fabs(tdata['RC_GALZ']) > .5)]
        txa, txm, txc= run_xd(dafe)
        print txm.flatten(), numpy.sqrt(txc.flatten()), txa[0]*len(dafe)
        bovy_plot.bovy_plot([txm[0,0],txm[0,0]],[0.,200.],'--',lw=2.,
                            color=colors[ii],overplot=True)
        bovy_plot.bovy_text(labelposx,labelposy[ii],
                            r'$%i < R \leq %i: \sigma = %.3f$' % (Rmins[ii],Rmaxs[ii],numpy.sqrt(txc[0,0,0])),
                            size=16.,color=colors[ii],ha='right')
        overplot= True
    bovy_plot.bovy_end_print(basefilename+'_afefehhist.%s' % _EXT)
    #High SN
    bovy_plot.bovy_print(fig_width=7.)
    overplot= False
    for ii in range(nRbins):
        tindx= (data['RC_GALR'] > Rmins[ii])\
            *(data['RC_GALR'] <= Rmaxs[ii])\
            *(numpy.fabs(data['RC_GALZ']) <= 3.)\
            *(data['SNR'] >= 150.)
        tdata= data[tindx]
        hindx= determine_highalpha(tdata)
        tdata= tdata[hindx]
        bovy_plot.bovy_hist(tdata['ALPHAFE']-ridge1(tdata['METALS']),
                            color=colors[ii],histtype='step',
                            bins=19,range=[-0.275,0.15],yrange=[0.,19.5],
                            weights=numpy.ones(len(tdata)),
                            xlabel=r'$\delta[\alpha\mathrm{/Fe}]$',
                            normed=True,
                            overplot=overplot)
        dafe= tdata['ALPHAFE']-ridge1(tdata['METALS'])
        dafe= dafe[(numpy.fabs(tdata['RC_GALZ']) > .5)]
        txa, txm, txc= run_xd(dafe)
        print txm.flatten(), numpy.sqrt(txc.flatten())
        bovy_plot.bovy_plot([txm[0,0],txm[0,0]],[0.,200.],'--',lw=2.,
                            color=colors[ii],overplot=True)
        overplot= True
    bovy_plot.bovy_end_print(basefilename+'_afefehhist_highsn.%s' % _EXT)
    return None

def ridge1(feh):
    return -0.2*(feh+0.6)+0.22

def determine_highalpha(data):
    """Determine which stars are in the high alpha sequence"""
    tafe= data['ALPHAFE']-ridge1(data['METALS'])
    niter= 3
    out= (data['METALS'] <= -0.2)*(data['METALS'] >= -0.6)
    return out
    nsig= [1.5,2.,3.]
    for ii in range(niter):
        tdisp= numpy.std(tafe[out])
        out= (numpy.fabs(tafe) <= nsig[ii]*tdisp)
    print numpy.mean(tafe[out]), tdisp, numpy.sum(out)
    return out

def run_xd(dafe):
    """Run XD on the delta afes"""
    ydata= numpy.empty((len(dafe),1))
    ycovar= numpy.zeros((len(dafe),1))
    ydata[:,0]= dafe
    xamp= numpy.array([0.5,0.5])
    xmean= numpy.empty((2,1))
    xcovar= numpy.empty((2,1,1))
    xmean[0,0]= 0.
    xmean[1,0]= -0.12
    xcovar[0,0,0]= 0.07
    xcovar[1,0,0]= 0.07
    extreme_deconvolution(ydata,ycovar,xamp,xmean,xcovar,fixmean=[False,True])
    return (xamp,xmean,xcovar)

def add_bovy_rix(axScatter,Rmin,Rmax,zmin,zmax):
    #Now calculate and plot our model
    fehs= numpy.linspace(-1.6,0.5,26)
    afes= numpy.linspace(-0.15,0.55,25)
    ourDist= numpy.zeros((len(fehs),len(afes)))
    for ii in range(len(fehs)):
        for jj in range(len(afes)):
            ourDist[ii,jj]= maps.abundanceDist(fehs[ii],afes[jj],
                                               z=[zmin*1000.,zmax*1000.],
                                               r=[Rmin,Rmax])
    _ADJUSTABUNDANCES= True
    _ADJUSTFEH= True
    if _ADJUSTABUNDANCES:
        if _ADJUSTFEH:
            fehs+= 0.15
        afes/= 1.5
        afes-= 0.025
    #Contour this in axScatter
    ourDist[numpy.isnan(ourDist)]= 0.
    sortindx= numpy.argsort(ourDist.flatten())[::-1]
    cumul= numpy.cumsum(numpy.sort(ourDist.flatten())[::-1])/numpy.sum(ourDist.flatten())
    cntrThis= numpy.zeros(numpy.prod(ourDist.shape))
    cntrThis[sortindx]= cumul
    cntrThis= numpy.reshape(cntrThis,ourDist.shape)
    pyplot.sca(axScatter)
    CS= pyplot.contour(fehs,afes,
                       cntrThis.T,levels=[0.68,0.95],
                       linewidths=3.,linestyles='dashed',
                       colors='r',zorder=10)

if __name__ == '__main__':
    plot_fehafe(sys.argv[1])

"""
python plot_fehafe.py ../figs/fehafe
python plot_fehafe.py ../figs/fehafe_wbovyrix
(first set _ADDBOVYRIX)
python plot_fehafe.py ../figs/fehafe_noline
(first set _NOLINE)

montage -background none fehafe_5R7_1.0z2.0.png fehafe_7R9_1.0z2.0.png fehafe_9R11_1.0z2.0.png fehafe_5R7_0.5z1.0.png fehafe_7R9_0.5z1.0.png fehafe_9R11_0.5z1.0.png fehafe_5R7_0.0z0.5.png fehafe_7R9_0.0z0.5.png fehafe_9R11_0.0z0.5.png -tile 3x3 -geometry +0+0 fehafe.png

montage -background none fehafe_wbovyrix_5R7_1.0z2.0.png fehafe_wbovyrix_7R9_1.0z2.0.png fehafe_wbovyrix_9R11_1.0z2.0.png fehafe_wbovyrix_5R7_0.5z1.0.png fehafe_wbovyrix_7R9_0.5z1.0.png fehafe_wbovyrix_9R11_0.5z1.0.png fehafe_wbovyrix_5R7_0.0z0.5.png fehafe_wbovyrix_7R9_0.0z0.5.png fehafe_wbovyrix_9R11_0.0z0.5.png  -tile 3x3 -geometry +0+0 fehafe_wbovyrix.png

montage -background none fehafe_noline_5R7_1.0z2.0.png fehafe_noline_7R9_1.0z2.0.png fehafe_noline_9R11_1.0z2.0.png fehafe_noline_5R7_0.5z1.0.png fehafe_noline_7R9_0.5z1.0.png fehafe_noline_9R11_0.5z1.0.png fehafe_noline_5R7_0.0z0.5.png fehafe_noline_7R9_0.0z0.5.png fehafe_noline_9R11_0.0z0.5.png  -tile 3x3 -geometry +0+0 fehafe_noline.png

montage -background none fehafe_0.0z0.5.png fehafe_0.5z1.0.png fehafe_1.0z2.0.png fehafe_2.0z4.0.png -tile 2x2 -geometry +0+0 fehafe_zonly.png

montage -background none fehafe_wbovyrix_0.0z0.5.png fehafe_wbovyrix_0.5z1.0.png fehafe_wbovyrix_1.0z2.0.png fehafe_wbovyrix_2.0z4.0.png -tile 2x2 -geometry +0+0 fehafe_wbovyrix_zonly.png



"""
