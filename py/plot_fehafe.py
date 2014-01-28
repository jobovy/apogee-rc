import sys
import numpy
from scipy import special
from galpy.util import bovy_plot
from matplotlib import pyplot
import monoAbundanceMW as maps
import apogee.tools.read as apread
_EXT='png'
_ADDLLOGGCUT= True
_SCATTERPLOT= False
_ADDBOVYRIX= True
def plot_fehafe(basefilename):
    #First read the sample
    data= apread.rcsample()
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
                add_bovy_rix(axScatter,Rmins[ii],Rmaxs[ii],zmins[jj],zmaxs[jj])
            #Overplot ridge line
            def ridge1(feh):
                return -0.2*(feh+0.6)+0.23
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
        def ridge1(feh):
            return -0.2*(feh+0.6)+0.23
        bovy_plot.bovy_plot([-0.8,0.3],[ridge1(-0.8),ridge1(0.3)],
                            '-',lw=2.,color='y',overplot=True)
        bovy_plot.bovy_text(r'$%.1f < |z| / \mathrm{kpc} \leq %.1f$' % (zmins[jj],zmaxs[jj]),
                                top_left=True,size=16.)
        bovy_plot.bovy_end_print(basefilename+'_%.1fz%.1f.' % (zmins[jj],zmaxs[jj])+_EXT)
    return None

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

montage -background none fehafe_5R7_1.0z2.0.png fehafe_7R9_1.0z2.0.png fehafe_9R11_1.0z2.0.png fehafe_5R7_0.5z1.0.png fehafe_7R9_0.5z1.0.png fehafe_9R11_0.5z1.0.png fehafe_5R7_0.0z0.5.png fehafe_7R9_0.0z0.5.png fehafe_9R11_0.0z0.5.png -tile 3x3 -geometry +0+0 fehafe.png

montage -background none fehafe_wbovyrix_5R7_1.0z2.0.png fehafe_wbovyrix_7R9_1.0z2.0.png fehafe_wbovyrix_9R11_1.0z2.0.png fehafe_wbovyrix_5R7_0.5z1.0.png fehafe_wbovyrix_7R9_0.5z1.0.png fehafe_wbovyrix_9R11_0.5z1.0.png fehafe_wbovyrix_9R11_0.5z1.0.png fehafe_wbovyrix_7R9_0.0z0.5.png fehafe_wbovyrix_9R11_0.0z0.5.png  -tile 3x3 -geometry +0+0 fehafe_wbovyrix.png

montage -background none fehafe_0.0z0.5.png fehafe_0.5z1.0.png fehafe_1.0z2.0.png fehafe_2.0z4.0.png -tile 2x2 -geometry +0+0 fehafe_zonly.png

montage -background none fehafe_wbovyrix_0.0z0.5.png fehafe_wbovyrix_0.5z1.0.png fehafe_wbovyrix_1.0z2.0.png fehafe_wbovyrix_2.0z4.0.png -tile 2x2 -geometry +0+0 fehafe_wbovyrix_zonly.png



"""
