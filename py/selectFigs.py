import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.select.apogeeSelect
_DEBUG= False
_EXT= 'ps'
if __name__ == '__main__':
    #Make various selection function figures
    #Load selection function
    if _DEBUG:
        #Only load a small subset for testing
        apo= apogee.select.apogeeSelect(sample='main',
                                        locations=[4240,4241,4242])
    else:
        apo= apogee.select.apogeeSelect(sample='main')
    #First plot observational progress
    cohorts= ['short','medium','long']
    for cohort in cohorts:
        bovy_plot.bovy_print(fig_width=8.,fig_height=4.)
        apo.plot_obs_progress(cohort=cohort,add_cohort_label=True)
        bovy_plot.bovy_end_print('../tex-catalog/obsprogress-%s.%s' % (cohort,_EXT))
    #Now plot the color-magnitude diagram
    bovy_plot.bovy_print()
    apo.plotColorMag(bins=101,specbins=51,onedhistsbins=201,onedhistsspecbins=101,cntrSmooth=.75)
    bovy_plot.bovy_end_print('../tex-catalog/colorMagSF.%s' % _EXT)
    #Now plot the selection function in real-space
    bovy_plot.bovy_print(fig_width=6.)
    apo.plot_selfunc_xy(vmax=15.)
    bovy_plot.bovy_end_print('../tex-catalog/sfxy.%s' % _EXT)
    bovy_plot.bovy_print(fig_width=6.)
    apo.plot_selfunc_xy(type='rz',vmax=15.)
    bovy_plot.bovy_end_print('../tex-catalog/sfrz.%s' % _EXT)
    #Now calculate the KS probabilities and plot the histogram
    ksshort= apo.check_consistency('short','short')
    ksmedium= apo.check_consistency('medium','medium')
    kslong= apo.check_consistency('long','long')
    xrange= [0.,1.]
    bovy_plot.bovy_print()
    bovy_plot.bovy_hist(ksshort,range=xrange,bins=21,
                        xlabel=r'$\mathrm{KS\ probability}$',#\ that\ the\ spectroscopic\ sample}$'+'\n'+r'$\mathrm{was\ drawn\ from\ the\ photometric\ sample}$'+'\n'+r'$\times\ \mathrm{the\ model\ selection\ function}$',
                        ylabel=r'$\mathrm{distribution}$',
                        yrange=[0.,4.],
                        ec='k',histtype='step',normed=True)
    bovy_plot.bovy_hist(ksmedium,range=xrange,normed=True,
                        bins=11,overplot=True,
                        ec='k',histtype='step',
                        ls='dashed')
    bovy_plot.bovy_hist(kslong,range=xrange,normed=True,
                        bins=11,overplot=True,
                        ec='k',histtype='step',
                        ls='dashdot')
    #Add proxies for legend
    from matplotlib.lines import Line2D
    lineshort= Line2D([0.],[0.],linestyle='-',color='k')
    linemedium= Line2D([0.],[0.],linestyle='--',color='k')
    linelong= Line2D([0.],[0.],linestyle='-.',color='k')
    #Add legend
    pyplot.legend((lineshort,linemedium,linelong),
                      (r'$\mathrm{short\ cohorts}$',
                       r'$\mathrm{medium\ cohorts}$',
                       r'$\mathrm{long\ cohorts}$'),
                  loc='upper left',#bbox_to_anchor=(.91,.375),
                  numpoints=2,
                  prop={'size':16},
                  frameon=False)
    bovy_plot.bovy_end_print('../tex-catalog/sfks.%s' % _EXT)
    #Also plot the cumulative distribution
    bovy_plot.bovy_print()
    lineshort= bovy_plot.bovy_plot(sorted(ksshort),
                                   numpy.linspace(1./len(ksshort),1.,len(ksshort)),
                                   'k-',
                                   xlabel=r'$\mathrm{KS\ probability}$',#\ that\ the\ spectroscopic\ sample}$'+'\n'+r'$\mathrm{was\ drawn\ from\ the\ photometric\ sample}$'+'\n'+r'$\times\ \mathrm{the\ model\ selection\ function}$',
                                   ylabel=r'$\mathrm{cumulative\ distribution}$',
                                   yrange=[0.,1.],
                                   zorder=3)
    linemedium= bovy_plot.bovy_plot(sorted(ksmedium),
                                    numpy.linspace(1./len(ksmedium),1.,len(ksmedium)),
                                    'k--',overplot=True,zorder=2)
    linelong= bovy_plot.bovy_plot(sorted(kslong),
                                  numpy.linspace(1./len(kslong),1.,len(kslong)),
                                  'k-.',overplot=True,zorder=1)
    bovy_plot.bovy_plot([0.,1.],[0.,1.],'-',color='0.5',lw=2.,overplot=True,
                        zorder=0)
    bovy_plot.bovy_end_print('../tex-catalog/sfks_cumul.%s' % _EXT)
