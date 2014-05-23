import os
import sys
import csv
import numpy
import esutil
from galpy.util import bovy_plot
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
import apogee.tools.read as apread
import rcmodel
from rcmodel import rcdist
import asciitable
def read_saga():
    """Read the SAGA data"""
    filename= os.path.join('..','data','star_prop_uvby.dat')
    readme= os.path.join('..','data','ReadMe.txt')
    sagadata= asciitable.read(filename,
                              readme=readme,
                              Reader=asciitable.cds.Cds,
                              guess=False,
                              fill_values=[('', '-999')])
    return sagadata

def match_apokasc_saga():
    sagadata= read_saga()
    ids= sagadata['KICID'].data
    apokasc= apread.apokasc()
    dists= numpy.zeros(len(apokasc))-1
    edists= numpy.zeros(len(apokasc))-1
    for ii in range(len(apokasc)):
        indx= ids == apokasc['KEPLER ID'][ii]
        if numpy.sum(indx) == 0: continue
        dists[ii]= sagadata['Dis'][indx]
        edists[ii]= 0.5*(sagadata['E_Dis'][indx]+sagadata['e_Dis'][indx])
    apokasc= esutil.numpy_util.add_fields(apokasc,[('DIST_SEISMO', float),
                                                   ('E_DIST_SEISMO', float),
                                                   ])
    apokasc['DIST_SEISMO']= dists/1000.
    apokasc['E_DIST_SEISMO']= edists/1000.
    return apokasc

def compare_seismic_distances(plotfilename):
    apokasc= match_apokasc_saga()
    #Perform RC selection
    logg= apokasc['KASC_RG_LOGG_SCALE_2']
    teff= apokasc['TEFF']
    z= 0.017*10.**apokasc['METALS']
    jk= apokasc['J0']-apokasc['K0']
    indx= (logg >= 1.8)\
        *(logg <= 0.0018*(teff+382.5*apokasc['METALS']-4607)+2.5)\
        *(jk < 0.8)\
        *(jk >= 0.5)\
        *(z <= 0.06)\
        *(z <= rcmodel.jkzcut(jk,upper=True))\
        *(z >= rcmodel.jkzcut(jk))#\
        #*(apokasc['SEISMO EVOL'] == 'CLUMP')
    print "Found %i RC stars in APOKASC" % numpy.sum(indx)
    rcdists= numpy.zeros(len(apokasc))-1
    rcd= rcdist('../data/rcmodel_mode_jkz_ks_parsec_newlogg.sav')
    rcdists[indx]= rcd(jk[indx],z[indx],apokasc['K0'][indx])
    pindx= (rcdists > 0.)*(apokasc['DIST_SEISMO'] > 0.)
    print "Found %i RC stars in APOKASC with seismic distances" % numpy.sum(pindx)
    #Setup plot
    bovy_plot.bovy_print(fig_height=7.)
    dx= 0.6
    left, bottom, width, height= 0.1, 0.9-dx, 0.8, dx
    axTop= pyplot.axes([left,bottom,width,height])
    fig= pyplot.gcf()
    fig.sca(axTop)
    bovy_plot.bovy_plot([0.,20.],[0.,20.],'k-',lw=2.,color='0.4',
                        overplot=True,zorder=0)
    bovy_plot.bovy_plot(rcdists[pindx],apokasc['DIST_SEISMO'][pindx],
                        'k.',overplot=True,zorder=10)
    if False:
        pyplot.errorbar(rcdists[pindx],
                        apokasc['DIST_SEISMO'][pindx],
                        xerr=0.05*rcdists[pindx],
                        yerr=apokasc['E_DIST_SEISMO'][pindx],
                        marker=',',color='k',
                        linestyle='none')
    thisax= pyplot.gca()
    thisax.set_ylim(0.,5.)
    pyplot.xlim(0.,5.)
    bovy_plot._add_ticks()
    nullfmt   = NullFormatter()         # no labels
    axTop.xaxis.set_major_formatter(nullfmt)
    bovy_plot._add_ticks()
    pyplot.ylabel(r'$\mathrm{seismic\ distance}\,(\mathrm{kpc})$')
    #Second plot
    left, bottom, width, height= 0.1, 0.1, 0.8, 0.8-dx
    thisax= pyplot.axes([left,bottom,width,height])
    fig.sca(thisax)
    bovy_plot.bovy_plot([0.,20.],[0.,0.],'k-',lw=2.,color='0.4',
                        overplot=True,zorder=0)
    bovy_plot.bovy_plot(rcdists[pindx],
                        (apokasc['DIST_SEISMO'][pindx]-rcdists[pindx])/rcdists[pindx],
                        'k.',overplot=True,zorder=10)
    thisax= pyplot.gca()
    thisax.set_ylim(-0.2,0.2)
    pyplot.xlim(0.,5.)
    bovy_plot._add_ticks()
    nullfmt   = NullFormatter()         # no labels
    bovy_plot._add_ticks()
    pyplot.ylabel(r'$\mathrm{relative\ differene}$')
    pyplot.xlabel(r'$\mathrm{RC\ distance}\,(\mathrm{kpc})$')
    medoffset= numpy.median((apokasc['DIST_SEISMO'][pindx]-rcdists[pindx])/rcdists[pindx])
    medsig= 1.4826*numpy.median(numpy.fabs((apokasc['DIST_SEISMO'][pindx]-rcdists[pindx])/rcdists[pindx]-medoffset))
    bovy_plot.bovy_text(2.75,-0.125,r'$\mathrm{diff} = %.3f\pm%.3f$' % \
                            (medoffset,medsig),size=14.)
    bovy_plot.bovy_end_print(plotfilename)

def write_matches(savefilename):
    apokasc= match_apokasc_saga()
    #Perform RC selection
    logg= apokasc['KASC_RG_LOGG_SCALE_2']
    teff= apokasc['TEFF']
    z= 0.017*10.**apokasc['METALS']
    jk= apokasc['J0']-apokasc['K0']
    indx= (logg >= 1.8)\
        *(logg <= 0.0018*(teff+382.5*apokasc['METALS']-4607)+2.5)\
        *(jk < 0.8)\
        *(jk >= 0.5)\
        *(z <= 0.06)\
        *(z <= rcmodel.jkzcut(jk,upper=True))\
        *(z >= rcmodel.jkzcut(jk))
    print "Found %i RC stars in APOKASC" % numpy.sum(indx)
    rcdists= numpy.zeros(len(apokasc))-1
    rcd= rcdist('../data/rcmodel_mode_jkz_ks_parsec_newlogg.sav')
    rcdists[indx]= rcd(jk[indx],z[indx],apokasc['K0'][indx])
    pindx= (rcdists > 0.)*(apokasc['DIST_SEISMO'] > 0.)
    apokasc= apokasc[pindx]
    rcdists= rcdists[pindx]
    print "Found %i RC stars in APOKASC with seismic distances" % numpy.sum(pindx)
    savefile= open(savefilename,'w')
    savefile.write('#ID,RCDIST(pc),SEISMODIST(SCALE,MO,pc),(SEISMODIST-RCDIST)/RCDIST\n')
    csvwriter = csv.writer(savefile, delimiter=',')
    for ii in range(len(apokasc)):
        csvwriter.writerow([apokasc['KEPLER ID'][ii],
                            rcdists[ii]*1000.,
                            apokasc['DIST_SEISMO'][ii]*1000.,
                            (apokasc['DIST_SEISMO'][ii]-rcdists[ii])/rcdists[ii]])
    savefile.close()

if __name__ == '__main__':
    if sys.argv[2] == 'match':
        write_matches(sys.argv[1])
    else:
        compare_seismic_distances(sys.argv[1])
        
