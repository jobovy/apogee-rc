import numpy
import esutil
from galpy.util import bovy_plot
import apogee.tools.read as apread
def match_apokasc_rc():
    #First read apokasc
    kascdata= apread.apokasc()
    #Then read rc
    rcdata= apread.rcsample()
    #Match
    h=esutil.htm.HTM()
    m1,m2,d12 = h.match(kascdata['RA'],kascdata['DEC'],
                        rcdata['RA'],rcdata['DEC'],
                        2./3600.,maxmatch=1)
    kascdata= esutil.numpy_util.add_fields(kascdata,[('RC', int)])
    kascdata['RC']= 0
    kascdata['RC'][m1]= 1
    return kascdata

if __name__ == '__main__':
    data= match_apokasc_rc()
    clumpseismo= data['SEISMO EVOL'] == 'CLUMP'
    noseismo= data['SEISMO EVOL'] == 'UNKNOWN'
    noclumpseismo= (data['SEISMO EVOL'] == 'RGB') \
        + (data['SEISMO EVOL'] == 'DWARF/SUBGIANT')
    rcclumpseismo= clumpseismo*(data['RC'] == 1)
    rcnoclumpseismo= noclumpseismo*(data['RC'] == 1)
    #Statistics using evolutionary state measurements
    print "%i APOKASC stars have evolutionary state measurements" % (numpy.sum(clumpseismo)+numpy.sum(noclumpseismo))
    print "%i APOKASC RC stars have evolutionary state measurements" % (numpy.sum(rcclumpseismo)+numpy.sum(rcnoclumpseismo))
    print "%i / %i = %i%% APOKASC CLUMP stars are in the RC catalog" % (numpy.sum(rcclumpseismo),numpy.sum(clumpseismo),float(numpy.sum(rcclumpseismo))/numpy.sum(clumpseismo)*100.)
    print "%i / %i = %i%% APOKASC non-CLUMP stars are in the RC catalog" % (numpy.sum(rcnoclumpseismo),numpy.sum(noclumpseismo),float(numpy.sum(rcnoclumpseismo))/numpy.sum(noclumpseismo)*100.)
    print "%i / %i = %i%% APOKASC non-CLUMP stars out of all stars with evolutionary measurements are in the RC catalog" % (numpy.sum(rcnoclumpseismo),numpy.sum(rcnoclumpseismo)+numpy.sum(rcclumpseismo),float(numpy.sum(rcnoclumpseismo))/(numpy.sum(rcnoclumpseismo)+numpy.sum(rcclumpseismo))*100.)
    bovy_plot.bovy_print()
    rcindx= data['RC'] == 1
    bovy_plot.bovy_plot(data['METALS'][rcnoclumpseismo],
                        data['ALPHAFE'][rcnoclumpseismo],
                        'ro',
                        xrange=[-1.,0.7],
                        yrange=[-0.15,0.35],
                        xlabel=r'$[\mathrm{Fe/H}]$',
                        ylabel=r'$[\alpha/\mathrm{Fe}]$',
                        zorder=1,ms=4.5,mec='none',
                        onedhists=True,onedhistxnormed=True,
                        onedhistynormed=True,onedhistec='r',bins=20)

    bovy_plot.bovy_plot(data['METALS'][rcindx],data['ALPHAFE'][rcindx],
                        'k.',overplot=True,zorder=0,
                        onedhists=True,onedhistxnormed=True,
                        onedhistynormed=True,onedhistec='k',bins=20)
    bovy_plot.bovy_end_print('apokasc_rc_metalsafe.png')
    #Statistics using simple seismo logg cut
    clumplogg= (data['KASC_RG_LOGG_SCALE_2'] > 1.8)\
        *(data['KASC_RG_LOGG_SCALE_2'] < 2.8)
    rcclumplogg= clumplogg*(data['RC'] == 1)
    print "%i / %i = %i%% APOKASC logg clump stars are in the RC catalog" % (numpy.sum(rcclumplogg),numpy.sum(clumplogg),float(numpy.sum(rcclumplogg))/numpy.sum(clumplogg)*100)
    rcnoclumplogg= (True-clumplogg)*(data['RC'] == 1)
    print "%i / %i = %i%% APOKASC logg non-clump stars are in the RC catalog" % (numpy.sum(rcnoclumplogg),numpy.sum(True-clumplogg),float(numpy.sum(rcnoclumplogg))/numpy.sum(True-clumplogg)*100.)
    print "%i / %i = %i%% APOKASC logg non-clump stars out of all stars are in the RC catalog" % (numpy.sum(rcnoclumplogg),numpy.sum(data['RC'] == 1),float(numpy.sum(rcnoclumplogg))/numpy.sum(data['RC'] == 1)*100.)
    bovy_plot.bovy_print()
    bovy_plot.bovy_plot(data['METALS'][rcnoclumplogg],
                        data['ALPHAFE'][rcnoclumplogg],
                        'ro',ms=4.5,
                        xrange=[-1.,0.7],
                        yrange=[-0.15,0.35],
                        xlabel=r'$[\mathrm{Fe/H}]$',
                        ylabel=r'$[\alpha/\mathrm{Fe}]$',onedhists=True,
                        mec='none',
                        onedhistxnormed=True,onedhistynormed=True,
                        onedhistcolor='r',zorder=1,onedhistec='r',
                        bins=30)
    bovy_plot.bovy_plot(data['METALS'][rcindx],data['ALPHAFE'][rcindx],
                        'k.',overplot=True,onedhists=True,
                        onedhistxnormed=True,onedhistynormed=True,
                        onedhistcolor='k',zorder=0,bins=30)
    bovy_plot.bovy_end_print('apokasc_rclogg_metalsafe.png')

