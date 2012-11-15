import sys
import numpy
import fitsio
from galpy.util import bovy_plot
import apogee
import isodist
import rcmodel
def readData(halo=False,bulge=False):
    p= apogee.tools.path.apallPath()
    data= fitsio.read(p)
    if halo:
        #Cuts, halo
        bmin= 16.
        data=data[(numpy.fabs(data['GLAT']) > bmin)]
    elif bulge:
        #Cuts, bulge
        bmax= 8.
        lmin= 357.
        lmax= 22.
        data=data[(numpy.fabs(data['GLAT']) <= bmax)]
        indx1= (data['GLON'] >= lmin)
        indx2= (data['GLON'] <= lmax)
        data= data[indx1+indx2]
    else:
        #Cuts, disk
        bmax= 16.
        lmin= 24.
        lmax= 240.
        data=data[(numpy.fabs(data['GLAT']) <= bmax)*(data['GLON'] >= lmin)\
                      *(data['GLON'] <= lmax)]
    #regular stars
    #indx= numpy.array(['STAR' in data['OBJTYPE'][ii] for ii in range(len(data))],dtype='bool')
    #data= data[indx]
    data= data[(data['PRE'] == 0)] #Pre-shutdown
    data= data[(data['AK_WISE'] != -9999.9999)]
    return data

def plot_apogee_jkz(plotfilename,sample):
    if sample == 'disk':
        data= readData()
    elif sample == 'halo':
        data= readData(halo=True)
    elif sample == 'bulge':
        data= readData(bulge=True)
    ntotal= len(data)
    #logg cut
    data= data[(data['LOGG_AVG'] < 2.8)*(data['LOGG_AVG'] > 1.8)]
    npass= len(data)
    #Deredden
    aj= data['AK_WISE']*2.5
    ah= data['AK_WISE']*1.55
    j0= data['JMAG']-aj
    h0= data['HMAG']-ah
    k0= data['KMAG']-data['AK_WISE']
    #calculate Zs
    zsolar= 0.019
    z= isodist.FEH2Z(data['METALS_AVG'])/zsolar
    #Plot
    bovy_plot.bovy_print()
    if sample == 'bulge':
        bovy_plot.bovy_plot(j0-k0,z,'ko',
                            xrange=[0.5,0.75],
                            yrange=[0.,0.03/zsolar],
                            xlabel=r'$(J-K_s)_0$',
                            ylabel=r'$Z/Z_\odot$',ms=5.)
    else:
        bovy_plot.scatterplot(j0-k0,z,'k,',
                              xrange=[0.5,0.75],
                              yrange=[0.,0.03/zsolar],
                              xlabel=r'$(J-K_s)_0$',
                              ylabel=r'$Z/Z_\odot$',
                              bins=21)
    #Overplot cuts
    jks= numpy.linspace(0.5,0.75,1001)
    bovy_plot.bovy_plot(jks,rcmodel.jkzcut(jks)/zsolar,
                        'k--',lw=2.,overplot=True)
    bovy_plot.bovy_plot(jks,rcmodel.jkzcut(jks,upper=True)/zsolar,
                        'k--',lw=2.,overplot=True)
    #Data between the cuts
    indx= (j0-k0 < 0.75)*(j0-k0 > 0.5)\
                   *(z <= rcmodel.jkzcut(j0-k0,upper=True)/zsolar)\
                   *(z >= rcmodel.jkzcut(j0-k0)/zsolar)
#    j0= j0[indx]
#    h0= h0[indx]
#    k0= k0[indx]
#    z= z[indx]
#    bovy_plot.bovy_plot(j0-k0,z,'w.',
#                        overplot=True,
#                        mec='w',
#                        ms=.5)
    bovy_plot.bovy_text(r'$\mathrm{APOGEE\ ' + sample + '\ sample}$',
                        title=True)
    bovy_plot.bovy_text(r'$%i/%i/%i$' % (numpy.sum(indx),npass,ntotal),
                        bottom_right=True,size=14.)
    bovy_plot.bovy_end_print(plotfilename)
    return None                        

if __name__ == '__main__':
    plot_apogee_jkz(sys.argv[1],sys.argv[2])
