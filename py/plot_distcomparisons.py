# Plot distance comparisons between APOKASC, RC, and various DR12 distances
import sys
import numpy
from galpy.util import bovy_plot
from matplotlib import pyplot
import apogee.tools.read as apread
def plot_distcomparisons(plotfilename):
    data= apread.allStar(adddist=True)
    plotKASCDiffs= []
    plotKASCDiffErrs= []
    plotRCDiffs= []
    plotRCDiffErrs= []
    # RC vs. APOKASC
    rcDiff= distDiff(data['APOKASC_DIST_DIRECT'],data['RC_DIST'])
    plotKASCDiffs.append(rcDiff[0])
    plotKASCDiffErrs.append(rcDiff[1])
    # RC vs. RC
    plotRCDiffs.append(0.)
    plotRCDiffErrs.append(0.)
    # BPG Dist1 vs. APOKASC
    apokascMinDist= 1.
    bpgDiff= distDiff(data['APOKASC_DIST_DIRECT'],data['BPG_DIST1_MEAN'],
                      minDist=apokascMinDist)
    plotKASCDiffs.append(bpgDiff[0])
    plotKASCDiffErrs.append(bpgDiff[1])
    # BPG vs. RC
    bpgDiff= distDiff(data['RC_DIST'],data['BPG_DIST1_MEAN'])
    plotRCDiffs.append(bpgDiff[0])
    plotRCDiffErrs.append(bpgDiff[1])
    # Hayden peak vs. APOKASC
    haydenDiff= distDiff(data['APOKASC_DIST_DIRECT'],data['HAYDEN_DIST_PEAK'],
                         minDist=apokascMinDist)
    plotKASCDiffs.append(haydenDiff[0])
    plotKASCDiffErrs.append(haydenDiff[1])
    # Hayden peak vs. RC
    haydenDiff= distDiff(data['RC_DIST'],data['HAYDEN_DIST_PEAK'])
    plotRCDiffs.append(haydenDiff[0])
    plotRCDiffErrs.append(haydenDiff[1])
    # Schultheis vs. APOKASC
    schultheisDiff= distDiff(data['APOKASC_DIST_DIRECT'],data['SCHULTHEIS_DIST'],
                             minDist=apokascMinDist)
    plotKASCDiffs.append(schultheisDiff[0])
    plotKASCDiffErrs.append(schultheisDiff[1])
    # Schultheis vs. RC
    schultheisDiff= distDiff(data['RC_DIST'],data['SCHULTHEIS_DIST'])
    plotRCDiffs.append(schultheisDiff[0])
    plotRCDiffErrs.append(schultheisDiff[1])
    # plot
    bovy_plot.bovy_print(fig_width=7.,
                         text_fontsize=20.,
                         legend_fontsize=24.,
                         xtick_labelsize=18.,
                         ytick_labelsize=18.,
                         axes_labelsize=24.)
    ms= 8.
    line1= bovy_plot.bovy_plot([1,2,3,4],plotKASCDiffs,'bo',ms=ms,
                        xrange=[0,5],
                        yrange=[-0.11,0.11],
                        ylabel=r'$\mathrm{distance\ modulus\ offset}$')
    pyplot.errorbar([1,2,3,4],plotKASCDiffs,yerr=plotKASCDiffErrs,
                    ls='none',marker='o',color='b',ms=ms)
    line2= bovy_plot.bovy_plot([1,2,3,4],plotRCDiffs,'ro',ms=ms,
                        overplot=True)
    pyplot.errorbar([1,2,3,4],plotRCDiffs,yerr=plotRCDiffErrs,
                    ls='none',marker='o',color='r',ms=ms)
    pyplot.legend([line1[0],line2[0]],
                  [r'$\mathrm{wrt\ APOKASC}$',
                   r'$\mathrm{wrt\ RC}$'],
                  loc='lower left',fontsize=14.,frameon=False,numpoints=1)
    #Put labels and rotate them
    pyplot.xticks([1,2,3,4],
                  [r'$\mathrm{RC}$',
                   r"$\mathrm{BPG\ dist1}$",
                   r"$\mathrm{Hayden\ peak}$",
                   r"$\mathrm{Schultheis}$"],size=16.,
                  rotation=45.)
    bovy_plot.bovy_end_print(plotfilename)   
    return None

def distDiff(ref,comp,minDist=1.):
    """Return the median distance difference and its error for distances > 1 kpc in ref"""
    compIndx= (ref > minDist)*(comp > -0.0001)
    medDiff= 5.*(numpy.log10(comp)-numpy.log10(ref))[compIndx]
    return (numpy.median(medDiff),
            1.4826*numpy.median(numpy.fabs(medDiff-numpy.median(medDiff)))/numpy.sqrt(float(numpy.sum(compIndx))))

if __name__ == '__main__':
    plot_distcomparisons(sys.argv[1])
