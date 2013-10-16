import math
import numpy
from galpy.df import dehnendf
from galpy.util import save_pickles
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
    ds= numpy.linspace(0.,10./8.,nds)
    dfc= dehnendf(profileParams=(3./8.,2.,36./218.),beta=0.,correct=False)
    vlos= numpy.zeros(nds)
    for jj in range(nds):
        print jj
        d,l= ds[jj], 90./180.*numpy.pi
        R,theta,d,l= safe_dl_to_rphi(d,l)
        vlos[jj]= dfc.meanvT(R)*math.sin(theta+l)
    save_pickles('vhelio_l90.sav',ds*8.,vlos*218.-218.)
