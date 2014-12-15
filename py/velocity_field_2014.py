import sys
import os, os.path
import shutil
import cPickle as pickle
import math as m
import numpy as nu
from optparse import OptionParser
import subprocess
import multiprocessing
from galpy.util import save_pickles, bovy_plot, bovy_coords, multi
from galpy.orbit import Orbit
from galpy.df import dehnendf, shudf
from galpy.potential import SteadyLogSpiralPotential, \
    LogarithmicHaloPotential, TransientLogSpiralPotential, \
    DehnenBarPotential, PowerSphericalPotential, \
    MovingObjectPotential, lindbladR, EllipticalDiskPotential
from galpy.df_src.evolveddiskdf import evolveddiskdf
from matplotlib import pyplot
from matplotlib import transforms
from plot_psd import _RCXMIN, _RCXMAX, _RCYMIN, _RCYMAX, _RCDX
def velocity_field(parser):
    """
    NAME:
       velocity_field
    PURPOSE:
       calculate the velocity field
    INPUT:
       parser - from optparse
    OUTPUT:
       stuff as specified by the options
    HISTORY:
       2011-04-02 - Written - Bovy (NYU)
    """
    (options,args)= parser.parse_args()
    if len(args) == 0:
        parser.print_help()
        sys.exit(-1)
    #Set up grid
    if options.galcoords:
        if options.res == 1: #Just do phi
            xgrid= nu.linspace(0.,m.pi, #ONLY ONE HALF OF TWOPI
                               options.phires)
            ygrid= nu.array([1.])
        else:
            if not options.phimax is None:
                xgrid= nu.linspace(-options.phimax,
                                    options.phimax,
                                    options.phires)
            else:
                xgrid= nu.linspace(0.,2.*m.pi*(1.-1./options.res/2.),
                                   2.*options.res)
            ygrid= nu.linspace(options.rmin,options.rmax,options.res)
    elif options.altrect:
        #Grid we do the RC analysis on
        xgrid= nu.linspace((_RCXMIN-2.25-8.)/8.+_RCDX/8./2.,
                           (_RCXMAX+2.25-8.)/8.-_RCDX/8./2.,
                           options.res)
        ygrid= nu.linspace(_RCYMIN/8.-2.25/8.+_RCDX/8./2.,
                           _RCYMAX/8.+2.25/8.-_RCDX/8./2.,
                           options.res)
    else:
        #Grid we do the RC analysis on
        xgrid= nu.linspace((_RCXMIN-8.)/8.+_RCDX/8./2.,
                           (_RCXMAX-8.)/8.-_RCDX/8./2.,
                           options.res)
        ygrid= nu.linspace(_RCYMIN/8.+_RCDX/8./2.,
                           _RCYMAX/8.-_RCDX/8./2.,
                           options.res)
    nx= len(xgrid)
    ny= len(ygrid)
    #Set up potential
    if options.beta == 0.:
        axip= LogarithmicHaloPotential(normalize=1.,q=0.95)
    else:
        axip= PowerSphericalPotential(alpha=2.-2.*options.beta,normalize=1.)
    pot= [axip]
    if options.bar:
        barp= DehnenBarPotential(alpha=options.baralpha,
                                 tform=options.bar_tform,
                                 tsteady=options.bar_tsteady,
                                 beta=options.beta,
                                 rolr=options.bar_olr,
                                 barphi=options.bar_angle/180.*nu.pi)
        pot.append(barp)
    if options.steadyspiral:
        omegas= options.steadyspiralomegas
        ts= 2.*m.pi/omegas #4:1 ILR at Solar circle
        steadyspiralp= SteadyLogSpiralPotential(A=-0.075,
                                                alpha=options.steadyspiralalpha,
                                                tform=options.steadyspiraltform/ts,
                                                tsteady=options.steadyspiraltsteady/ts,
                                                omegas=omegas,
                                                m=options.steadyspiralm,
                                                gamma=options.steadyspiralgamma)
        pot.append(steadyspiralp)
    if options.singletransientspiral:
        transientsp= TransientLogSpiralPotential(A=-0.035,
                                                 alpha=-7.,
                                                 omegas=.65,
                                                 sigma=options.transientspiralsigma,
                                                 to=options.transientspiralto)
        pot.append(transientsp)
    if options.elliptical:
        elp= EllipticalDiskPotential(cp=options.el_cp,
                                     sp=options.el_sp,
                                     p=options.el_p,
                                     tform=options.el_tform,
                                     tsteady=options.el_tsteady)
        pot.append(elp)
    if options.to is None:
        if options.bar and not options.steadyspiral \
                and not options.movingobject \
                and not options.transientspirals:
            to= barp.tform()
        else:
            to= -40.*2.*m.pi #Default= 40 Galactic rotations
    else:
        to= options.to
    if options.movingobject:
        o= Orbit([0.000001,.39,0.,6.,0.,m.pi/2.])
        o.integrate(nu.linspace(0,-to,10001),axip)
        o._orb.t+= to
        mp= MovingObjectPotential(o,GM=0.06,softening_model='plummer',
                                  softening_length=0.1)
        pot.append(mp)
    #Set-up df
    if options.dftype.lower() == 'dehnen':
        dfc= dehnendf(beta=options.beta,correct=True,niter=20,
                      profileParams=(options.rd,options.rs,options.so),
                      savedir='./')
    elif options.dftype.lower() == 'shu':
        dfc= shudf(beta=options.beta,correct=True,niter=20,
                   profileParams=(options.rd,options.rs,options.so),
                   savedir='./')
    edf= evolveddiskdf(dfc,pot,to=to)
    if options.nt == 1:
        evalts= nu.array([0.])
    else:
        evalts= nu.linspace(0.,to,options.nt)
    #Set up savefile
    if os.path.exists(options.savefilename):
        savefile= open(options.savefilename,'rb')
        surfmass= pickle.load(savefile)
        meanvr= pickle.load(savefile)
        meanvt= pickle.load(savefile)
        sigmar2= pickle.load(savefile)
        sigmat2= pickle.load(savefile)
        sigmart= pickle.load(savefile)
        vertexdev= pickle.load(savefile)
        surfmass_init= pickle.load(savefile)
        meanvt_init= pickle.load(savefile)
        sigmar2_init= pickle.load(savefile)
        sigmat2_init= pickle.load(savefile)
        ii= pickle.load(savefile)
        jj= pickle.load(savefile)
        grids= pickle.load(savefile)
        savefile.close()
    else:
        ii, jj= 0, 0
        surfmass= nu.zeros((nx,ny,options.nt))
        meanvr= nu.zeros((nx,ny,options.nt))
        meanvt= nu.zeros((nx,ny,options.nt))
        sigmar2= nu.zeros((nx,ny,options.nt))
        sigmat2= nu.zeros((nx,ny,options.nt))
        sigmart= nu.zeros((nx,ny,options.nt))
        vertexdev= nu.zeros((nx,ny,options.nt))
        surfmass_init= nu.zeros((nx,ny))
        meanvt_init= nu.zeros((nx,ny))
        sigmar2_init= nu.zeros((nx,ny))
        sigmat2_init= nu.zeros((nx,ny))
        grids= []
    #Loop through, saving
    while ii < nx:
        if not options.multi is None:
            sys.stdout.write('\r'+"X gridpoint %i out of %i" % \
                                 (ii+1,nx))
            sys.stdout.flush()
            multOut= multi.parallel_map(\
                lambda x: _calc_meanvel_single(x,ii,options,evalts,xgrid,ygrid,edf),
                range(ny),
                numcores=nu.amin([ny,multiprocessing.cpu_count(),
                                  options.multi]))
            for jj in range(ny):
                surfmass[ii,jj,:]= multOut[jj][0]
                meanvr[ii,jj,:]= multOut[jj][1]
                meanvt[ii,jj,:]= multOut[jj][2]
                sigmar2[ii,jj,:]= multOut[jj][3]
                sigmat2[ii,jj,:]= multOut[jj][4]
                sigmart[ii,jj,:]= multOut[jj][5]
                vertexdev[ii,jj,:]= multOut[jj][6]
                grids.append(multOut[jj][7])
                surfmass_init[ii,jj]= multOut[jj][8]
                meanvt_init[ii,jj]= multOut[jj][9]
                sigmar2_init[ii,jj]= multOut[jj][10]
                sigmat2_init[ii,jj]= multOut[jj][11]
            save_output(options.savefilename,surfmass,meanvr,meanvt,sigmar2,
                        sigmat2,sigmart,vertexdev,ii,jj,grids,
                        surfmass_init,meanvt_init,sigmar2_init,sigmat2_init)
            ii+= 1
            continue
        while jj < ny:
            if options.print_vprogress:
                sys.stdout.write("Spatial gridpoint %i out of %i\n" % \
                                     (jj+ii*ny+1,nx*ny))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r'+"Gridpoint %i out of %i" % \
                                     (jj+ii*ny+1,nx*ny))
                sys.stdout.flush()
            #Calculate R and phi
            if options.galcoords:
                R, phi= ygrid[jj], xgrid[ii]
            else:
                R= nu.sqrt((1.+xgrid[ii])**2.+ygrid[jj]**2.)
                phi= nu.arcsin(ygrid[jj]/R)
            #Calculate surfmass etc.
            smass, grid= edf.vmomentsurfacemass(R,0,0,grid=True,phi=phi,
                                                returnGrid=True,t=evalts,
                                                gridpoints=options.grid,
                                                nsigma=options.nsigma,
                                                hierarchgrid=options.hierarchgrid,
                                                nlevels=options.nlevels,
                                                print_progress=options.print_vprogress,
                                                integrate_method=options.integrate_method)
            if not options.dontsavegrid: grids.append(grid)
            surfmass[ii,jj,:]= smass
            meanvr[ii,jj,:]= edf.meanvR(R,phi=phi,grid=grid,t=evalts,
                                        surfacemass=surfmass[ii,jj,:],
                                        nsigma=options.nsigma,
                                        hierarchgrid=options.hierarchgrid,
                                        nlevels=options.nlevels)
            meanvt[ii,jj,:]= edf.meanvT(R,phi=phi,grid=grid,t=evalts,
                                        surfacemass=surfmass[ii,jj,:],
                                        nsigma=options.nsigma,
                                        hierarchgrid=options.hierarchgrid,
                                        nlevels=options.nlevels)
            sigmar2[ii,jj,:]= edf.sigmaR2(R,phi=phi,grid=grid,t=evalts,
                                          surfacemass=surfmass[ii,jj,:],
                                          meanvR=meanvr[ii,jj,:],
                                          nsigma=options.nsigma,
                                          hierarchgrid=options.hierarchgrid,
                                          nlevels=options.nlevels)
            sigmat2[ii,jj,:]= edf.sigmaT2(R,phi=phi,grid=grid,t=evalts,
                                          surfacemass=surfmass[ii,jj,:],
                                          meanvT=meanvt[ii,jj,:],
                                          nsigma=options.nsigma,
                                          hierarchgrid=options.hierarchgrid,
                                          nlevels=options.nlevels)
            sigmart[ii,jj,:]= edf.sigmaRT(R,phi=phi,grid=grid,t=evalts,
                                          surfacemass=surfmass[ii,jj,:],
                                          meanvR=meanvr[ii,jj,:],
                                          meanvT=meanvt[ii,jj,:],
                                          nsigma=options.nsigma,
                                          hierarchgrid=options.hierarchgrid,
                                          nlevels=options.nlevels)
            vertexdev[ii,jj,:]= edf.vertexdev(R,phi=phi,grid=grid,t=evalts,
                                              sigmaR2=sigmar2[ii,jj,:],
                                              sigmaT2=sigmat2[ii,jj,:],
                                              sigmaRT=sigmart[ii,jj,:],
                                              nsigma=options.nsigma,
                                              hierarchgrid=options.hierarchgrid,
                                              nlevels=options.nlevels)
            #Also calculate initial, non-trivial values
            use_init_grid= True #Much faster
            if use_init_grid:
                smass_init, init_grid= edf.vmomentsurfacemass(R,0,0,phi=phi,
                                                              t=edf._to,
                                                              grid=True,
                                                              gridpoints=options.grid,
                                                              returnGrid=True)
            else:
                smass_init= edf.vmomentsurfacemass(R,0,0,phi=phi,
                                                   t=edf._to)
                init_grid= False
            surfmass_init[ii,jj]= smass_init
            meanvt_init[ii,jj]= edf.meanvT(R,phi=phi,t=edf._to,
                                           grid=init_grid,
                                           surfacemass=surfmass_init[ii,jj],
                                           nsigma=options.nsigma)
            sigmar2_init[ii,jj]= edf.sigmaR2(R,phi=phi,t=edf._to,
                                             grid=init_grid,
                                             surfacemass=surfmass_init[ii,jj],
                                             meanvR=0.,nsigma=options.nsigma)
            sigmat2_init[ii,jj]= edf.sigmaT2(R,phi=phi,t=edf._to,
                                             grid=init_grid,
                                             surfacemass=surfmass_init[ii,jj],
                                             meanvT=meanvt_init[ii,jj],
                                             nsigma=options.nsigma)
            jj+= 1
            #Save
            save_output(options.savefilename,surfmass,meanvr,meanvt,sigmar2,
                        sigmat2,sigmart,vertexdev,ii,jj,grids,
                        surfmass_init,meanvt_init,sigmar2_init,sigmat2_init)
        ii+= 1
        jj= 0
    sys.stdout.write('\n')
    #Calculate Oort constants?
    if not options.oort is None and os.path.exists(options.oort):
        savefile= open(options.oort,'rb')
        oorta= pickle.load(savefile)
        oortb= pickle.load(savefile)
        oortc= pickle.load(savefile)
        oortk= pickle.load(savefile)
        oorta_init= pickle.load(savefile)
        oortb_init= pickle.load(savefile)
        oortc_init= pickle.load(savefile)
        oortk_init= pickle.load(savefile)
        oort_ii= pickle.load(savefile)
        oort_jj= pickle.load(savefile)
        derivRGrids= pickle.load(savefile)
        derivphiGrids= pickle.load(savefile)
        savefile.close()
    elif not options.oort is None:
        oort_ii, oort_jj= 0, 0
        oorta= nu.zeros((nx,ny,options.nt))
        oortb= nu.zeros((nx,ny,options.nt))
        oortc= nu.zeros((nx,ny,options.nt))
        oortk= nu.zeros((nx,ny,options.nt))
        oorta_init= nu.zeros((nx,ny))
        oortb_init= nu.zeros((nx,ny))
        oortc_init= nu.zeros((nx,ny))
        oortk_init= nu.zeros((nx,ny))
        derivRGrids= []
        derivphiGrids= []
    else: #to skip the loop and not mess up plotting
        oort_ii, oort_jj= nx, ny
        oorta= None
        oortb= None
        oortc= None
        oortk= None
        oorta_init= None
        oortb_init= None
        oortc_init= None
        oortk_init= None
        derivRGrids= []
        derivphiGrids= []
    #Loop through, saving
    while oort_ii < nx:
        while oort_jj < ny:
            if options.print_vprogress:
                sys.stdout.write("Spatial gridpoint %i out of %i\n" % \
                                 (oort_jj+oort_ii*ny+1,nx*ny))
                sys.stdout.flush()
            else:
                sys.stdout.write('\r'+"Gridpoint %i out of %i" % \
                                     (oort_jj+oort_ii*ny+1,nx*ny))
                sys.stdout.flush()
            #Calculate R and phi
            if options.galcoords:
                R, phi= ygrid[oort_jj], xgrid[oort_ii]
            else:
                R= nu.sqrt((1.+xgrid[oort_ii])**2.+ygrid[oort_jj]**2.)
                phi= nu.arcsin(ygrid[oort_jj]/R)
            if len(grids) > 0: grid= grids[oort_jj+oort_ii*ny]
            else: grid= True
            oorta[oort_ii,oort_jj,:], grid,derivrgrid, derivphigrid=\
                                      edf.oortA(R,phi=phi,
                                                grid=grid,
                                                derivRGrid=True,
                                                derivphiGrid=True,
                                                returnGrids=True,t=evalts,
                                                gridpoints=options.grid,
                                                derivGridpoints=options.grid,
                                                nsigma=options.nsigma,
                                                hierarchgrid=options.hierarchgrid,
                                                nlevels=options.nlevels,
                                                integrate_method=options.integrate_method)
            if not options.dontsavegrid:
                derivRGrids.append(derivrgrid)
                derivphiGrids.append(derivphigrid)
            oortb[oort_ii,oort_jj,:]= edf.oortB(R,phi=phi,t=evalts,
                                                grid=grid,
                                                derivRGrid=derivrgrid,
                                                derivphiGrid=derivphigrid,
                                                returnGrids=False,
                                                nsigma=options.nsigma,
                                                hierarchgrid=options.hierarchgrid,
                                                nlevels=options.nlevels)
            oortc[oort_ii,oort_jj,:]= edf.oortC(R,phi=phi,t=evalts,
                                                grid=grid,
                                                derivRGrid=derivrgrid,
                                                derivphiGrid=derivphigrid,
                                                returnGrids=False,
                                                nsigma=options.nsigma,
                                                hierarchgrid=options.hierarchgrid,
                                                nlevels=options.nlevels)
            oortk[oort_ii,oort_jj,:]= edf.oortK(R,phi=phi,t=evalts,
                                                grid=grid,
                                                derivRGrid=derivrgrid,
                                                derivphiGrid=derivphigrid,
                                                returnGrids=False,
                                                nsigma=options.nsigma,
                                                hierarchgrid=options.hierarchgrid,
                                                nlevels=options.nlevels)
            #Also calculate initial, non-trivial values
            use_init_grid= True #Much faster
            if use_init_grid:
                oorta_init[oort_ii,oort_jj], init_grid, init_derivrgrid,init_derivphigrid=\
                                             edf.oortA(R,phi=phi,
                                                       t=edf._to,
                                                       grid=True,
                                                       derivRGrid=True,
                                                       derivphiGrid=True,
                                                       derivGridpoints=options.grid,
                                                       gridpoints=options.grid,
                                                       returnGrids=True,
                                                       nsigma=options.nsigma)
                oortb_init[oort_ii,oort_jj]=\
                                              edf.oortB(R,phi=phi,
                                                       t=edf._to,
                                                       grid=init_grid,
                                                       derivRGrid=init_derivrgrid,
                                                       derivphiGrid=init_derivphigrid,
                                                        gridpoints=options.grid,
                                                       returnGrids=False,
                                                       nsigma=options.nsigma)
                oortc_init[oort_ii,oort_jj]=\
                                              edf.oortC(R,phi=phi,
                                                       t=edf._to,
                                                       grid=init_grid,
                                                       derivRGrid=init_derivrgrid,
                                                       derivphiGrid=init_derivphigrid,
                                                        gridpoints=options.grid,
                                                       returnGrids=False,
                                                       nsigma=options.nsigma)
                oortk_init[oort_ii,oort_jj]=\
                                              edf.oortK(R,phi=phi,
                                                       t=edf._to,
                                                       grid=init_grid,
                                                       derivRGrid=init_derivrgrid,
                                                       derivphiGrid=init_derivphigrid,
                                                        gridpoints=options.grid,
                                                       returnGrids=False,
                                                       nsigma=options.nsigma)
            else:
                oorta_init[oort_ii,oort_jj]= edf._initdf.oortA(R,nsigma=options.nsigma)
                oortb_init[oort_ii,oort_jj]= edf._initdf.oortB(R,nsigma=options.nsigma)
                oortc_init[oort_ii,oort_jj]= edf._initdf.oortC(R,nsigma=options.nsigma)
                oortk_init[oort_ii,oort_jj]= edf._initdf.oortK(R,nsigma=options.nsigma)
            oort_jj+= 1
            #Save
            save_pickles(options.oort,oorta,oortb,oortc,oortk,
                        oorta_init,oortb_init,oortc_init,oortk_init,
                        oort_ii,oort_jj,derivRGrids,derivphiGrids)
        oort_ii+= 1
        oort_jj= 0
    #print oorta, oortb, oortc, oortk
    if options.plot_resonance:
        resonance= []
        #Figure out where the resonances are
        if options.bar:
            OmegaP= barp.OmegaP()
            resonance.append(lindbladR(axip,OmegaP,m=-2.))
        if options.steadyspiral:
            OmegaP= steadyspiralp.OmegaP()
            resonance.append((lindbladR(axip,OmegaP,m=2.),False)) #False for no corot
            resonance.append((lindbladR(axip,OmegaP,m=-2.),False))
            resonance.append((lindbladR(axip,OmegaP,m=4.),False))
            resonance.append((lindbladR(axip,OmegaP,m=-4.),False))
            resonance.append((lindbladR(axip,OmegaP,m='corot'),True))
    else:
        resonance= None
    if options.dontplot:
        return None
    if options.movie:
        create_field_movie(options,surfmass,meanvr,meanvt,sigmar2,
                           sigmat2,sigmart,vertexdev,ii,jj,surfmass_init,
                           meanvt_init,sigmar2_init,sigmat2_init,evalts,to,
                           xgrid,ygrid,args[0],resonance,
                           oorta,oorta_init,
                           oortb,oortb_init,
                           oortc,oortc_init,
                           oortk,oortk_init)
    else:
        bovy_plot.bovy_print()
        plot_velocity_field(0,options,surfmass,meanvr,meanvt,sigmar2,
                            sigmat2,sigmart,vertexdev,ii,jj,surfmass_init,
                            meanvt_init,sigmar2_init,sigmat2_init,xgrid,ygrid,
                            resonance,
                            oorta,oorta_init,
                            oortb,oortb_init,
                            oortc,oortc_init,
                            oortk,oortk_init)
        bovy_plot.bovy_end_print(args[0])
        return None

def _calc_meanvel_single(jj,ii,options,evalts,xgrid,ygrid,edf):
    #Calculate R and phi
    if options.galcoords:
        R, phi= ygrid[jj], xgrid[ii]
    else:
        R= nu.sqrt((1.+xgrid[ii])**2.+ygrid[jj]**2.)
        phi= nu.arcsin(ygrid[jj]/R)
    #Calculate surfmass etc.
    smass, grid= edf.vmomentsurfacemass(R,0,0,grid=True,phi=phi,
                                        returnGrid=True,t=evalts,
                                        gridpoints=options.grid,
                                        nsigma=options.nsigma,
                                        hierarchgrid=options.hierarchgrid,
                                        nlevels=options.nlevels,
                                        print_progress=options.print_vprogress,
                                        integrate_method=options.integrate_method)
    surfmass= smass
    meanvr= edf.meanvR(R,phi=phi,grid=grid,t=evalts,
                       surfacemass=surfmass,
                       nsigma=options.nsigma,
                       hierarchgrid=options.hierarchgrid,
                       nlevels=options.nlevels)
    meanvt= edf.meanvT(R,phi=phi,grid=grid,t=evalts,
                       surfacemass=surfmass,
                       nsigma=options.nsigma,
                       hierarchgrid=options.hierarchgrid,
                       nlevels=options.nlevels)
    sigmar2= edf.sigmaR2(R,phi=phi,grid=grid,t=evalts,
                         surfacemass=surfmass,
                         meanvR=meanvr,
                         nsigma=options.nsigma,
                         hierarchgrid=options.hierarchgrid,
                         nlevels=options.nlevels)
    sigmat2= edf.sigmaT2(R,phi=phi,grid=grid,t=evalts,
                         surfacemass=surfmass,
                         meanvT=meanvt,
                         nsigma=options.nsigma,
                         hierarchgrid=options.hierarchgrid,
                         nlevels=options.nlevels)
    sigmart= edf.sigmaRT(R,phi=phi,grid=grid,t=evalts,
                         surfacemass=surfmass,
                         meanvR=meanvr,
                         meanvT=meanvt,
                         nsigma=options.nsigma,
                         hierarchgrid=options.hierarchgrid,
                         nlevels=options.nlevels)
    vertexdev= edf.vertexdev(R,phi=phi,grid=grid,t=evalts,
                             sigmaR2=sigmar2,
                             sigmaT2=sigmat2,
                             sigmaRT=sigmart,
                             nsigma=options.nsigma,
                             hierarchgrid=options.hierarchgrid,
                             nlevels=options.nlevels)
    #Also calculate initial, non-trivial values
    use_init_grid= True #Much faster
    if use_init_grid:
        smass_init, init_grid= edf.vmomentsurfacemass(R,0,0,phi=phi,
                                                      t=edf._to,
                                                      grid=True,
                                                      gridpoints=options.grid,
                                                      returnGrid=True)
    else:
        smass_init= edf.vmomentsurfacemass(R,0,0,phi=phi,
                                           t=edf._to)
        init_grid= False
    surfmass_init= smass_init
    meanvt_init= edf.meanvT(R,phi=phi,t=edf._to,
                            grid=init_grid,
                            surfacemass=surfmass_init,
                            nsigma=options.nsigma)
    sigmar2_init= edf.sigmaR2(R,phi=phi,t=edf._to,
                              grid=init_grid,
                              surfacemass=surfmass_init,
                              meanvR=0.,nsigma=options.nsigma)
    sigmat2_init= edf.sigmaT2(R,phi=phi,t=edf._to,
                              grid=init_grid,
                              surfacemass=surfmass_init,
                              meanvT=meanvt_init,
                              nsigma=options.nsigma)
    return [surfmass,meanvr,meanvt,sigmar2,sigmat2,sigmart,vertexdev,grid,
            surfmass_init,meanvt_init,sigmar2_init,sigmat2_init]

def create_field_movie(options,surfmass,meanvr,meanvt,sigmar2,
                       sigmat2,sigmart,vertexdev,ii,jj,surfmass_init,
                       meanvt_init,sigmar2_init,sigmat2_init,evalts,to,
                       xgrid,ygrid,moviefilename,resonance,
                       oorta,oorta_init,
                       oortb,oortb_init,
                       oortc,oortc_init,
                       oortk,oortk_init):
    import tempfile
    #First create all of the plots
    if options.tmpdir is None:
        tmpdir= '/tmp/'
    else:
        tmpdir= options.tmpdir
    tempdir= tempfile.mkdtemp(dir=tmpdir) #Temporary directory
    tmpfiles= []
    file_length= int(m.ceil(m.log10(options.nt)))
    #Create all frames
    if options.nt > 50: #Agg bug
        #Pickle
        picklethis= (options,surfmass,meanvr,meanvt,
                     sigmar2,
                     sigmat2,sigmart,vertexdev,ii,jj,surfmass_init,
                     meanvt_init,sigmar2_init,sigmat2_init,
                     xgrid,ygrid,resonance,
                     oorta,oorta_init,
                     oortb,oortb_init,
                     oortc,oortc_init,
                     oortk,oortk_init)
        picklefile= open(os.path.join(tempdir,'plotpickle.sav'),'wb')
        pickle.dump(picklethis,picklefile)
        pickle.dump(evalts,picklefile)
        pickle.dump(to,picklefile)
        picklefile.close()
    _NFRAMES= 50
    for ii in range(nu.amin([options.nt,_NFRAMES])):
        sys.stdout.write('\r'+"Working on frame %i out of %i" % \
                             (ii+1,options.nt))
        sys.stdout.flush()
        tmpfiles.append(os.path.join(tempdir,
                                     str(ii).zfill(file_length)))
        bovy_plot.bovy_print()
        plot_velocity_field(options.nt-ii-1,options,surfmass,meanvr,meanvt,
                            sigmar2,
                            sigmat2,sigmart,vertexdev,ii,jj,surfmass_init,
                            meanvt_init,sigmar2_init,sigmat2_init,
                            xgrid,ygrid,resonance,
                            oorta,oorta_init,
                            oortb,oortb_init,
                            oortc,oortc_init,
                            oortk,oortk_init)
        #Add time
        if options.plotpolar:
            bovy_plot.bovy_text(m.pi/4*6.,options.rmax*1.75,
                                r'$\mathrm{Age}\ = %4.1f \approx %4.1f\ \mathrm{Gyr}$' \
                                    % (evalts[options.nt-ii-1]-to,
                                       (evalts[options.nt-ii-1]-to)/2./m.pi/4.))
        else:
            bovy_plot.bovy_text(r'$\mathrm{Age}\ = %4.1f \approx %4.1f\ \mathrm{Gyr}$' \
                                    % (evalts[options.nt-ii-1]-to,
                                       (evalts[options.nt-ii-1]-to)/2./m.pi/4.),
                                top_left=True)

        bovy_plot.bovy_end_print(tmpfiles[ii]+'.png')
        #Convert to jpeg
        try:
            subprocess.check_call(['convert',
                                   tmpfiles[ii]+'.png',
                                   tmpfiles[ii]+'.jpg'])
        except subprocess.CalledProcessError:
            print "'convert' failed"
            raise subprocess.CalledProcessError
    #Do the rest
    extraSets= (options.nt-1)/_NFRAMES
    ii= _NFRAMES
    while extraSets > 0:
        thisend= nu.amin([options.nt,ii+_NFRAMES])
        #Call external routine
        subprocess.check_call([os.getenv('PYTHON'),
                               'plot_movie_frame.py',
                               os.path.join(tempdir,'plotpickle.sav'),
                               tempdir,
                               str(ii),
                               str(thisend),
                               str(file_length)])
        ii+= _NFRAMES
        extraSets-= 1
    sys.stdout.write('\n')
    #turn them into a movie
    try:
        #first pass
        subprocess.check_call(['ffmpeg',
                               '-pass','1',
                               '-r',str(options.framerate),
                               '-b', str(options.bitrate),
                               '-i',
                               os.path.join(tempdir,
                                            '%'+'0%id.jpg' % file_length),
                               '-y',
                               moviefilename])
        #2nd pass
        subprocess.check_call(['ffmpeg',
                               '-pass','2',
                               '-r',str(options.framerate),
                               '-b', str(options.bitrate),
                               '-i',
                               os.path.join(tempdir,
                                            '%'+'0%id.jpg' % file_length),
                               '-y',
                               moviefilename])
        """
        if thumbnail:
            thumbnameTemp= re.split(r'\.',filename)
            thumbnameTemp= thumbnameTemp[0:len(thumbnameTemp)-1]
            thumbname= ''
            for t in thumbnameTemp:
                thumbname+= t
            thumbname+= '.jpg'
            subprocess.check_call(['ffmpeg',
                                   '-itsoffset','-4','-y',
                                   '-i',filename,
                                   '-vcodec',
                                   'mjpeg',
                                   '-vframes','1',
                                   '-an',
                                   '-f',
                                   'rawvideo',
                                   '-s', '%ix%i' % (thumbsize,thumbsize),
                                   thumbname])
        """
    except subprocess.CalledProcessError:
        print "'ffmpeg' failed, scroll up to see ffmpeg's error message"
        raise subprocess.CalledProcessError
    finally:
        _cleanupMovieTempdir(tempdir)

def _cleanupMovieTempdir(tempdir):
    shutil.rmtree(tempdir)

def plot_velocity_field(indx,options,surfmass,meanvr,meanvt,sigmar2,
                        sigmat2,sigmart,vertexdev,ii,jj,surfmass_init,
                        meanvt_init,sigmar2_init,sigmat2_init,
                        xgrid,ygrid,resonance,
                        oorta,oorta_init,
                        oortb,oortb_init,
                        oortc,oortc_init,
                        oortk,oortk_init):
    if options.plottype.lower() == 'surfmass':
        if options.absolute:
            plotthis= surfmass[:,:,indx]
            plotthis_azavg= nu.sum(plotthis,axis=0)
            zeropoint= None
            #Also calculate maximum and minimum
            zlabel= r'$\Sigma$'
            z2label= r'$\langle\Sigma\langle$'
            vmin, vmax= 0., .2
        else:
            plotthis= surfmass[:,:,indx]/surfmass_init[:,:]
            plotthis_azavg= nu.sum(surfmass[:,:,indx],axis=0)/\
                nu.sum(surfmass_init[:,:],axis=0)
            zeropoint= [1.,1.]
            #Also calculate maximum and minimum
            vmin, vmax= .5, 1.5
            zlabel= r'$\Sigma / \Sigma^0$'
            z2label= r'$\langle\Sigma\rangle / \Sigma^0$'
    elif options.plottype.lower() == 'vr':
        plotthis= meanvr[:,:,indx]
        vmin, vmax= -20./235., 20./235.
        zlabel=r'$\bar{v}_R / v_0$'
        plotthis_azavg= nu.sum(plotthis*surfmass[:,:,indx],axis=0)/\
            nu.sum(surfmass[:,:,indx],axis=0)
        zeropoint= [0.,0.]
        z2label=r'$\langle\bar{v}_R\rangle / v_0$'
    elif options.plottype.lower() == 'vt':
        if options.absolute:
            plotthis= meanvt[:,:,indx]
            vmin, vmax= 1.-20./235.,1+20./235.
            zlabel=r'$\bar{v}_T / v_0$'
            plotthis_azavg= nu.sum(plotthis*surfmass[:,:,indx],axis=0)/\
                nu.sum(surfmass[:,:,indx],axis=0)
            zeropoint= [0.,0.]
            z2label=r'$\langle\bar{v}_T\rangle / v_0$'
        else:
            plotthis= meanvt[:,:,indx]-meanvt_init[:,:]
            vmin, vmax= -20./235.,20./235.
            zlabel=r'$\bar{v}_T - \bar{v}_T^0$'
            z2label=r'$\langle\bar{v}_T\rangle - \bar{v}_T^0$'
            plotthis_azavg= nu.sum(plotthis*surfmass[:,:,indx],axis=0)/\
                nu.sum(surfmass[:,:,indx],axis=0)
            zeropoint= [0.,0.]
    elif options.plottype.lower() == 'sigmar':
        if options.absolute:
            plotthis= nu.sqrt(sigmar2[:,:,indx])
            vmin, vmax= 0.,0.5
            zlabel=r'$\sigma_R / v_0$'
            z2label=r'$\langle\sigma_R\rangle / v_0$'
            plotthis_azavg= nu.sqrt((nu.sum(sigmar2[:,:,indx]\
                                                *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         +nu.sum(meanvr**2.[:,:,indx]\
                                                     *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         -(nu.sum(meanvr[:,:,indx]*surfmass[:,:,indx],axis=0)/\
                                               nu.sum(surfmass[:,:,indx],axis=0))**2.))
            zeropoint= None
        else:
            plotthis= nu.sqrt(sigmar2[:,:,indx]/sigmar2_init[:,:])
            vmin, vmax= .8, 1.2
            zlabel=r'$\sigma_R / \sigma_R^0$'
            z2label=r'$\langle\sigma_R\rangle / \sigma_R^0$'
            plotthis_azavg= nu.sqrt((nu.sum(sigmar2[:,:,indx]\
                                                *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         +nu.sum(meanvr[:,:,indx]**2.\
                                                     *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         -(nu.sum(meanvr[:,:,indx]*surfmass[:,:,indx],axis=0)/\
                                               nu.sum(surfmass[:,:,indx],axis=0))**2.)/\
                                        nu.mean(sigmar2_init[:,:],axis=0))
            zeropoint= [1.,1.]
    elif options.plottype.lower() == 'sigmat':
        if options.absolute:
            plotthis= nu.sqrt(sigmat2[:,:,indx])
            vmin, vmax= 0.,.5
            zlabel=r'$\sigma_T / v_0$'
            zlabel=r'$\langle\sigma_T\rangle / v_0$'
            plotthis_azavg= nu.sqrt((nu.sum(sigmat2[:,:,indx]\
                                                *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         +nu.sum(meanvt**2.[:,:,indx]\
                                                     *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         -(nu.sum(meanvt[:,:,indx]*surfmass[:,:,indx],axis=0)/\
                                               nu.sum(surfmass[:,:,indx],axis=0))**2.))
            zeropoint= None
        else:
            plotthis= nu.sqrt(sigmat2[:,:,indx]/sigmat2_init[:,:])
            vmin, vmax= .8,1.2
            zlabel=r'$\sigma_T / \sigma_T^0$'
            z2label=r'$\langle\sigma_T\rangle / \sigma_T^0$'
            plotthis_azavg= nu.sqrt((nu.sum(sigmat2[:,:,indx]\
                                                *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         +nu.sum(meanvt[:,:,indx]**2.\
                                                     *surfmass[:,:,indx],axis=0)/\
                                         nu.sum(surfmass[:,:,indx],axis=0)\
                                         -(nu.sum(meanvt[:,:,indx]*surfmass[:,:,indx],axis=0)/\
                                               nu.sum(surfmass[:,:,indx],axis=0))**2.)/\
                                        nu.mean(sigmat2_init[:,:],axis=0))
            zeropoint= [1.,1.]
    elif options.plottype.lower() == 'sigmart':
        plotthis= sigmart[:,:,indx]/nu.sqrt(sigmar2[:,:,indx]\
                                                *sigmat2[:,:,indx])
        vmin, vmax= -1.,1.
        zlabel=r'$\\rho_{RT} / v_0$'
    elif options.plottype.lower() == 'vertexdev':
        plotthis= vertexdev[:,:,indx]
        vmin, vmax= -45.,45.
        zlabel=r'$l_V\ [\mathrm{deg}]$'
    elif options.plottype.lower() == 'x2':
        if options.absolute:
            plotthis= sigmat2[:,:,indx]/sigmar2[:,:,indx]
            vmin, vmax= 0.3,0.7
            zlabel=r'$X^2$'
            z2label=r'$\langle\sigma_T/\sigma_R\rangle$'
            pltthis_azavg= None
            zeropoint= None
        else:
            plotthis= sigmat2[:,:,indx]/sigmar2[:,:,indx]\
                      *sigmar2_init[:,:]/sigmat2_init[:,:]
            vmin, vmax= .8, 1.2
            zlabel=r'$X^2/X_0^2$'
            z2label= None
            plotthis_azavg= None
            zeropoint= None
    elif options.plottype.lower() == 'oorta':
        if options.absolute:
            plotthis= oorta[:,:,indx]
            zeropoint= None
            #Also calculate maximum and minimum
            zlabel= r'$A$'
            #z2label= r'$\langle\Sigma\langle$'
            vmin, vmax= 0., 1.
        else:
            plotthis= oorta[:,:,indx]/oorta_init[:,:]
            zeropoint= [1.,1.]
            #Also calculate maximum and minimum
            vmin, vmax= .5, 1.5
            zlabel= r'$A / A^0$'
            #z2label= r'$\langle\Sigma\rangle / \Sigma^0$'
    elif options.plottype.lower() == 'oortb':
        if options.absolute:
            plotthis= oortb[:,:,indx]
            zeropoint= None
            #Also calculate maximum and minimum
            zlabel= r'$B$'
            #z2label= r'$\langle\Sigma\langle$'
            vmin, vmax= 0., 1.
        else:
            plotthis= oortb[:,:,indx]/oortb_init[:,:]
            zeropoint= [1.,1.]
            #Also calculate maximum and minimum
            vmin, vmax= .5, 1.5
            zlabel= r'$B / B^0$'
            #z2label= r'$\langle\Sigma\rangle / \Sigma^0$'
    elif options.plottype.lower() == 'oortc':
        plotthis= oortc[:,:,indx]
        vmin, vmax= -0.2, 0.2
        zlabel=r'$C$'
        #z2label=r'$\langle\bar{v}_R\rangle / v_0$'
    elif options.plottype.lower() == 'oortk':
        plotthis= oortk[:,:,indx]
        vmin, vmax= -0.2, 0.2
        zlabel=r'$K$'
        #z2label=r'$\langle\bar{v}_R\rangle / v_0$'
    elif options.plottype.lower() == 'vlos':
        #Build phi and l arrays
        phis= nu.zeros(meanvr[:,:,indx].shape)
        ls= nu.zeros(meanvr[:,:,indx].shape)
        for ii in range(len(xgrid)):
            phis[ii,:]= xgrid[ii]
            for jj in range(len(ygrid)):
                ls[ii,jj]= bovy_coords.rphi_to_dl_2d(ygrid[jj],
                                                     xgrid[ii],
                                                     degree=False,
                                                     ro=1.,
                                                     phio=0.)[1]
        #Then calculate vlos
        if options.absolute:
            plotthis= meanvt[:,:,indx]*nu.sin(phis+ls)\
                -meanvr[:,:,indx]*nu.cos(phis+ls)
            vmin, vmax= 1.-20./235.,1+20./235.
            zlabel=r'$\bar{v}_{\mathrm{los}} / v_0$'
        else:
            plotthis= (meanvt[:,:,indx]-meanvt_init[:,:])*nu.sin(phis+ls)\
                -meanvr[:,:,indx]*nu.cos(phis+ls)
            vmin, vmax= -20./235.,20./235.
            zlabel=r'$\bar{v}_{\mathrm{los}} / v_0-\bar{v}_{\mathrm{los}}^0 / v_0$'
    plotthis_azavg_min= nu.amin(plotthis,axis=0)
    plotthis_azavg_max= nu.amax(plotthis,axis=0)
    if not options.vmin is None:
        vmin= options.vmin
    if not options.vmax is None:
        vmax= options.vmax
    if not options.movie: print nu.amin(plotthis), nu.amax(plotthis)
    if options.galcoords and options.plotpolar:
        if options.azavg:
            fig= pyplot.figure()
            left, bottom = 0.1, 0.3
            width, height= 0.8, 0.6
            ax= fig.add_axes([left,bottom,width,height],
                             projection='galpolar') #galpolar is in bovy_plot
            #ax= pyplot.subplot(211,projection='galpolar')
        else:
            ax= pyplot.subplot(111,projection='polar')#galpolar is in bovy_plot
            ax.set_theta_direction(-1) #clockwise
        if not resonance is None:
            for r,corot in resonance:
                if r < 0.5 or r > 2.:
                    color= 'k'
                else:
                    color='w'
                if corot:
                    ls='-.'
                else:
                    ls= '--'
                ax.plot(nu.linspace(0.,2.*m.pi,501,),
                        nu.zeros(501)+r,ls=ls,color=color,zorder=5)
        plotxgrid= nu.linspace(xgrid[0]-(xgrid[1]-xgrid[0])/2.,
                               xgrid[-1]+(xgrid[1]-xgrid[0])/2.,
                               len(xgrid)+1)
        plotygrid= nu.linspace(ygrid[0]-(ygrid[1]-ygrid[0])/2.,
                               ygrid[-1]+(ygrid[1]-ygrid[0])/2.,
                               len(ygrid)+1)
        out= ax.pcolor(plotxgrid,plotygrid,plotthis.T,cmap='rainbow',
                       vmin=vmin,vmax=vmax)
        shrink= 0.8
        CB1= pyplot.colorbar(out,shrink=shrink)
        bbox = CB1.ax.get_position().get_points()
        CB1.ax.set_position(transforms.Bbox.from_extents(bbox[0,0]+0.025,
                                                         bbox[0,1],
                                                         bbox[1,0],
                                                         bbox[1,1]))
        CB1.set_label(zlabel)
        pyplot.ylim(0.,1.15*options.rmax)
        #Azimuthally averaged
        if options.azavg:
            left, bottom = 0.1, 0.1
            width, height= 0.8, 0.15
            ax2= fig.add_axes([left,bottom,width,height])
            #ax2= pyplot.subplot(210)
            bovy_plot.bovy_plot(ygrid,plotthis_azavg,'k-',overplot=True)
            if not zeropoint is None:
                bovy_plot.bovy_plot([0.,ygrid[-1]],zeropoint,'k--',
                                    overplot=True)
            bovy_plot.bovy_plot(ygrid,plotthis_azavg_min,color='0.5',ls='-',
                                overplot=True)
            bovy_plot.bovy_plot(ygrid,plotthis_azavg_max,color='0.5',ls='-',
                                overplot=True)
            if not resonance is None:
                for r in resonance:
                    bovy_plot.bovy_plot([r,r],[vmin,vmax],overplot=True,
                                        ls='--',color='0.5')
            pyplot.xlim(0.,ygrid[-1])
            pyplot.ylim(vmin,vmax)
            pyplot.xlabel(r'$R / R_0$')
            pyplot.ylabel(z2label)
            bovy_plot._add_ticks()
            #Set current axes back to polar ones for movie label
            pyplot.sca(ax)
    elif options.res == 1: #Just plot phi
        bovy_plot.bovy_plot(xgrid/nu.pi*180.,plotthis[:,0],'k-',
                            yrange=[vmin,vmax],
                            xrange=[0.,180.],
                            xlabel=r'$\mathrm{azimuth}\ [\mathrm{deg}]$',
                            ylabel=zlabel)
    else:
        bovy_plot.bovy_dens2d(plotthis.T,origin='lower',cmap='rainbow',#cmap='gist_yarg',
                              xrange=[-0.25,0.25],yrange=[-0.25,0.25],
                              xlabel=r'$X$',ylabel=r'$Y$',
                              interpolation='nearest',
                              vmin=vmin,vmax=vmax,contours=False,
                              colorbar=True,zlabel=zlabel,shrink=0.77)
        if not resonance is None:
            #phimin
            phimin= m.atan(0.25/.75)
            phis= nu.linspace(-phimin,phimin,1001)
            for r in resonance:
                bovy_plot.bovy_plot(1.-r*nu.cos(phis),r*nu.sin(phis),'w--',
                                    overplot=True)
def read_oort_output(savefilename):
    savefile= open(savefilename,'rb')
    oorta= pickle.load(savefile)
    oortb= pickle.load(savefile)
    oortc= pickle.load(savefile)
    oortk= pickle.load(savefile)
    oorta_init= pickle.load(savefile)
    oortb_init= pickle.load(savefile)
    oortc_init= pickle.load(savefile)
    oortk_init= pickle.load(savefile)
    oort_ii= pickle.load(savefile)
    oort_jj= pickle.load(savefile)
    derivRGrids= pickle.load(savefile)
    derivphiGrids= pickle.load(savefile)
    savefile.close()
    return (oorta,oortb,oortc,oortk,
            oorta_init,oortb_init,oortc_init,oortk_init,
            oort_ii,oort_jj,derivRGrids,derivphiGrids)
    
def read_output(savefilename):
    savefile= open(savefilename,'rb')
    surfmass= pickle.load(savefile)
    meanvr= pickle.load(savefile)
    meanvt= pickle.load(savefile)
    sigmar2= pickle.load(savefile)
    sigmat2= pickle.load(savefile)
    sigmart= pickle.load(savefile)
    vertexdev= pickle.load(savefile)
    surfmass_init= pickle.load(savefile)
    meanvt_init= pickle.load(savefile)
    sigmar2_init= pickle.load(savefile)
    sigmat2_init= pickle.load(savefile)
    ii= pickle.load(savefile)
    jj= pickle.load(savefile)
    grids= pickle.load(savefile)
    savefile.close()
    return (surfmass,meanvr,meanvt,sigmar2,sigmat2,sigmart,vertexdev,
            surfmass_init,meanvt_init,sigmar2_init,sigmat2_init,
            ii,jj,grids)

def save_output(savefilename,surfmass,meanvr,meanvt,sigmar2,
                sigmat2,sigmart,vertexdev,ii,jj,grids,
                surfmass_init,meanvt_init,sigmar2_init,sigmat2_init):
    saving= True
    interrupted= False
    while saving:
        try:
            savefile= open(savefilename,'wb')
            pickle.dump(surfmass,savefile)
            pickle.dump(meanvr,savefile)
            pickle.dump(meanvt,savefile)
            pickle.dump(sigmar2,savefile)
            pickle.dump(sigmat2,savefile)
            pickle.dump(sigmart,savefile)
            pickle.dump(vertexdev,savefile)
            pickle.dump(surfmass_init,savefile)
            pickle.dump(meanvt_init,savefile)
            pickle.dump(sigmar2_init,savefile)
            pickle.dump(sigmat2_init,savefile)
            pickle.dump(ii,savefile)
            pickle.dump(jj,savefile)
            pickle.dump(grids,savefile)
            savefile.close()
            saving= False
            if interrupted:
                raise KeyboardInterrupt
        except KeyboardInterrupt:
            if not saving:
                raise
            print "KeyboardInterrupt ignored while saving pickle ..."
            interrupted= True
    return None

def get_options():
    usage = "usage: %prog [options] <plotfilename>\n\nplotfilename= name of the file that the figure will be saved to"
    parser = OptionParser(usage=usage)
    parser.add_option("-s","--savefile",dest="savefilename",
                      default=None,
                      help="Name of the file that the field will be saved to")
    parser.add_option("--bar",action="store_true", 
                      default=False, dest="bar",
                      help="Use bar")
    parser.add_option("--steadyspiral",action="store_true", 
                      default=False, dest="steadyspiral",
                      help="Use steady spiral")
    parser.add_option("--transientspirals",action="store_true", 
                      default=False, dest="transientspirals",
                      help="Use transient spirals (NOT IMPLEMENTED YET)")
    parser.add_option("--singletransientspiral",action="store_true", 
                      default=False, dest="singletransientspiral",
                      help="Use a single transient spiral")
    parser.add_option("--movingobject",action="store_true", 
                      default=False, dest="movingobject",
                      help="Use a Moving Object like Quillen's")
    parser.add_option("--elliptical",action="store_true", 
                      default=False, dest="elliptical",
                      help="Use EllipticalDiskPotential")
    parser.add_option("-b","--beta",dest="beta",type='float',
                      default=0.,
                      help="Rotation curve power law index")
    parser.add_option("--to",dest="to",type='float',
                      default=None,
                      help="Time to when df was steady-state, default=bar formation")
    parser.add_option("-r",dest="res",type='int',
                      default=26,
                      help="Resolution of spatial grid")
    parser.add_option("--phires",dest="phires",type='int',
                      default=5,
                      help="Resolution of phi spatial grid (if phimax is set)")
    parser.add_option("--phimax",dest="phimax",type='float',
                      default=None,
                      help="Maximum (=-minimum) phi in rad (default=whole 2pi)")
    parser.add_option("--rmin",dest="rmin",type='float',
                      default=0.5,
                      help="Minimum Galactocentric radius")
    parser.add_option("--rmax",dest="rmax",type='float',
                      default=2.,
                      help="Maximum Galactocentric radius")
    parser.add_option("-g",dest="grid",type='int',
                      default=26,
                      help="Resolution of velocity grid")
    parser.add_option("--vmin",dest="vmin",type='float',
                      default=None,
                      help="Minimum value for density plot")
    parser.add_option("--vmax",dest="vmax",type='float',
                      default=None,
                      help="Maximum value for density plot")
    parser.add_option("--dftype", dest="dftype",
                      default='dehnen',
                      help="Type of DF to use")
    parser.add_option("--plottype", dest="plottype",
                      default='vr',
                      help="Type of plot to make (vr, vt, sigmar2, sigmat2, sigmart, vertexdev)")
    parser.add_option("--absolute",action="store_true", 
                      default=False, dest="absolute",
                      help="Plot 'absolute' values, otherwise relative to unevolved DF")
    parser.add_option("--altrect",action="store_true", 
                      default=False, dest="altrect",
                      help="Make the grid in an alternative rectangular grid")
    parser.add_option("--galcoords",action="store_true", 
                      default=False, dest="galcoords",
                      help="Make the grid in (R,azimuth) rather than (X,Y)")
    parser.add_option("--plotpolar",action="store_true", 
                      default=False, dest="plotpolar",
                      help="Plot in polar coordinates (with galcoords)")
    parser.add_option("--rd",dest="rd",type='float',
                      default=1./3.,
                      help="Initial disk scale-length")
    parser.add_option("--rs",dest="rs",type='float',
                      default=2.,
                      help="Initial disk sigma_R scale-length")
    parser.add_option("--so",dest="so",type='float',
                      default=31.4/220.,
                      help="Initial disk sigma_R at R_0")
    parser.add_option("--baralpha",dest="baralpha",type='float',
                      default=0.01,
                      help="Bar strength alpha parameter")
    parser.add_option("--bar_tform",dest="bar_tform",type='float',
                      default=-4.,
                      help="Bar formation time in bar periods")
    parser.add_option("--bar_tsteady",dest="bar_tsteady",type='float',
                      default=2.,
                      help="Bar steady time in bar periods")
    parser.add_option("--bar_olr",dest="bar_olr",type='float',
                      default=.9,
                      help="Bar radius of OLR")
    parser.add_option("--bar_angle",dest="bar_angle",type='float',
                      default=25.,
                      help="Angle between Sun--GC line and bar semi-major axis (deg)")
    parser.add_option("--el_cp",dest="el_cp",type='float',
                      default=0.05,
                      help="EllipticalDiskPotential cp")
    parser.add_option("--el_sp",dest="el_sp",type='float',
                      default=0.0,
                      help="EllipticalDiskPotential sp")
    parser.add_option("--el_p",dest="el_p",type='float',
                      default=0.0,
                      help="EllipticalDiskPotential p")
    parser.add_option("--dontsavegrid",action="store_true", 
                      default=False, dest="dontsavegrid",
                      help="Don't save the grids on which the velocity distribution is calculated")
    parser.add_option("--nt",dest="nt",type='int',
                      default=1,
                      help="Number of times to evaluate the DF at")
    parser.add_option("--movie",action="store_true", 
                      default=False, dest="movie",
                      help="create movie instead")
    parser.add_option("--tmpdir",dest="tmpdir",
                      default=None,
                      help="Temporary directory for movie frames")
    parser.add_option("--framerate",dest="framerate",type='int',
                      default=25,
                      help="Framerate for movie")
    parser.add_option("--bitrate",dest="bitrate",
                      default='512k',
                      help="bitrate for movie")
    parser.add_option("--nsigma",dest="nsigma",type='float',
                      default=None,
                      help="nsigma for df integration")
    parser.add_option("--transientspiralsigma",dest="transientspiralsigma",
                      default=5.,type='float',
                      help="'duration' of spiral wave (Gaussian sigma)")
    parser.add_option("--transientspiralto",dest="transientspiralto",
                      default=0.,type='float',
                      help="time at which the transient spiral peaks")
    parser.add_option("--steadyspiraltform",dest="steadyspiraltform",
                      default=None,type='float',
                      help="time at which the steady spiral forms (NOT in spiral periods)")
    parser.add_option("--steadyspiraltsteady",dest="steadyspiraltsteady",
                      default=None,type='float',
                      help="time at which the steady spiral becomes steady (NOT in spiral periods)")
    parser.add_option("--steadyspiralalpha",dest="steadyspiralalpha",
                      default=-12.5,type='float',
                      help="Alpha parameter of the steady spiral")
    parser.add_option("--steadyspiralomegas",dest="steadyspiralomegas",
                      default=0.65,type='float',
                      help="Pattern speed of the steady spiral")
    parser.add_option("--steadyspiralgamma",dest="steadyspiralgamma",
                      default=1.2,type='float',
                      help="Gamma parameter of the steady spiral")
    parser.add_option("--steadyspiralm",dest="steadyspiralm",
                      default=2,type='int',
                      help="Gamma parameter of the steady spiral")
    parser.add_option("--el_tform",dest="el_tform",
                      default=None,type='float',
                      help="time at which the elliptical disk forms)")
    parser.add_option("--el_tsteady",dest="el_tsteady",
                      default=None,type='float',
                      help="Time at which the elliptical disk becomes steady")
    parser.add_option("--hierarchgrid",action="store_true", 
                      default=False, dest="hierarchgrid",
                      help="Use a hierarchical grid")
    parser.add_option("--nlevels",dest="nlevels",type='int',
                      default=5,
                      help="Number of hierarchical levels (if --hierarchgrid is set)")
    parser.add_option("--azavg",action="store_true", 
                      default=False, dest="azavg",
                      help="Azimuthally average and plot as well")
    parser.add_option("--print_vprogress",action="store_true", 
                      default=False, dest="print_vprogress",
                      help="Print progress in calculating the velocity distribution for each spatial grid point")
    parser.add_option("--integrate_method",dest="integrate_method",
                      default="dopr54_c",
                      help="Orbit integration method")
    parser.add_option("--plot_resonance",action="store_true", 
                      default=False, dest="plot_resonance",
                      help="Plot the location of the main resonances")
    parser.add_option("--dontplot",action="store_true", 
                      default=False, dest="dontplot",
                      help="Don't plot")
    parser.add_option("--oort",dest="oort",
                      default=None,
                      help="(also) calculate the Oort constants, save them to this file in a manner similar to the first savefile")
    parser.add_option("-m","--multi",dest="multi",type='int',
                      default=None,
                      help="If set, use multiprocessing")
    return parser

if __name__ == '__main__':
    velocity_field(get_options())
