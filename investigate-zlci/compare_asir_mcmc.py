#!/usr/bin/env python

import pypsignifit as psi
import pylab as pl
import numpy as np
import scipy.stats
import pypsignifit.psigobservers as Obs
from pypsignifit.psigniplot import plotJoint

sigmoid = "logistic"
core    = "mw0.1"
priors  = ("Gauss(0,100)","invGamma(.5,.5)","Beta(2,30)")
Beta = scipy.stats.beta
nafc = 2

def compare_on_dataset ( d, nafc=2 ):
    for k in xrange ( 1 ):
        A = psi.ASIRInference ( d, nafc=nafc, priors=priors, sigmoid=sigmoid, core=core )

        Am025,Am975 = pl.prctile ( A.mcestimates[:,0], (2.5,97.5) )
        Aw025,Aw975 = pl.prctile ( A.mcestimates[:,1], (2.5,97.5) )
        Amm,Awm = A.mcestimates[:,:2].mean(0)
        Ams,Aws = A.mcestimates[:,:2].std(0)

        M = psi.BayesInference ( d, nafc=nafc, priors=priors, sigmoid=sigmoid, core=core )
        M.sample( start=M.farstart )
        M.sample( start=M.farstart )

        Mm025,Mm975 = pl.prctile ( M.mcestimates[:,0], (2.5,97.5) )
        Mw025,Mw975 = pl.prctile ( M.mcestimates[:,1], (2.5,97.5) )
        Mmm,Mwm = M.mcestimates[:,:2].mean(0)
        Mms,Mws = M.mcestimates[:,:2].std(0)
        MR = M.Rhat ()

        print "ASIR:",Am025,"MCMC:",Mm025, "R^:",MR

    psi.ConvergenceMCMC ( M, 0 )
    psi.ConvergenceMCMC ( M, 1 )
    psi.ConvergenceMCMC ( M, 2 )
    plotJoint ( A )

def sampling_scheme ( observer, nblocks ):
    if observer.model["nafc"] == 2:
        B = Beta ( 1.5,.6 )
    elif observer.model["nafc"] == 1:
        B = Beta ( .5,.5 )
    Fx = B.ppf ( np.mgrid[.025:.975:nblocks*1j] )
    return O.getlevels ( Fx )

if __name__ == "__main__":
    nblocks, ntrials = 6,20
    O = Obs.Observer ( 4,2,.03, sigmoid=sigmoid, core=core, nafc=nafc )
    d = np.array ( O.DoAnExperiment ( sampling_scheme(O,nblocks), ntrials ) )

    compare_on_dataset ( d, O.model["nafc"] )

    pl.show()
