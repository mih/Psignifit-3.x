#!/usr/bin/env python

import reader
import sys
import pypsignifit as pf
import pylab as pl
import numpy as np
import swignifit.swignifit_raw as sfr
# import pypsignifit.psigniplot as pp

d,s = reader.read_data ( sys.argv[1] )

model = {'nafc':1, 'sigmoid':"logistic", 'core':'mw0.1'}
m = 4.0
w = 4.0
l = 0.05
g = 0.02

priors = ["Gauss(%f,%f)" % (m, m), "Gauss(%f,%f)" % (w, w*2), "Beta(2,50)",
"Beta(1,50)"]
mcmc = pf.BayesInference(d, sigmoid=model['sigmoid'], core=model['core'],
        nafc=model['nafc'], priors=priors, verbose=True)

# th = np.mgrid[3:5:20j,0:10:20j,0:.1:20j,0:.1:20j]
# th_ = np.reshape(th,(4,-1))
# post = np.zeros ( th_.shape[-1] )
# lapl = np.zeros ( th_.shape[-1] )
# 
# def laplace ( prm ):
#     s = mcmc._steps
#     m = mcmc.mapestimate
#     return np.sum(((prm-m)/s)**2)
# 
# for i in xrange ( len(post) ):
#     lapl[i] = laplace ( th_[:,i] )
#     prm = sfr.vector_double ( th_[:,i] )
#     post[i] = mcmc._pmf.neglpost ( prm, mcmc._data )
# post = np.reshape ( post, (20,20,20,20) )
# lapl = np.reshape ( lapl, (20,20,20,20) )
# 
# def plt ( slices ):
#     pl.subplot(211)
#     pl.imshow ( post[slices] )
#     pl.subplot(212)
#     pl.imshow ( lapl[slices] )
# 
# plt ( (10,10,slice(0,-1),slice(0,-1)) )
# pl.show()

print mcmc.approx

# pl.subplot(411)
# pl.plot ( mcmc.debug_samples[0][:,0] )
# pl.subplot(412)
# pl.plot ( mcmc.debug_samples[0][:,1] )
# pl.subplot(413)
# pl.plot ( mcmc.debug_samples[0][:,2] )
# pl.subplot(414)
# pl.plot ( mcmc.debug_samples[0][:,3] )
# pl.show()

# pf.ConvergenceMCMC(mcmc,0)
# pf.ConvergenceMCMC(mcmc,1)
# pf.ConvergenceMCMC(mcmc,2)
# pf.show()
