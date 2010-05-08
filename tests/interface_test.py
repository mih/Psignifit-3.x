#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" Unit Tests for interface to swig wrapped methods """

import numpy as np
import unittest as ut
import swignifit as sf
import swignifit.interface as inter

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

class TestBootstrap(ut.TestCase):

    def test_basic(self):
        inter.bootstrap(data)

    def test_old_doctest(self):

        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('flat','flat','Uniform(0,0.1)')
        sf.set_seed(1)
        samples,est,D,thres,bias,acc,Rkd,Rpd,out,influ = inter.bootstrap(d,nsamples=2000,priors=priors)
        self.assertAlmostEqual( np.mean(est[:,0]), 2.7273945991794095)
        self.assertAlmostEqual( np.mean(est[:,1]), 1.3939511033770027)


    def test_start(self):
        inter.bootstrap(data, nsamples=25, start=[0.1, 0.2, 0.3])

    def test_nsamples(self):
        inter.bootstrap(data, nsamples=666)

    def test_nafc(self):
        inter.bootstrap(data, nafc=23)

    def test_sigmoid(self):
        inter.bootstrap(data, nsamples=25, sigmoid='gumbel_l')

    def test_core(self):
        inter.bootstrap(data, nsamples=25, core='linear')

    def test_prior(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        inter.bootstrap(data, nsamples=25, priors=priors)

    def test_cuts(self):
        inter.bootstrap(data, nsamples=25, cuts=[0.5,0.6,0.75])

    def test_parameteric(self):
        inter.bootstrap(data, nsamples=25, parametric=False)

class TestMCMC(ut.TestCase):

    def test_old_doctest(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')
        stepwidths = (1.,1.,0.01)
        sf.set_seed(1)
        (estimates, deviance, posterior_predictive_data,
        posterior_predictive_deviances, posterior_predictive_Rpd,
        posterior_predictive_Rkd, logposterior_ratios) = inter.psimcmc(d,nsamples=10000,priors=priors,stepwidths=stepwidths)
        self.assertAlmostEqual( np.mean(estimates[:,0]), 2.5304815981388971)
        self.assertAlmostEqual( np.mean(estimates[:,1]), 1.6707238984255586)


if __name__ == "__main__":
    ut.main()

