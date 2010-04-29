#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" Unit Tests for interface to swig wrapped methods """

import numpy, pylab
import unittest as ut
import swignifit as sf
import swignifit.interface as inter

class TestBootstrap(ut.TestCase):
    def setUp(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        self.d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

    def test_basic(self):
        inter.bootstrap(self.d)

    def test_old_doctest(self):

        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('flat','flat','Uniform(0,0.1)')
        samples,est,D,thres,bias,acc,Rkd,Rpd,out,influ = inter.bootstrap(d,nsamples=2000,priors=priors)
        self.assertAlmostEqual( numpy.mean(est[:,0]), 2.7762481672120902)
        self.assertAlmostEqual( numpy.mean(est[:,1]), 1.4243919674602623)


    def test_start(self):
        inter.bootstrap(self.d, nsamples=25, start=[0.1, 0.2, 0.3])

    def test_nsamples(self):
        inter.bootstrap(self.d, nsamples=666)

    def test_nafc(self):
        inter.bootstrap(self.d, nafc=23)

    def test_sigmoid(self):
        inter.bootstrap(self.d, nsamples=25, sigmoid='gumbel_l')

    def test_core(self):
        inter.bootstrap(self.d, nsamples=25, core='linear')

    def test_prior(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        inter.bootstrap(self.d, nsamples=25, priors=priors)

    def test_cuts(self):
        inter.bootstrap(self.d, nsamples=25, cuts=[0.5,0.6,0.75])

    def test_parameteric(self):
        inter.bootstrap(self.d, nsamples=25, parametric=False)

if __name__ == "__main__":
    ut.main()

