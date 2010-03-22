#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" Unit Tests for raw swig wrapper """

import numpy, pylab
import unittest as ut
import swignifit as sf

class TestSigmoid(ut.TestCase):
    """ test that all sigmoids have been wrapped and can be executed """

    def all_methods(self, sigmoid):
        sigmoid.f(0.0)
        sigmoid.df(0.0)
        sigmoid.ddf(0.0)
        sigmoid.inv(0.1)

    def test_cauchy(self):
        self.all_methods(sf.PsiCauchy())

    def test_exponential(self):
        self.all_methods(sf.PsiExponential())

    def test_gauss(self):
        self.all_methods(sf.PsiGauss())

    def test_gumbell(self):
        self.all_methods(sf.PsiGumbelL())

    def test_gumbelr(self):
        self.all_methods(sf.PsiGumbelR())

    def test_logistic(self):
        self.all_methods(sf.PsiLogistic())

    def test_exponential_exception(self):
        s = sf.PsiExponential()
        self.assertRaises(ValueError, s.inv, 0)
        self.assertRaises(ValueError, s.inv, 1)

class TestData(ut.TestCase):

    @staticmethod
    def generate_test_dataset():
        x = sf.vector_double([0.,2.,4.,6., 8., 10.])
        k = sf.vector_int([24, 32, 40,48, 50,48])
        n = sf.vector_int(6*[50])
        return sf.PsiData(x, n, k, 2)

    def test_data(self):
        data = TestData.generate_test_dataset()

class TestCore(ut.TestCase):

    data = TestData.generate_test_dataset()

    def all_methods(self, core):
        params = sf.vector_double([1.0,1.0])
        core.g(0.0, params)
        core.dg(0.0,params,0)
        core.dg(0.0,params,1)
        core.ddg(0.0,params,0,0)
        core.ddg(0.0,params,0,1)
        core.ddg(0.0,params,1,0)
        core.ddg(0.0,params,1,1)
        core.inv(0.0,params)
        core.dinv(0.0,params,0)
        core.dinv(0.0,params,1)
        core.transform(2,1.0,1.0)

    def test_ab_core(self):
        self.all_methods(sf.abCore())

    def test_linear_core(self):
        self.all_methods(sf.linearCore())

    def test_log_core(self):
        self.all_methods(sf.logCore(TestCore.data))

    def test_mw_core(self):
        # mwCore constructor is a bit different than the rest
        self.all_methods(sf.mwCore(1))

    def test_poly_core(self):
        self.all_methods(sf.polyCore(TestCore.data))

    def test_weibull_core(self):
        self.all_methods(sf.weibullCore(TestCore.data))

    def test_exceptions(self):
        c = sf.logCore(TestCore.data)
        params = sf.vector_double([1.0,1.0])
        self.assertRaises(ValueError, c.g, -1.0, params)
        c = sf.weibullCore(TestCore.data)
        self.assertRaises(ValueError, c.dg, -1.0, params, 0)
        self.assertRaises(ValueError, c.ddg, -1.0, params, 0, 1)


class TestPsychometric(ut.TestCase):

    def test_pschometric(self):
        core = sf.abCore()
        sigmoid = sf.PsiLogistic()
        psi = sf.PsiPsychometric(2,core,sigmoid)
        params = sf.vector_double([0.5,0.5,0.01])
        data = TestData.generate_test_dataset()

        psi.evaluate(0.0,params)
        psi.negllikeli(params,data)
        psi.neglpost(params, data)
        psi.leastfavourable(params, data, 0.0)
        psi.deviance(params, data)
        psi.ddnegllikeli(params, data)
        psi.dnegllikeli(params, data)
        psi.getCore()
        psi.getSigmoid()

class TestPriors(ut.TestCase):

    def all_methods(self, prior):
        prior.pdf(0.0)
        prior.dpdf(0.0)
        prior.rand()

    def test_beta_prior(self):
        self.all_methods(sf.BetaPrior(1.5, 3))

    def test_gamma_prior(self):
        self.all_methods(sf.GammaPrior(1.5, 3))

    def test_ngamma_prior(self):
        self.all_methods(sf.nGammaPrior(1.5, 3))

    def test_gauss_prior(self):
        self.all_methods(sf.GaussPrior(1.5, 3))

    def test_uniform_prior(self):
        self.all_methods(sf.UniformPrior(1.5, 3))

class TestBootstrap(ut.TestCase):

    def test_bootstrap(self):
        # IMPORTANT we need to assign local variables for these objects.  If we
        # simply do pmf.setPrior(2,sf.UniformPrior(0.0, 0.1)) the
        # UniformPrior will be out of scope and deleted by the time we reach
        # bootstrap, resulting in a segfault.
        data = TestData.generate_test_dataset()
        core = sf.abCore()
        prior = sf.UniformPrior(0.0, 0.1)
        sigmoid = sf.PsiLogistic()
        cuts = sf.vector_double([1, 0.5])

        pmf = sf.PsiPsychometric(2, core, sigmoid)
        pmf.setPrior(2, prior)
        bs_list = sf.bootstrap(999, data, pmf, cuts)
        print bs_list


#x = numpy.arange(0,10,0.1)
#y = numpy.zeros(len(x))
#for i,val in enumerate(x):
#    y[i] = psi.evaluate(val,params)
#pylab.plot(x,y)
#pylab.show()

if __name__ == "__main__":
    ut.main()

