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
        data.setNcorrect(sf.vector_int([24, 32, 40,48, 50,48]))

        data.getIntensities()
        data.getNtrials()
        data.getNcorrect()
        data.getPcorrect()

        blocks = data.getNblocks()

        for i in range(blocks):
            data.getIntensity(i)
            data.getNtrials(i)
            data.getNcorrect(i)
            data.getPcorrect(i)
            data.getNoverK(i)

        data.getNalternatives()
        data.nonasymptotic()

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

    @staticmethod
    def generate_test_model():
        # IMPORTANT: here we can use the fact that PsiPsychometic manages its
        # own memory, and we don't need to hang on to th sigmoid, core, and
        # prior.
        return sf.PsiPsychometric(2, sf.abCore(), sf.PsiLogistic())

    def test_pschometric(self):
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        params = sf.vector_double([0.5,0.5,0.01])

        pr = sf.UniformPrior(0,1)
        psi.setPrior(0,pr)
        psi.setPrior(1,pr)
        psi.setPrior(2,pr)

        psi.evaluate(0.0,params)
        psi.negllikeli(params,data)
        psi.neglpost(params, data)
        psi.leastfavourable(params, data, 0.0)
        psi.deviance(params, data)
        psi.ddnegllikeli(params, data)
        psi.dnegllikeli(params, data)
        psi.getCore()
        psi.getSigmoid()

    def test_memory_management(self):
        core = sf.abCore()
        sigmoid = sf.PsiLogistic()
        psi = sf.PsiPsychometric(2,core,sigmoid)

    def test_exceptions(self):
        psi = TestPsychometric.generate_test_model()
        # for 2AFC we have 3 paramters with indices [0,1,2]
        self.assertRaises(ValueError, psi.setPrior,3, sf.UniformPrior(0,1))

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
        data = TestData.generate_test_dataset()
        psi = TestPsychometric.generate_test_model()
        cuts = sf.vector_double([1, 0.5])
        bs_list = sf.bootstrap(999, data, psi, cuts)


#x = numpy.arange(0,10,0.1)
#y = numpy.zeros(len(x))
#for i,val in enumerate(x):
#    y[i] = psi.evaluate(val,params)
#pylab.plot(x,y)
#pylab.show()

if __name__ == "__main__":
    ut.main()

