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
import swignifit.swignifit_raw as sfr
import swignifit.interface as inter

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

class TestUtility(ut.TestCase):

    def test_make_dataset(self):
        dataset = inter.make_dataset(data, 1)
        self.assertTrue((np.array(x) == np.array(dataset.getIntensities())).all())
        self.assertTrue((np.array(k) == np.array(dataset.getNcorrect())).all())
        self.assertTrue((np.array(n) == np.array(dataset.getNtrials())).all())
        self.assertEqual(1, dataset.getNalternatives())

    def get_core(self):
        sigmoid = sfr.PsiLogistic()
        dataset = inter.make_dataset(data)
        ab = inter.get_core("ab", dataset, sigmoid)
        assertEqual("abCore", ab.__class__.__name__)
        mw = inter.get_core("mw", dataset, sigmoid)
        assertEqual("mwCore", mw.__class__.__name__)
        assertEqual(0.1, mw.getAlpha())
        mw = inter.get_core("mw0.2", dataset, sigmoid)
        assertEqual("mwCore", mw.__class__.__name__)
        assertEqual(0.2, mw.getAlpha())
        linear = inter.get_core("linear", dataset, sigmoid)
        assertEqual("linearCore", linear.__class__.__name__)
        log = inter.get_core("log", dataset, sigmoid)
        assertEqual("logCore", log.__class__.__name__)
        weibull = inter.get_core("weibull", dataset, sigmmoid)
        assertEqual("weibullCore", weibull.__class__.__name__)
        poly = inter.get_core("poly", dataset, sigmoid)
        assertEqual("polyCore", poly.__class__.__name__)

    def test_get_prior(self):
        uniform = inter.get_prior("Uniform(1,2)")
        self.assertEqual("UniformPrior", uniform.__class__.__name__)
        gauss = inter.get_prior("Gauss(0,1)")
        self.assertEqual("GaussPrior", gauss.__class__.__name__)
        beta = inter.get_prior("Beta(1.5, 3)")
        self.assertEqual("BetaPrior", beta.__class__.__name__)
        gamma = inter.get_prior("Gamma(1.5, 3)")
        self.assertEqual("GammaPrior", gamma.__class__.__name__)
        ngamma = inter.get_prior("nGamma(1.5,3)")
        self.assertEqual("nGammaPrior", ngamma.__class__.__name__)
        flat = inter.get_prior("flat")
        self.assertEqual(None, flat)
        unconstrained = inter.get_prior("unconstrained")
        self.assertEqual(None, unconstrained)

class TestBootstrap(ut.TestCase):

    def test_basic(self):
        inter.psibootstrap(data)

    def test_old_doctest(self):

        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('flat','flat','Uniform(0,0.1)')
        sfr.set_seed(1)
        samples,est,D,thres,bias,acc,Rkd,Rpd,out,influ = inter.psibootstrap(d,nsamples=2000,priors=priors)
        self.assertAlmostEqual( np.mean(est[:,0]), 2.7273945991794095)
        self.assertAlmostEqual( np.mean(est[:,1]), 1.3939511033770027)


    def test_start(self):
        inter.psibootstrap(data, nsamples=25, start=[0.1, 0.2, 0.3])

    def test_nsamples(self):
        inter.psibootstrap(data, nsamples=666)

    def test_nafc(self):
        inter.psibootstrap(data, nafc=23)

    def test_sigmoid(self):
        inter.psibootstrap(data, nsamples=25, sigmoid='gumbel_l')

    def test_core(self):
        inter.psibootstrap(data, nsamples=25, core='linear')

    def test_prior(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        inter.psibootstrap(data, nsamples=25, priors=priors)

    def test_cuts(self):
        inter.psibootstrap(data, nsamples=25, cuts=[0.5,0.6,0.75])

    def test_parameteric(self):
        inter.psibootstrap(data, nsamples=25, parametric=False)

class TestMCMC(ut.TestCase):

    def test_old_doctest(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')
        stepwidths = (1.,1.,0.01)
        sfr.set_seed(1)
        (estimates, deviance, posterior_predictive_data,
        posterior_predictive_deviances, posterior_predictive_Rpd,
        posterior_predictive_Rkd, logposterior_ratios) = inter.psimcmc(d,nsamples=10000,priors=priors,stepwidths=stepwidths)
        self.assertAlmostEqual( np.mean(estimates[:,0]), 2.5304815981388971)
        self.assertAlmostEqual( np.mean(estimates[:,1]), 1.6707238984255586)

class TestMapestimate(ut.TestCase):

    def test_old_doctest(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        priors = ('flat','flat','Uniform(0,0.1)')
        estimate, fisher, thres, deviance = inter.psimapestimate ( d, priors=priors )
        for i,value in enumerate([ 2.75183178, 1.45728231, 0.01555514]):
            self.assertAlmostEqual(value, estimate[i])
        print fisher
        self.assertAlmostEqual(2.75183178, thres)
        self.assertAlmostEqual(8.0713313969, deviance)

    def test_start(self):
        estimate, fisher, thres, deviance = inter.psimapestimate (data,
                start=[0.1, 0.2, 0.3])

class TestDiagnostics(ut.TestCase):

    def test_old_doctest(self):
        x = [float(2*k) for k in xrange(6)]
        k = [34,32,40,48,50,48]
        n = [50]*6
        d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
        prm = [2.75, 1.45, 0.015]
        pred,di,D,thres,Rpd,Rkd = inter.psidiagnostics(d,prm)
        self.assertAlmostEqual(8.07484858608, D)
        self.assertAlmostEqual(1.68932796526, di[0])
        self.assertAlmostEqual(-0.19344675783032761, Rpd)

    def test_intensities_only(self):
        prm = [2.75, 1.45, 0.015]
        predicted = inter.psidiagnostics(x, prm)

if __name__ == "__main__":
    ut.main()

