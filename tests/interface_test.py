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
import swignifit.utility as sfu

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

class TestUtility(ut.TestCase):

    def test_make_dataset(self):
        dataset = sfu.make_dataset(data, 1)
        self.assertTrue((np.array(x) == np.array(dataset.getIntensities())).all())
        self.assertTrue((np.array(k) == np.array(dataset.getNcorrect())).all())
        self.assertTrue((np.array(n) == np.array(dataset.getNtrials())).all())
        self.assertEqual(1, dataset.getNalternatives())

    def test_get_sigmoid(self):
        logistic = sfr.PsiLogistic()
        self.assertEqual("PsiLogistic", logistic.__class__.__name__)
        gauss = sfr.PsiGauss()
        self.assertEqual("PsiGauss", gauss.__class__.__name__)
        gumbel_l = sfr.PsiGumbelL()
        self.assertEqual("PsiGumbelL", gumbel_l.__class__.__name__)
        gumbel_r = sfr.PsiGumbelR()
        self.assertEqual("PsiGumbelR", gumbel_r.__class__.__name__)
        cauchy = sfr.PsiCauchy()
        self.assertEqual("PsiCauchy", cauchy.__class__.__name__)
        exponential = sfr.PsiExponential()
        self.assertEqual("PsiExponential", exponential.__class__.__name__)


    def test_get_core(self):
        sigmoid = sfr.PsiLogistic()
        dataset = sfu.make_dataset(data, 1)
        ab = sfu.get_core("ab", dataset, sigmoid)
        self.assertEqual("abCore", ab.__class__.__name__)
        mw = sfu.get_core("mw", dataset, sigmoid)
        self.assertEqual("mwCore", mw.__class__.__name__)
        self.assertEqual(0.1, mw.getAlpha())
        mw = sfu.get_core("mw0.2", dataset, sigmoid)
        self.assertEqual("mwCore", mw.__class__.__name__)
        self.assertEqual(0.2, mw.getAlpha())
        linear = sfu.get_core("linear", dataset, sigmoid)
        self.assertEqual("linearCore", linear.__class__.__name__)
        log = sfu.get_core("log", dataset, sigmoid)
        self.assertEqual("logCore", log.__class__.__name__)
        weibull = sfu.get_core("weibull", dataset, sigmoid)
        self.assertEqual("weibullCore", weibull.__class__.__name__)
        poly = sfu.get_core("poly", dataset, sigmoid)
        self.assertEqual("polyCore", poly.__class__.__name__)

    def test_get_prior(self):
        uniform = sfu.get_prior("Uniform(1,2)")
        self.assertEqual("UniformPrior", uniform.__class__.__name__)
        gauss = sfu.get_prior("Gauss(0,1)")
        self.assertEqual("GaussPrior", gauss.__class__.__name__)
        beta = sfu.get_prior("Beta(1.5, 3)")
        self.assertEqual("BetaPrior", beta.__class__.__name__)
        gamma = sfu.get_prior("Gamma(1.5, 3)")
        self.assertEqual("GammaPrior", gamma.__class__.__name__)
        ngamma = sfu.get_prior("nGamma(1.5,3)")
        self.assertEqual("nGammaPrior", ngamma.__class__.__name__)
        flat = sfu.get_prior("flat")
        self.assertEqual(None, flat)
        unconstrained = sfu.get_prior("unconstrained")
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

    def test_basic(self):
        inter.psimcmc(data)

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

    def test_start(self):
        inter.psimcmc(data,nsamples=25, start=[0.1,0.2,0.3])

    def test_nsamples(self):
        inter.psimcmc(data,nsamples=666)

    def test_nafc(self):
        inter.psimcmc(data,nsamples=25, nafc=23)

    def test_sigmoid(self):
        inter.psimcmc(data,nsamples=25, sigmoid='gumbel_r')

    def test_core(self):
        inter.psimcmc(data, nsamples=25, core='ab')

    def test_prior(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        inter.psimcmc(data, nsamples=25, priors=priors)

    def test_stepwidth(self):
        inter.psimcmc(data, nsamples=25, stepwidths=[0.1, 0.2, 0.3])

class TestMapestimate(ut.TestCase):

    def test_basic(self):
        inter.psimapestimate(data)

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

    def test_nafc(self):
        inter.psimapestimate(data, nafc=23)

    def test_sigmoid(self):
        inter.psimapestimate(data, sigmoid='gauss')

    def test_core(self):
        inter.psimapestimate(data, core='mw0.2')

    def test_priors(self):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        inter.psimapestimate(data, priors=priors)

    def test_cuts(self):
        inter.psimapestimate(data, cuts=[0.5, 0.75, 0.85])

    def test_start(self):
        estimate, fisher, thres, deviance = inter.psimapestimate (data,
                start=[0.1, 0.2, 0.3])

class TestDiagnostics(ut.TestCase):

    prm = [2.75, 1.45, 0.015]

    def test_basic(self):
        inter.psidiagnostics(data, TestDiagnostics.prm)

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

    def test_nafc(self):
        inter.psidiagnostics(data, TestDiagnostics.prm, nafc=23)

    def test_sigmoid(self):
        inter.psidiagnostics(data, TestDiagnostics.prm, sigmoid='logistic')

    def test_core(self):
        inter.psidiagnostics(data, TestDiagnostics.prm, core='linear')

    def test_cuts(self):
        inter.psidiagnostics(data, TestDiagnostics.prm, cuts=[0.5, 0.75, 0.85])

    def test_intensities_only(self):
        predicted = inter.psidiagnostics(x, TestDiagnostics.prm)

if __name__ == "__main__":
    ut.main()

