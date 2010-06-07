#!/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" Script to compare psipy and swignifit for correctness and speed """

import _psipy as psipy
import swignifit.swignifit_raw as sfr
import swignifit.interface as sfi
import unittest as ut
import time
import nose.tools as nt
from numpy.testing import assert_array_almost_equal as aaae

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

def print_time(name, sec):
    """ print seconds in terms of hours, minutes and secondd """
    hours, remainder = divmod(sec, 3600)
    minutes, seconds = divmod(remainder, 60) 
    print '%s took %d hours %d minutes and %d seconds' % (name, hours, minutes, seconds)

def assert_output_equal(o1, o2):
    nt.assert_equal(len(o1), len(o2))
    for i in xrange(len(o1)):
        #print i
        #print o1[i]
        #print o2[i]
        aaae(o1[i], o2[i])

def compare_wrappers(test_method, output_description):
    sfi_output = test_method(sfi)
    psipy_output = test_method(psipy)
    nt.assert_equal(len(sfi_output), len(output_description))
    nt.assert_equal(len(psipy_output), len(output_description))
    compare_output(sfi_output, psipy_output, output_description)

def compare_output(output1, output2, output_description):
    for i,name in enumerate(output_description):
        aaae(output1[i], output2[i], err_msg="two outputs "+\
                "differ for output index: "+ str(i) + " named: "+name)

class TestBootstrap(ut.TestCase):

    output_description = ["samples", "est", "D", "thres", "bias", "acc",
            "Rkd", "Rpd", "out", "influ"]

    @staticmethod
    def basic_helper(wrapper):
        priors = ('flat','flat','Uniform(0,0.1)')
        sfr.setSeed(1)
        return wrapper.bootstrap(data,nsamples=2000,priors=priors)

    @staticmethod
    def extended_helper(wrapper):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        sfr.setSeed(1)
        return wrapper.bootstrap(data, start=[0.1, 0.2, 0.3], nsamples=100, nafc=4,
                sigmoid="gumbel_l", core="linear", priors=priors,
                cuts=[0.5,0.6,0.75])


    def test_basic(self):
        compare_wrappers(TestBootstrap.basic_helper,
                TestBootstrap.output_description)

    def test_extended(self):
        compare_wrappers(TestBootstrap.extended_helper,
                TestBootstrap.output_description)

class TestMCMC(ut.TestCase):
    output_description = ["estimates", "deviance",
            "posterior_predictive_data", "posterior_predictive_deviances",
            "posterior_predictive_Rpd", "posterior_predictive_Rkd",
            "logposterior_ratios"]

    @staticmethod
    def basic_helper(wrapper):
        priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')
        stepwidths = (1.,1.,0.01)
        sfr.setSeed(1)
        return wrapper.mcmc(data, nsamples=1000, priors=priors, stepwidths=stepwidths)

    @staticmethod
    def extended_helper(wrapper):
        stepwidths = (1.,1.,0.01)
        sfr.setSeed(1)
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        return wrapper.mcmc(data, start=[0.1,0.2,0.3], nsamples=1000, nafc=4,
                sigmoid="gumbel_r", core="ab", priors=priors,
                stepwidths=stepwidths)

    def test_basic(self):
        def helper(wrapper):
            sfr.setSeed(1)
            return wrapper.mcmc(data, nsamples=20)
        compare_wrappers(helper, TestMCMC.output_description)

    def test_extended(self):
        compare_wrappers(TestMCMC.extended_helper,
                TestMCMC.output_description)

    # we have so many tests here since we were once debugging a nasty error in
    # the RNG the cause of which was long unknowen

    def test_with_more_params(self):
        compare_wrappers(TestMCMC.basic_helper, TestMCMC.output_description)

    def test_with_more_samples(self):
        def helper(wrapper):
            sfr.setSeed(1)
            return wrapper.mcmc(data, nsamples=20000)
        print "Testing mcmc with 20000, this may take a few seconds"
        compare_wrappers(helper, TestMCMC.output_description)

    def test_two_same_psipy(self):
        psipy_output1 = TestMCMC.basic_helper(psipy)
        psipy_output2 = TestMCMC.basic_helper(psipy)
        compare_output(psipy_output1, psipy_output2, TestMCMC.output_description)

    def test_two_same_swignifit(self):
        sfi_output1 = TestMCMC.basic_helper(sfi)
        sfi_output2 = TestMCMC.basic_helper(sfi)
        compare_output(sfi_output1, sfi_output2, TestMCMC.output_description)

    def test_different_seed(self):
        def helper(wrapper):
            sfr.setSeed(6)
            return wrapper.mcmc(data, nsamples=20)
        compare_wrappers(helper, TestMCMC.output_description)

    def test_order_fail(self):
        # here we just switch the order
        sfi_output = TestMCMC.basic_helper(sfi)
        psipy_output = TestMCMC.basic_helper(psipy)
        compare_output(sfi_output, psipy_output, TestMCMC.output_description)

        psipy_output = TestMCMC.basic_helper(psipy)
        sfi_output = TestMCMC.basic_helper(sfi)
        compare_output(psipy_output, sfi_output, TestMCMC.output_description)

class TestMapestimate(ut.TestCase):

    output_description = ["estimate", "fisher", "thres", "deviance"]

    @staticmethod
    def basic_helper(wrapper):
        priors = ('flat','flat','Uniform(0,0.1)')
        return wrapper.mapestimate (data, priors=priors )

    @staticmethod
    def extended_helper(wrapper):
        priors = ('Gauss(0,10)', 'Gamma(2,3)', 'Uniform(1,5)')
        return wrapper.mapestimate(data, sigmoid='gauss', core='mw0.2',
                priors=priors, cuts=[0.5, 0.75, 0.85], start=[0.1, 0.2, 0.3])

    def test_basic(self):
        compare_wrappers(TestMapestimate.basic_helper, TestMapestimate.output_description)

    def test_extended(self):
        compare_wrappers(TestMapestimate.extended_helper, TestMapestimate.output_description)

class TestDiagnostics(ut.TestCase):

    output_description = ["pred" ,"di" ,"D", "thres", "Rpd", "Rkd"]
    prm = [2.75, 1.45, 0.015]

    @staticmethod
    def basic_helper(wrapper):
        return wrapper.diagnostics(data, TestDiagnostics.prm)

    @staticmethod
    def extended_helper(wrapper):
        return wrapper.diagnostics(data, TestDiagnostics.prm, nafc=4,
                sigmoid="cauchy", core="weibull", cuts=[0.5, 0.75, 0.85])

    @staticmethod
    def intensities_helper(wrapper):
        return wrapper.diagnostics(x, TestDiagnostics.prm)

    def test_basic(self):
        compare_wrappers(TestDiagnostics.basic_helper,
                TestDiagnostics.output_description)

    def test_intensities(self):
        compare_wrappers(TestDiagnostics.intensities_helper,
                TestDiagnostics.output_description)

if __name__ == "__main__":
    #s = ut.TestSuite()
    #s.addTest(TestMCMC("test_fail_two_same_psipy"))
    #s.addTest(TestMCMC("test_fail_two_same_swignifit"))
    #ut.TextTestRunner(verbosity=2).run(s)
    ut.main()

