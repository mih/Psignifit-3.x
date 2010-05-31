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
import timeit
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
    for i,name in enumerate(output_description):
        aaae(sfi_output[i], psipy_output[i], err_msg="swignifit and psipy "+\
                "differ for method: "+str(test_method) + "\n"+\
                "and output index: "+ str(i) + " named: "+name)

class TestBootstrap(ut.TestCase):

    @staticmethod
    def basic_helper(wrapper):
        priors = ('flat','flat','Uniform(0,0.1)')
        sfr.setSeed(1)
        return wrapper.bootstrap(data,nsamples=2000,priors=priors)

    def test_basic_correct(self):
        compare_wrappers(TestBootstrap.basic_helper,
            ["samples", "est", "D", "thres", "bias", "acc", "Rkd", "Rpd",
                "out", "influ"])

    def no_basic_time(self):
        t = timeit.Timer("pvs.TestBootstrap.basic_helper(pvs.sfi)", "import psipy_vs_swignifit as pvs")
        print 'swignifit:', t.timeit(number=5)
        t = timeit.Timer("pvs.TestBootstrap.basic_helper(pvs.psipy)", "import psipy_vs_swignifit as pvs")
        print 'psipy:', t.timeit(number=5)
        gc.enable()

class TestMCMC(ut.TestCase):

    @staticmethod
    def basic_helper(wrapper):
        priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')
        stepwidths = (1.,1.,0.01)
        # use a seed of 6 it fails sometimes
        # use a seed of 1 it won't fail, but instaed BOTH versions will give
        # different output occasionally
        sfr.setSeed(1)
        return wrapper.mcmc(data, nsamples=1000, priors=priors, stepwidths=stepwidths)

    def test_basic_correct(self):
        compare_wrappers(TestMCMC.basic_helper, ["estimates", "deviance",
            "posterior_predictive_data", "posterior_predictive_deviances",
            "posterior_predictive_Rpd", "posterior_predictive_Rkd",
            "logposterior_ratios"])
        #(estimates, deviance, posterior_predictive_data,
        #posterior_predictive_deviances, posterior_predictive_Rpd,
        #posterior_predictive_Rkd, logposterior_ratios)
        return wrapper.mcmc(data, nsamples=20) #, stepwidths=stepwidths)

    def test_fail_two_same_psipy(self):
        psipy_output1 = TestMCMC.basic_helper(psipy)
        psipy_output2 = TestMCMC.basic_helper(psipy)
        print psipy_output1[0]
        print psipy_output2[0]
        assert_output_equal(psipy_output1, psipy_output2)

    def test_fail_two_same_swignifit(self):
        sfi_output1 = TestMCMC.basic_helper(sfi)
        sfi_output2 = TestMCMC.basic_helper(sfi)
        print sfi_output1[0]
        print sfi_output2[0]
        assert_output_equal(sfi_output1, sfi_output2)

############################################################
# i have tried here to elaborate on some of the errors i have been getting
# but unfortunately this just makes the code undeterministic
# some of the tests fails sometimes

    def test_working(self):
        def helper(self):
            sfr.setSeed(5)
            return wrapper.mcmc(data, nsamples=20)
        sfi_output = TestMCMC.basic_helper(sfi)
        psipy_output = TestMCMC.basic_helper(psipy)
        assert_output_equal(sfi_output, psipy_output)

    def test_order_fail(self):
        def helper(wrapper):
            sfr.setSeed(6)
            return wrapper.mcmc(data, nsamples=20)
        # here we just switch the order
        sfi_output = helper(sfi)
        psipy_output = helper(psipy)
        assert_output_equal(sfi_output, psipy_output)

    def test_simple_fail(self):
        def helper(wrapper):
            sfr.setSeed(6)
            return wrapper.mcmc(data, nsamples=20)

        sfi_output = helper(sfi)
        psipy_output = helper(psipy)
        assert_output_equal(sfi_output, psipy_output)

    def no_basic_time(self):
        t = timeit.Timer("pvs.TestMCMC.basic_helper(pvs.sfi)", "import psipy_vs_swignifit as pvs")
        print 'swignifit per execution:', t.timeit(number=5)/5
        t = timeit.Timer("pvs.TestMCMC.basic_helper(pvs.psipy)", "import psipy_vs_swignifit as pvs")
        print 'psipy per execution:', t.timeit(number=5)/5
        #gc.enable()

class TestMapestimate(ut.TestCase):

    @staticmethod
    def basic_helper(wrapper):
        priors = ('flat','flat','Uniform(0,0.1)')
        return wrapper.mapestimate (data, priors=priors )

    def test_basic_correct(self):
        sfi_output = TestMapestimate.basic_helper(sfi)
        psipy_output = TestMapestimate.basic_helper(psipy)
        assert_output_equal(sfi_output, psipy_output)


if __name__ == "__main__":
    ut.main()

