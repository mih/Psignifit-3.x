#!/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

import timeit
import _psipy as psipy
import swignifit.interface as sfi
import swignifit.swignifit_raw as sfr

x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
data = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]

def compare_time(name, sfi_time, psipy_time, repeat):
    print "time: " + name
    print "------------------------------"
    print 'swignifit total time: \t\t' , sfi_time
    print 'swignifit time per execution: \t' , sfi_time/repeat
    print 'psipy total time: \t\t' , psipy_time
    print 'psipy time per execution: \t' , psipy_time/repeat
    print 'percent swignfit slower: \t' , sfi_time/psipy_time * 100 - 100

def bootstrap_helper(wrapper):
    priors = ('flat','flat','Uniform(0,0.1)')
    sfr.setSeed(1)
    return wrapper.bootstrap(data,nsamples=10000,priors=priors)

def mcmc_helper(wrapper):
    priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')
    stepwidths = (1.,1.,0.01)
    sfr.setSeed(1)
    return wrapper.mcmc(data, nsamples=20000, priors=priors, stepwidths=stepwidths)

def time_bootstrap():
    repeat = 5
    t = timeit.Timer("pvs.bootstrap_helper(pvs.sfi)", "import psipy_vs_swignifit_time as pvs")
    sfi_time = t.timeit(number=repeat)
    t = timeit.Timer("pvs.bootstrap_helper(pvs.psipy)", "import psipy_vs_swignifit_time as pvs")
    psipy_time = t.timeit(number=repeat)
    compare_time('bootstrap', sfi_time, psipy_time, repeat)

def time_mcmc():
    repeat = 5
    t = timeit.Timer("pvs.mcmc_helper(pvs.sfi)", "import psipy_vs_swignifit_time as pvs")
    sfi_time = t.timeit(number=repeat)
    t = timeit.Timer("pvs.mcmc_helper(pvs.psipy)", "import psipy_vs_swignifit_time as pvs")
    psipy_time = t.timeit(number=repeat)
    compare_time('mcmc', sfi_time, psipy_time, repeat)

if __name__ == "__main__":
    "Will now compare execution time of psipy and swignifit, this may take a"+\
    "while"
    time_bootstrap()
    time_mcmc()
