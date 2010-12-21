#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__docformat__ = "restructuredtext"

import sys
import subprocess

# This is the interface to psi++
import swignifit as interface

# This is "real" psignifit
from psignidata import *
from psigniplot import *

# This is to enable display of graphics
from pylab import show

interface.set_seed( 0 )

def set_seed(value):
    interface.set_seed(value)

# add git hash to options, to automate the process of storing it
# note that when not running this via PYTHONPATH from the git repo
# this will fail, instead it should query pypsignifit for a proper version
# number, however proper version numbers were not implemented at the time of
# writing (see git blame for the exact date.)
def __get_command_output(command_string, env):
    """ Execute arbitrary commands.

    Parameters
    ----------
    command_list : strings
        the command and its arguments
    env: mapping
        mapping ov environment variables to values

    Returns
    -------
    output : string
        the raw output of the command executed
    """

    command_list = command_string.split(' ')
    p = subprocess.Popen(command_list, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, env=env)
    p.wait()
    return p.returncode, p.stdout.read()


def version():
    # yes yes yes, its an evil hack, if you know of a better way, refactor it!
    git_dir = sys.modules['pypsignifit'].__path__[0].replace('pypsignifit','')
    ret_code, HEAD_SHA, = __get_command_output("git rev-parse HEAD",
            {"GIT_DIR" : ("%s.git" % git_dir)})
    if ret_code == 128:
        print "Unlikely to be running coverage script from a Git repository."
        return None
    else:
        return HEAD_SHA.strip()

def __test__ ( ):
    "If we call the file directly, we perform a test run"
    import numpy as N
    import pylab as p
    import sys

    if len(sys.argv) == 1 or sys.argv[1] == "bootstrap":
        bootstrap = True
    elif sys.argv[1] == "bayes":
        bootstrap = False

    x = [float(2*k) for k in xrange(6)]
    k = [34,32,40,48,50,48]
    n = [50]*6
    d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    d = N.array(zip(x,k,n))
    priors = ("flat","flat","Uniform(0,0.1)")

    if bootstrap:
        b = BootstrapInference ( d, sample=2000, priors=priors )
        GoodnessOfFit(b)
        ParameterPlot(b)
    else:
        priors = ("Gauss(0,4)","Gamma(1,3)","Beta(2,30)")
        mcmc = BayesInference ( d, sigmoid="cauchy", priors=priors )
        mcmc.sample(start=(6,4,.3))
        mcmc.sample(start=(1,1,.1))
        print "Posterior Intervals",mcmc.getCI()
        print "Model Evidence", mcmc.evidence
        print "Rhat (m):",mcmc.Rhat ()
        print "Nsamples:",mcmc.nsamples
        print "DIC:",mcmc.DIC
        print "pD:", mcmc.pD

        GoodnessOfFit(mcmc)
        for prm in xrange(3):
            ConvergenceMCMC ( mcmc, parameter=prm )
        print mcmc

    p.show()

if __name__ == "__main__":
    __test__()
