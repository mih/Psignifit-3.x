#!/usr/bin/env python

import unittest as ut
import pypsignifit.psignidata as pd
import swignifit.swignifit_raw as sft

import numpy as np

def approximatedly_equal ( x, y, eps=1e-7 ):
    return abs ( x-y ) > eps

class TestPsiInference ( ut.TestCase ):
    def setUp ( self ):
        self.pinf = pd.PsiInference ()
        self.pinf.estimate = [2.,1.,.02]
    def test_evaluate ( self ):
        evaluated = self.pinf.evaluate ( [1,2],[2.,1.,.02] )
        self.assertEqual ( approximatedly_equal ( evaluated[0], 0.62909188225759771), 0 )
        self.assertEqual ( approximatedly_equal ( evaluated[1], 0.74), 0 )
        evaluated = self.pinf.evaluate ( [1,2] )
        self.assertEqual ( approximatedly_equal ( evaluated[0], 0.62909188225759771), 0 )
        self.assertEqual ( approximatedly_equal ( evaluated[1], 0.74), 0 )
    def test_getThres ( self ):
        self.assertRaises ( NotImplementedError, self.pinf.getThres, .5 )
        self.pinf.data = [[1,2,3]]
        evaluated = self.pinf.getThres ( .5 )
        self.assertEqual ( approximatedly_equal ( evaluated, 2.), 0 )
        evaluated = self.pinf.getThres ( .3 )
        self.assertEqual ( approximatedly_equal ( evaluated, 1.1527021396127963), 0 )
    def test_repr ( self ):
        self.assertEqual ( self.pinf.__repr__(), "< PsiInference object >" )
    def test_properties ( self ):
        self.assertEqual ( self.pinf.desc, 'sigmoid: logistic\ncore: ab\nnAFC: 2' )
        self.assertEqual ( self.pinf.label, "Psychometric function fit" )
        self.assertEqual ( self.pinf.color, "b" )
        self.assertEqual ( self.pinf.linestyle, "-" )
        self.assertEqual ( self.pinf.marker, "o" )
        self.assertEqual ( self.pinf.linewidth, 1 )

class TestBootstrapInference ( ut.TestCase ):
    def setUp ( self ):
        sft.setSeed(0)
        nafc = 2
        stimulus_intensities = [0.0,2.0,4.0,6.0,8.0,10.0]
        number_of_correct = [34,32,40,48,50,48]
        number_of_trials  = [50]*len(stimulus_intensities)
        data = zip(stimulus_intensities,number_of_correct,number_of_trials)
        self.parametric    = pd.BootstrapInference ( data, priors=("","","Beta(2,30)"), parametric=True )
        self.nonparametric = pd.BootstrapInference ( data, priors=("","","Beta(2,30)"), parametric=False )
    def test_map ( self ):
        map1 = self.parametric.estimate
        map2 = self.nonparametric.estimate
        should_be = [ 2.7375488579273224, 1.4039342533896364, 0.020320093764199146 ]
        for val1,val2,val3 in zip(map1,map2,should_be):
            self.assertEqual ( approximatedly_equal ( val1, val2 ), 0 )
            self.assertEqual ( approximatedly_equal ( val1, val3 ), 0 )

    def test_boots ( self ):
        self.parametric.sample ()
        self.nonparametric.sample ()
        parci = self.parametric.getCI(1)
        nprci = self.nonparametric.getCI(1)
        self.assertEqual ( approximatedly_equal ( parci[0], 1.6100394 ), 0 )
        self.assertEqual ( approximatedly_equal ( parci[1], 3.8797201 ), 0 )
        self.assertEqual ( approximatedly_equal ( nprci[0], 1.01129861 ), 0 )
        self.assertEqual ( approximatedly_equal ( nprci[1], 4.0417433 ), 0 )

        self.assertEqual ( self.parametric.nsamples, 2000 )
        self.assertEqual ( self.nonparametric.nsamples, 2000 )

        self.assertEqual ( approximatedly_equal ( self.parametric.deviance,  8.1689126711025022 ), 0 )
        self.assertEqual ( approximatedly_equal ( self.nonparametric.deviance,  8.1689126711025022 ), 0 )

    def test_sensitivity ( self ):
        self.parametric.sensitivity_analysis (verbose=False)
        parci = [ 1.52167887, 3.8797201 ]
        extci = self.parametric.getCI(1)
        for par,ext in zip(parci,extci):
            self.assertEqual ( approximatedly_equal ( par, ext ), 0 )

class TestBayesInference ( ut.TestCase ):
    def setUp ( self ):
        sft.setSeed(0)
        nafc = 2
        stimulus_intensities = [0.0,2.0,4.0,6.0,8.0,10.0]
        number_of_correct = [34,32,40,48,50,48]
        number_of_trials  = [50]*len(stimulus_intensities)
        data = zip(stimulus_intensities,number_of_correct,number_of_trials)
        self.mcmc = pd.BayesInference ( data, priors=("Gauss(0,100)","Gamma(1.01,200)","Beta(2,30)") )
    def test_all ( self ):
        mapest = self.mcmc.mapestimate
        meanest = self.mcmc.estimate
        map_target = [ 2.73973931, 6.15554732, 0.02034599]
        mean_target =[ 2.58517722, 6.85365569, 0.02872328]
        steps_made = self.mcmc._steps
        steps = [ 0.77238055, 2.56652649, 0.01368064]
        burnin = 12
        thinning = 1
        nsamples = 1997

        for k in xrange ( 3 ):
            self.assertEqual ( approximatedly_equal ( mapest[k], map_target[k] ), 0 )
            self.assertEqual ( approximatedly_equal ( meanest[k], mean_target[k] ), 0 )
            self.assertEqual ( approximatedly_equal ( steps_made[k], steps[k] ), 0 )

        self.assertEqual ( approximatedly_equal ( self.mcmc.bayesian_p(),0.119679519279 ), 0 )
        self.assertEqual ( burnin, self.mcmc.burnin )
        self.assertEqual ( thinning, self.mcmc.thin )
        self.assertEqual ( nsamples, self.mcmc.nsamples )


if __name__ == "__main__":
    ut.main()
