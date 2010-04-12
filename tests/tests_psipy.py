#!/usr/bin/env python

###################################################
# Check that we can call all functions from _psipy

import _psipy
import unittest as ut

#################### 2afc #########################

def makedata ( neg=False ):
    if neg:
        return zip([1,2,3,4,5,6],[49,50,45,29,26,25],[50]*6)
    else:
        return zip([1,2,3,4,5,6],[25,26,29,45,50,49],[50]*6)

class TestPsipy_2afc ( ut.TestCase ):
    def test_bootstrap ( self ):
        """Call bootstrap"""
        priors = ("","","Uniform(0,.1)")
        _psipy.bootstrap ( makedata(), priors=priors )
    def test_mcmc ( self ):
        """Call mcmc"""
        priors = ("Gauss(0,100)","Gamma(1,4)","Beta(2,50)")
        _psipy.mcmc ( makedata(), start=(4,1,.02), priors=priors )
    def test_map ( self ):
        """Call mapestimator"""
        priors = ("","","Uniform(0,.1)")
        _psipy.mapestimate ( makedata(), priors=priors )
    def test_diagnostics ( self ):
        """Call diagnostics"""
        _psipy.diagnostics ( makedata(), (4,1,.02) )
    def test_sigmoids ( self ):
        """Try different sigmoids (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)")
        for s in ["logistic","cauchy","gauss","rgumbel","lgumbel","exp"]:
            _psipy.mapestimate ( makedata(), priors=priors, sigmoid=s )
    def test_cores ( self ):
        """Try different cores (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)")
        for c in ["ab","mw0.1","poly","weibull","log","linear"]:
            _psipy.mapestimate ( makedata(), priors=priors, core=c )
    def test_weibull ( self ):
        """Try different parameterizations of the weibull (assumes working mapestimator)"""
        _psipy.mapestimate ( makedata(), priors = ("","","Uniform(0,.1)"), sigmoid="exp", core="poly" )
        _psipy.mapestimate ( makedata(), priors = ("","","Uniform(0,.1)"), sigmoid="rgumbel", core="weibull" )
        _psipy.mapestimate ( makedata(), priors = ("","","Uniform(0,.1)"), sigmoid="rgumbel", core="log" )
    def test_priors ( self ):
        """Try different combinations of priors (assumes working mapestimator)"""
        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","Gamma(1,4)","Uniform(0,5)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    _psipy.mapestimate ( makedata(), priors = (mprior,wprior,lprior), core="mw0.1" )

        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","nGamma(1,4)","Uniform(-5,0)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    _psipy.mapestimate ( makedata(True), priors = (mprior,wprior,lprior), core="mw0.1" )

##################### yes/no ########################

def makedata_yn ( neg=False ):
    if neg:
        return zip ( [1,2,3,4,5,6],[40,46,27,25,14,1],[50]*6)
    else:
        return zip ( [1,2,3,4,5,6],[1,14,25,27,46,49],[50]*6)

class TestPsipy_yn ( ut.TestCase ):
    def test_bootstrap ( self ):
        """Call bootstrap"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        _psipy.bootstrap ( makedata_yn(), priors=priors, nafc=1 )
    def test_mcmc ( self ):
        """Call mcmc"""
        priors = ("Gauss(0,100)","Gamma(1,4)","Beta(2,50)","Beta(2,50)")
        _psipy.mcmc ( makedata_yn(), start=(4,1,.02,.02), priors=priors, nafc=1 )
    def test_map ( self ):
        """Call mapestimator"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        _psipy.mapestimate ( makedata_yn(), priors=priors, nafc=1 )
    def test_diagnostics ( self ):
        """Call diagnostics"""
        _psipy.diagnostics ( makedata_yn(), (4,1,.02,.02), nafc=1 )
    def test_sigmoids ( self ):
        """Try different sigmoids (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        for s in ["logistic","cauchy","gauss","rgumbel","lgumbel","exp"]:
            _psipy.mapestimate ( makedata_yn(), priors=priors, sigmoid=s, nafc=1 )
    def test_cores ( self ):
        """Try different cores (assumes working mapestimator)"""
        priors = ("","","Uniform(0,.1)","Uniform(0,.1)")
        for c in ["ab","mw0.1","poly","weibull","log","linear"]:
            _psipy.mapestimate ( makedata_yn(), priors=priors, core=c, nafc=1 )
    def test_weibull ( self ):
        """Try different parameterizations of the weibull (assumes working mapestimator)"""
        _psipy.mapestimate ( makedata_yn(), priors = ("","","Uniform(0,.1)","Uniform(0,.1)"), nafc=1, sigmoid="exp", core="poly" )
        _psipy.mapestimate ( makedata_yn(), priors = ("","","Uniform(0,.1)","Uniform(0,.1)"), nafc=1, sigmoid="rgumbel", core="weibull" )
        _psipy.mapestimate ( makedata_yn(), priors = ("","","Uniform(0,.1)","Uniform(0,.1)"), nafc=1, sigmoid="rgumbel", core="log" )
    def test_priors ( self ):
        """Try different combinations of priors (assumes working mapestimator)"""
        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","Gamma(1,4)","Uniform(0,5)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    for gprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                        _psipy.mapestimate ( makedata_yn(), priors = (mprior,wprior,lprior,gprior), core="mw0.1", nafc=1 )

        for mprior in ["Gauss(0,100)","Uniform(-30,30)",""]:
            for wprior in ["Gauss(0,100)","nGamma(1,4)","Uniform(-5,0)",""]:
                for lprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                    for gprior in ["Uniform(0,.1)","Beta(2,30)","Gauss(0.05,0.01)",""]:
                        _psipy.mapestimate ( makedata_yn(), priors = (mprior,wprior,lprior,gprior), core="mw0.1", nafc=1 )

if __name__ == '__main__':
    ut.main()
