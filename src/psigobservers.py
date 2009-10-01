#!/usr/bin/env python

import numpy as N
import _psipy

__doc__ = """
When we want to know how well a psychometric function can describe an observers behavior, we
may want to simulate an observer. This module implements a number of simulated observers.
The basic observer does not violate any assumptions. However more elaborated observers
violate some of the assumptions that are typical when fitting psychometric functions.
"""

class Observer ( object ):
    def __init__ ( self, *params, **model ):
        """A stationary binomial observer

        :Parameters:
            *params* :
                a list of parameters in the model. For nAFC tasks the parameters are a,b,lapse
                for a Yes/No task the parameters are a,b,lapse,guess
            *model* :
                a list of keyword arguments specifying the model. These are the same as in
                the psignidata module
        """
        if model.setdefault( "nafc", 2 ) == 1:
            self.a,self.b,self.lapse,self.guess = params
        else:
            self.a,self.b,self.lapse = params
            self.guess = 1./model["nafc"]
        self.model = {
                "sigmoid": model.setdefault ( "sigmoid", "logistic" ),
                "core":    model.setdefault ( "core",    "ab" ),
                "nafc":    model.setdefault ( "nafc",    2 )
                }
