#!/usr/bin/env python

import numpy as N
import _psipy

from psignidata import Property

__doc__ = """
When we want to know how well a psychometric function can describe an observers behavior, we
may want to simulate an observer. This module implements a number of simulated observers.
The basic observer does not violate any assumptions. However more elaborated observers
violate some of the assumptions that are typical when fitting psychometric functions.
"""

class Observer ( object ):
    def __init__ ( self, *params, **model ):
        """A stationary binomial observer

        This is the observer, we all want: No interdependencies between trials, no learning,
        no fluctuations in attention or motivation. Perfectly binomial responses in accordance
        with the psychometric function shape you supply.

        :Parameters:
            *params* :
                a list of parameters in the model. For nAFC tasks the parameters are a,b,lapse
                for a Yes/No task the parameters are a,b,lapse,guess
            *model* :
                a list of keyword arguments specifying the model. These are the same as in
                the psignidata module

        :Example:
        >>> O = Observer ( 4, 1, .02 )
        >>> O.seed ( 0 )
        >>> O.DoATrial ( 3 )
        1
        >>> O.DoABlock ( 4, 40 )
        28
        >>> O.DoABlock ( 6, 40 )
        37
        >>> O.DoAnExperiment ( [2,4,6], 50 )
        [[2, 27, 50], [4, 38, 50], [6, 46, 50]]
        >>> O.data
        [[3, 1, 1], [4, 28, 40], [6, 37, 40], [2, 27, 50], [4, 38, 50], [6, 46, 50]]
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

    def DoATrial ( self, stimulus_intensity=1 ):
        """Simulate a single trial

        :Parameters:
            *stimulus_intensity* :
                stimulus intensity to be presented
        """
        prob = float( _psipy.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

        return int ( N.random.rand()<prob )

    def DoABlock ( self, stimulus_intensity=1, ntrials=50 ):
        """Simulate a block of trials

        :Parameters:
            *stimulus_intensity* :
                stimulus intensity to be presented
            *ntrials* :
                number of trials in the block
        """
        prob = float( _psipy.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )
        return N.random.binomial ( ntrials, prob )


    @Property
    def params ():
        "parameters of the model"
        def fget ( self ):
            if self.model["nafc"] < 2:
                return [self.a,self.b,self.lapse,self.guess]
            else:
                return [self.a,self.b,self.lapse]


if __name__ == "__main__":
    O = Observer ( 4,.8,.02 )
    print O.DoABlock ( 3, 50 )
