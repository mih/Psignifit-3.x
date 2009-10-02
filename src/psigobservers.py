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

__all__ = ["Observer","LinearSystemLearner"]

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

        # Make sure, a,b,lapse and guess are floats
        self.a = float(self.a)
        self.b = float(self.b)
        self.lapse = float(self.lapse)
        self.guess = float(self.guess)

        self.model = {
                "sigmoid": model.setdefault ( "sigmoid", "logistic" ),
                "core":    model.setdefault ( "core",    "ab" ),
                "nafc":    model.setdefault ( "nafc",    2 )
                }
        self.data = []

    def DoATrial ( self, stimulus_intensity=1 ):
        """Simulate a single trial

        :Parameters:
            *stimulus_intensity* :
                stimulus intensity to be presented

        :Output:
            The response in the trail (1/0 coding for Yes/No in Yes/No-Tasks or for
            Correct/Incorrect in nAFC tasks)
        """
        prob = float( _psipy.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

        resp = int ( N.random.rand()<prob )

        if len(self.data) == 0 or not stimulus_intensity==self.data[-1][0]:
            self.data.append ( [stimulus_intensity,resp,1] )
        else:
            self.data[-1][1] += resp
            self.data[-1][2] += 1

        return resp


    def DoABlock ( self, stimulus_intensity=1, ntrials=50 ):
        """Simulate a block of trials

        :Parameters:
            *stimulus_intensity* :
                stimulus intensity to be presented
            *ntrials* :
                number of trials in the block

        :Output:
            The number of Yes-responses (in Yes/No) or the number of correct responses (nAFC)
        """
        prob = float( _psipy.diagnostics ( [stimulus_intensity], self.params,
            sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

        resp = N.random.binomial ( ntrials, prob )

        self.data.append ( [stimulus_intensity, resp, ntrials] )

        return resp

    def DoAnExperiment ( self, stimulus_intensities, ntrials=50 ):
        """Simulate a whole experiment

        :Parameters:
            *stimulus_intensities* :
                a sequence of stimulus intensities in the order they should be presented
            *ntrials* :
                Either an integer or a sequence. If this is an integer, it is interpreted as a constant
                number of trials that is presented in each trial, otherwise, the sequence is expected to
                contain a number of trials for each trials.

        :Output:
            A list of lists. Each element of the list has the structure
            [stimulus_intensity, number_of_correct_or_yes_responses,number_of_trials]
        """
        if isinstance ( ntrials, int ):
            ntrials = [ntrials] * len(stimulus_intensities)

        data = []

        for s,n in zip ( stimulus_intensities, ntrials ):
            data.append ( [s, self.DoABlock ( s, n ), n] )

        return data

    def seed ( self, seed ):
        """Seed the underlying random number generator to a defined value"""
        N.random.seed ( seed )

    @Property
    def params ():
        "parameters of the model"
        def fget ( self ):
            if self.model["nafc"] < 2:
                return [self.a,self.b,self.lapse,self.guess]
            else:
                return [self.a,self.b,self.lapse]

class LinearSystemLearner ( Observer ):
    def __init__( self, *params, **model ):
        """A nonstationary observer that learns like a linear system in one or more parameters

        For this observer, the parameters of the psychometric function change: The parameters given
        to the constructor are initial values. In addition, a dictionary is given containing the learning
        parameters.

        The constructor requires a dictionary that describes the changes of the parameters. The dictionary
        has keys 'a', 'b', and/or 'lapse'. For each key, there are two parameters:
        the rate of decay and the asymptote where the parameter converges. Thus
        {'a' : ( 40, 3)}
        means that after every trial, the a parameter is modified according to
        a -= (a-3)/4

        :Parameters:
            *params* :
                parameters of the model a, b, lapse and in Yes/No tasks also guess. In addition, a dictionary
                is required to describe the changes of the parameters (see above)
            *model* :
                a list of keyword arguments to describe the model. These are the same as in the psignidata
                module.

        :Example:
        >>> O = LinearSystemLearner ( 7, 2, .02, {'a': (40,3)} )
        >>> O.seed ( 0 )
        >>> O.a
        7.0
        >>> O.DoATrial ( 3 )
        1
        >>> O.DoABlock ( 4, 50 )
        34
        >>> O.DoAnExperiment ( [4,2,8,10,6], 50 )
        [[4, 43, 50], [2, 32, 50], [8, 48, 50], [10, 49, 50], [6, 38, 50]]
        >>> O.data
        [[3, 1, 1], [4, 34, 50], [4, 43, 50], [2, 32, 50], [8, 48, 50], [10, 49, 50], [6, 38, 50]]
        >>> O.a
        3.0019608723226945
        """
        # TODO: An other example for sharper discriminations (b) and for fatigue (lapse)
        Observer.__init__ ( self, *(params[:-1]), **model )
        self.learn = params[-1]

    def DoATrial ( self, stimulus_intensity=1 ):
        """Simulate a single trial with learning

        :Parameters:
            *stimulus_intensity* :
                intensity of the presented stimulus

        :Output:
            either 1 or 0 indicating Yes/No in Yes/No-Tasks or Correct/Incorrect in nAFC
        """
        # Get the response:
        resp = Observer.DoATrial ( self, stimulus_intensity )

        # Modify parameters:
        for prm,vals in self.learn.iteritems():
            if prm=="a":
                self.a -= (self.a-vals[1])/vals[0]
            elif prm=="b":
                self.b -= (self.b-vals[1])/vals[0]
            elif prm=="lapse":
                self.lapse -= (self.self.lapse-vals[1])/vals[0]
            else:
                # This should issue a warning
                pass

        return resp

    def DoABlock ( self, stimulus_intensity=1, ntrials=50 ):
        """Simulate a block of trials with learning

        :Parameters:
            *stimulus_intenstiy* :
                intensity of the presented stimulus
            *ntrials* :
                number of repetitions of the stimulus

        :Output:
            the number of Yes-responses in a Yes/No-task or the number of correct responses in nAFC
        """
        self.data.append ( [stimulus_intensity, 0, 0] )
        resp = 0
        for k in xrange(ntrials):
            resp += self.DoATrial ( stimulus_intensity )
        return resp

if __name__ == "__main__":
    O = Observer ( 4,.8,.02 )
    print O.DoABlock ( 3, 50 )
