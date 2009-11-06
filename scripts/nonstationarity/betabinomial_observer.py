#!/usr/bin/env python

from multiprocessing import Process
import analyzerun
import pypsignifit.psigobservers as observers
import numpy as N
import sys,os


O = observers.BetaBinomialObserver ( 4, 2, .02, 10 )
pos = N.array([-2,-1,0.,1,2,3],'d')
constraints = ("","","Uniform(0,.1)")

k = 0
sys.stderr.write("\n")
for stimpat,stimuli in enumerate( [
        [.3,.4,.48,.52,.6,.7],
        [.1,.3,.4,.6,.7,.9],
        [.3,.44,.7,.8,.9,.98],
        [.1,.2,.3,.4,.5,.6],
        [.08,.18,.28,.7,.85,.99],
        [.3,.4,.5,.6,.7,.99],
        [.34,.44,.54,.8,.9,.98]
        ] ):
    for n in [20,40,60,80,100]:
        for nblocksperposition in [1,2]:
            k+=1
            sys.stderr.write ( "\rSubmitting job %d" % (k,) )
            sys.stderr.flush()
            blocks = O.getlevels ( stimuli )
            fname = os.path.join ( "data","betabinomial","logistic_ab_k%d_n%d_m%d.dat" % (stimpat,n,nblocksperposition) )
            Process ( target=analyzerun.performbootstrap, args=(fname,), kwargs={"stimuli": blocks, "ntrials": n/nblocksperposition, "observer": O } ).start()

sys.stderr.write( " Submitted all jobs\n" )
sys.stderr.flush()

sys.stderr.write ( "Done\n" )
