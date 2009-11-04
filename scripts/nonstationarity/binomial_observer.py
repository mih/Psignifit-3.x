#!/usr/bin/env python

from multiprocessing import Process
import analyzerun
import pypsignifit.psigobservers as observers
import numpy as N
import sys,os

O = observers.Observer ( 4, 2, .02 )
pos = N.array([-2,-1,0.,1,2,3],'d')
constraints = ("","","Uniform(0,.1)")

k = 0
sys.stderr.write("\n")
for stimpat,stimuli in enumerate([ 4+pos, 4+ 0.5*N.sign(pos)*pos**2, 4+0.2*pos**3]):
    for n in [20,40,60,80,100]:
        for nblocksperposition in [1,2]:
            k+=1
            sys.stderr.write ( "\rSubmitting job %d" % (k,) )
            sys.stderr.flush()
            blocks = stimuli.tolist()*2
            fname = os.path.join ( "data","binomial","logistic_ab_k%d_n%d_m%d.dat" % (stimpat,n,nblocksperposition) )
            Process ( target=analyzerun.performbootstrap, name="pypsignifit_%d" % (k,), args=(fname,), kwargs={"stimuli": blocks, "ntrials": n, "observer": O } ).start()

sys.stderr.write ( "\rSubmitted all jobs\n" )
sys.stderr.flush()


sys.stderr.write ( "Done\n" )
