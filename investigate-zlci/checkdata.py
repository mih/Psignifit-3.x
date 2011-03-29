#!/usr/bin/env python

import reader
import sys
import pypsignifit as pp

d,s = reader.read_data ( sys.argv[1] )

model = {'nafc':1, 'sigmoid':"logistic", 'core':'mw0.1'}
m = 4.0
w = 4.0
l = 0.05
g = 0.02

priors = ["Gauss(%f,%f)" % (m, m), "Gauss(%f,%f)" % (w, w*2), "Beta(2,50)",
"Beta(1,50)"]
mcmc = pp.BayesInference(d, sigmoid=model['sigmoid'], core=model['core'],
        nafc=model['nafc'], priors=priors, verbose=True)

print mcmc.approx
