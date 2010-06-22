#!/usr/bin/env python

import os

submission_text = """universe = vanilla
notification = Always
should_transfer_files = NO
initialdir = /home/ingo/Hacking/Psi++/scripts/nonstationarity/
executable = coverage.py
environment = "PYTHONPATH=/home/ingo/lib/python/"

"""

sim_root = "simulation_data"

os.makedirs ( sim_root )

# Binomial
for nblocks in [4,6,8,12,24]:
    for blocksize in [10,20,40]:
        for sigm in ["weibull","logistic","rgumbel"]:
            for nalternatives in [1,2]:
                if sigm=="logistic":
                    sigm_args = ""
                elif sigm=="rgumbel":
                    sigm_args = "--ana-sigmoid=rgumbel"
                elif sigm=="weibull":
                    sigm_args = "--ana-sigmoid=exp --ana-core=poly"
                submission_text += "error = %s/%dafc/%s/binomial_%d_%d.log\n" % (sim_root,nalternatives,sigm,blocksize,nblocks)
                submission_text += "arguments = --blocksize=%d --nblocks=%d " % (blocksize,nblocks)
                submission_text += " --gen-nafc=%d --ana-nafc=%d %s" % (nalternatives,nalternatives,sigm_args)
                submission_text += " --datareduce"
                submission_text += " --output=%s/%dafc/%s/binomial_%d_%d.dat\nqueue\n\n" % (sim_root,nalternatives,sigm,blocksize,nblocks)
                try:
                    os.makedirs ( os.path.join ( sim_root, "%dafc"%(nalternatives,), "%s" % (sigm,) ) )
                except:
                    pass

# Betabinomial
for nblocks in [4,6,8,12,24]:
    for blocksize in [10,20,40]:
        for sigm in ["weibull","logistic","rgumbel"]:
            for nalternatives in [1,2]:
                if sigm=="logistic":
                    sigm_args = ""
                elif sigm=="rgumbel":
                    sigm_args = "--ana-sigmoid=rgumbel"
                elif sigm=="weibull":
                    sigm_args = "--ana-sigmoid=exp --ana-core=poly"
                submission_text += "error = %s/%dafc/%s/betabinomial_%d_%d.log\n" % (sim_root,nalternatives,sigm,blocksize,nblocks)
                submission_text += "arguments = --blocksize=%d --nblocks=%d --gen-observer=betabinomial " % (blocksize,nblocks)
                submission_text += " --gen-nafc=%d --ana-nafc=%d %s" % (nalternatives,nalternatives,sigm_args)
                submission_text += " --datareduce"
                submission_text += " --output=%s/%dafc/%s/betabinomial_%d_%d.dat\nqueue\n\n" % (sim_root,nalternatives,sigm,blocksize,nblocks)
                try:
                    os.makedirs ( os.path.join ( sim_root, "%dafc"%(nalternatives,), "%s" % (sigm,) ) )
                except:
                    pass

# Learning
for nblocks in [4,6,8,12,24]:
    for blocksize in [10,20,40]:
        for sigm in ["weibull","logistic","rgumbel"]:
            for nalternatives in [1,2]:
                if sigm=="logistic":
                    sigm_args = ""
                elif sigm=="rgumbel":
                    sigm_args = "--ana-sigmoid=rgumbel"
                elif sigm=="weibull":
                    sigm_args = "--ana-sigmoid=exp --ana-core=poly"
                submission_text += "error = %s/%dafc/%s/learning_%d_%d.log\n" % (sim_root,nalternatives,sigm,blocksize,nblocks)
                submission_text += "arguments = --blocksize=%d --nblocks=%d --gen-observer=learning " % (blocksize,nblocks)
                submission_text += " --gen-nafc=%d --ana-nafc=%d %s" % (nalternatives,nalternatives,sigm_args)
                submission_text += " --datareduce"
                submission_text += " --output=%s/%dafc/%s/learning_%d_%d.dat\nqueue\n\n" % (sim_root,nalternatives,sigm,blocksize,nblocks)
                try:
                    os.makedirs ( os.path.join ( sim_root, "%dafc"%(nalternatives,), "%s" % (sigm,) ) )
                except:
                    pass

print submission_text
