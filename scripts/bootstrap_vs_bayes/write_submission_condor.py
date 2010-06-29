#/usr/bin/env python

""" script to generate condor submission file to run simulations for Valentin's
master thesis """

import os

output = """universe = vanilla
notification = Always
notify_user = val
should_transfer_files = NO
initialdir = /home/val/git-working/psignifit/scripts/bootstrap_vs_bayes/
executable =/home/val/git-working/psignifit/scripts/nonstationarity/coverage.py
environment ="PYTHONPATH=/home/val/git-working/psignifit/build/lib.linux-x86_64-2.5"

"""

sim_root = "/home/val/thesis/simulation_dir"
data_root = os.path.join(sim_root, "simulation_data")
log_root =  os.path.join(sim_root, "simulation_logs")

if os.path.exists(sim_root):
    raise Exception("simulation data directory exists, DON'T OVERWRITE IT!!!")
else:
    os.makedirs (data_root)
    os.makedirs (log_root)

# Simulations

job_number = 0

for (nblocks, blocksize, lapse_rate, width, nalt, gen_sigmoid, ana_sigmoid) in \
    ((nblocks, blocksize, lapse_rate, width, nalt, gen_sigmoid, ana_sigmoid)
            for nblocks     in map(str, [4, 8, 16, 32])
            for blocksize   in map(str, [10, 20, 40])
            for lapse_rate  in map(str, [0.01, 0.05])
            for width       in map(str, [0.1, 0.2, 0.3])
            for nalt        in map(str, [1, 2])
            for gen_sigmoid in ["logistic", "rgumbel"]
            for ana_sigmoid in ["logistic", "rgumbel"]):
    output += "queue\n"
    output += "#job id = %d\n" % (job_number)
    output += "error = %s/job%d.log\n" % (log_root, job_number)
    output += "arguments = --nblocks=" + (nblocks)+\
                         " --blocksize=" + (blocksize)+\
                         " --gen-prm3=" + lapse_rate +\
                         " --gen-prm2=" + width +\
                         " --gen-nafc=" + nalt +\
                         " --ana-nafc=" + nalt +\
                         " --gen-sigmoid=" + gen_sigmoid +\
                         " --ana-sigmoid=" + ana_sigmoid +\
                         "\n\n"
    job_number += 1

output += "#total number of jobs = %d" % (job_number)
print output

