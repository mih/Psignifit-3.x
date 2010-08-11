#/usr/bin/env python

""" script to generate condor submission file to run simulations for Valentin's
master thesis """

import os
import time
import pypsignifit.psigobservers as po
from numpy import mgrid


output = """universe = vanilla
notification = Always
notify_user = val
should_transfer_files = NO
initialdir = /home/val/git-working/psignifit/scripts/bootstrap_vs_bayes/
executable =/home/val/git-working/psignifit/scripts/nonstationarity/coverage.py
environment ="PYTHONPATH=/home/val/git-working/psignifit/"

"""

sim_root = "/home/val/thesis/simulation_dir"
data_root = os.path.join(sim_root, "simulation_data")
log_root =  os.path.join(sim_root, "simulation_logs")

data_suffix = ".data"
log_suffix = ".log"
error_suffix = ".error"


if os.path.exists(sim_root):
    raise Exception("simulation data directory exists, DON'T OVERWRITE IT!!!")
else:
    os.makedirs (data_root)
    os.makedirs (log_root)

# create levels

# Fx = mgrid[.1:.99:1j*(nblocks)]
# model ={"sigmoid" : "logistic", "core" : "mw0.1", "nafc" : 1}
# params = [4, 4.0, 0.01, 0.02]

#levels = [[ 3.        ,  3.80913488,  4.3712636 ,  6.09132917],
#          [ 2.        ,  2.88539834,  3.45362578,  3.93235111,  4.4016988 ,
#            4.931913  ,  5.67413566,  8.18265834],
#          [ 2.        ,  2.48609305,  2.8408528 ,  3.13126402,  3.38540599,
#            3.61826975,  3.83938243,  4.05584552,  4.27390645,  4.50009965,
#            4.7425272 ,  5.01299204,  5.33173313,  5.74110516,  6.36386887,
#            8.18265834],
#          [ 2.        ,  2.259244  ,  2.47302265,  2.65706792,  2.8203285 ,
#            2.96840501,  3.10505633,  3.23294705,  3.35405383,  3.46990263,
#            3.58171565,  3.69050733,  3.79715025,  3.90242298,  4.00704705,
#            4.11171774,  4.21713221,  4.32401778,  4.43316336,  4.54545756,
#            4.66193847,  4.78386279,  4.91280722,  5.05082464,  5.20069835,
#            5.36638121,  5.55381231,  5.77258678,  6.03984177,  6.39121983,
#            6.92385846,  8.18265834]]

params = [4.0, 2.0, 0.05, 0.01]
model ={"sigmoid" : "logistic", "core" : "mw0.1", "nafc" : 1}
ob = po.Observer(*params, **model)

def make_levels(nblocks):
    """ will make levels for the desired number of blocks """
    Fx = mgrid[.1:.99:1j*(nblocks)]
    x = ob.getlevels(Fx)
    return x

# Simulations

def str_no_space(l):
    """ concatenate list of numbers into string without spaces """
    return "["+",".join([str(i) for i in l])+"]"

job_number = 0

for (nblocks, blocksize, lapse_rate, width, nalt, gen_sigmoid, ana_sigmoid) in \
    ((nblocks, blocksize, lapse_rate, width, nalt, gen_sigmoid, ana_sigmoid)
            for nblocks     in map(str, [4, 8, 16, 32, 64, 128])
            for blocksize   in map(str, [10, 20, 40])
            for lapse_rate  in map(str, [0.01, 0.05])
            for width       in map(str, [0.5, 2.0, 4])
            for nalt        in map(str, [1, 2])
            for gen_sigmoid in ["logistic", "gumbel_r"]
            for ana_sigmoid in ["logistic", "gumbel_r"]
            if int(nblocks)*int(blocksize) <= 1280):
    output += "#job id = %d\n" % (job_number)
    output += "output = " + os.path.join(log_root, "job"+ str(job_number) + log_suffix) + "\n"
    output += "error = "  + os.path.join(log_root, "job"+ str(job_number) + error_suffix) + "\n"
    output += "arguments = --nblocks=" + nblocks+\
                         " --fixed-levels=" + str_no_space(make_levels(int(nblocks)))+\
                         " --blocksize=" + (blocksize)+\
                         " --gen-prm3=" + lapse_rate +\
                         " --gen-prm2=" + width +\
                         " --gen-nafc=" + nalt +\
                         " --ana-nafc=" + nalt +\
                         " --gen-sigmoid=" + gen_sigmoid +\
                         " --ana-sigmoid=" + ana_sigmoid +\
                         " --seed=" + str(job_number) +\
                         " --output=" + os.path.join(data_root,
                                 "job"+str(job_number)+data_suffix) +\
                         "\n"
    output += "queue\n\n"
    job_number += 1

output += "#total number of jobs = %d" % (job_number)
print output

