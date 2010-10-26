#!/usr/bin/env python

from optparse import OptionParser
from pypsignifit import psigobservers
from pypsignifit import psigcorrect
import pypsignifit
from numpy import array,arange,mgrid,clip,zeros,sort, random
import os,sys
import pylab
import operator
import time

__helptext__ = """
Determine coverage of confidence intervals for a given combination of analysis/generating parameters.
"""

# Analyze command line options
parser = OptionParser ( usage=__helptext__ )

# Generating parameters
parser.add_option ( "--gen-observer", dest="observer", default="binomial",
        help="model for the variance of the observers responses. Valid choices are 'binomial', 'betabinomial', or 'learning'. Default: 'binomial'" )
parser.add_option ( "--gen-observer-params", dest="observer_params", default="",
        help="parameters of the variance model for the generating observer. If the observer is 'binomial', no parameter is required, the 'betabinomial' observer requires an additional parameter parameter m given the dispersion of the beta distribution used to model the observers response variance. If the 'learning' observer is selected, a string with the following parameters should be given: 'L,r', where L is the learning range both prm1 and prm2 vary from gen_prm+L to gen_prm-L. r is the learning rate. That is, both parameters vary from one trial to the next as prm(n+1)=prm(n)+(prm_gen-L-prm(n))/r. Thus, low values of r indicate faster learning. Default: '' (binomial observer), '10' (betabinomial observer), '0.7,40' (learning observer)" )
parser.add_option ( "--gen-sigmoid",  dest="gen_sigmoid", default="logistic",
        help="sigmoid used when generating the data. Default: 'logistic'" )
parser.add_option ( "--gen-core",     dest="gen_core",    default="mw0.1",
        help="core type used when generating the data. Default: 'mw0.1'" )
parser.add_option ( "--gen-prm1",     dest="gen_prm1",    default=4,    type="float",
        help="first parameter of the psychometric function used to generate the data (typically called alpha, or m). Default: 4" )
parser.add_option ( "--gen-prm2",     dest="gen_prm2",    default=2,    type="float",
        help="second parameter of the psychometric function used to generate the data (typically called beta, or w). Default: 2" )
parser.add_option ( "--gen-prm3",     dest="gen_prm3",    default=0.02, type="float",
        help="third parameter of the psychometric function used to generate the data (this is the lapse rate). Default: 0.02" )
parser.add_option ( "--gen-prm4",     dest="gen_prm4",    default=0.02, type="float",
        help="fourth parameter of the psychometric function used to generate the data (this is the guessing rate if a yes-no paradigm is employed). In nAFC, this parameter is not used. Default: 0.02" )
parser.add_option ( "--gen-nafc",     dest="gen_nafc",    default=2,    type="int",
        help="number of alternatives in the generating task (1 indicates a yes-no task). Default: 2" )

# Analyzing parameters
parser.add_option ( "--ana-sigmoid",  dest="ana_sigmoid", default="logistic",
        help="sigmoid used when analyzing the data. Default: logistic" )
parser.add_option ( "--ana-core",     dest="ana_core",    default="mw0.1",
        help="core type used when analyzing the data. Default: 'mw0.1'" )
parser.add_option ( "--ana-nafc",     dest="ana_nafc",    default=None,    type="int",
        help="number of alternatives assumed when analyzing the data. Default is to take the same as gen_nafc" )
parser.add_option ( "--nbootstrap",   dest="nbootstrap",  default=2000, type="int",
        help="number of bootstrap repetitions when determining confidence intervals or goodness of fit statistics. Default: 2000" )

# Simulation options
parser.add_option ( "--nsimulations", dest="nsimulations", default=1000, type="int",
        help="number of full simulation runs used to determine coverage information. Default: 1000" )
parser.add_option ( "--blocksize",    dest="blocksize",    default=10,   type="int",
        help="number of trials per block in the simulated experiment. Default: 10" )
parser.add_option ( "--nblocks",      dest="nblocks",      default=5,    type="int",
        help="number of blocks in the simulated experiment. Default: 5" )
parser.add_option ( "--fixed-levels", dest="fixed_levels", default=None,
        type="string", help="list of stimulus levels, if desired, if None,"+\
                "psychometric function will be sampled. Default: None")
parser.add_option ( "--seed", dest="seed", default="fixed",
        type="string", help="seed for simulation, can be 'fixed, 'time', or an"+\
                "integer value.")

parser.add_option ( "-o", "--output", dest="outputfile", default="test.log",
        help="name of the output file in which the data should be stored. By default, no output file is used" )

parser.add_option ( "--datareduce", dest="datareduce", action="store_true",
        help="reduce data based on the estimated nu parameter" )

options,arguments = parser.parse_args()

print "writing output to",options.outputfile

############################################################
#                                                          #
# Create useful variables from command line arguments      #
#                                                          #
############################################################

# Observer fixed parameters
if options.gen_nafc == 1:
    gen_prm = [ options.gen_prm1, options.gen_prm2, options.gen_prm3, options.gen_prm4 ]
elif options.gen_nafc > 1:
    gen_prm = [ options.gen_prm1, options.gen_prm2, options.gen_prm3 ]
else:
    raise IOError, "gen_nafc should be > 0, but is %d" % ( options.gen_nafc, )

if options.ana_nafc is None:
    options.ana_nafc = options.gen_nafc

# Observer variable parameters
gen_kwargs = { "nafc": options.gen_nafc,
        "sigmoid": options.gen_sigmoid,
        "core":    options.gen_core }

# Create the desired Observer
if options.observer == "binomial":

    def create_new_observer ():
        return psigobservers.Observer ( *gen_prm, **gen_kwargs )

elif options.observer == "betabinomial":

    if options.observer_params=="":
        M = 10
    else:
        M = float(options.observer_params)
    def create_new_observer ():
        return psigobservers.BetaBinomialObserver ( *(gen_prm+[M]), **gen_kwargs )

elif options.observer == "learning":

    if options.observer_params == "":
        L,r = .7,40.
    else:
        L,r = [float(x) for x in options.observer_params.split(",")]
    end1 = options.gen_prm1 - L
    end2 = options.gen_prm2 - L
    start_prm = list(gen_prm)
    start_prm[0] += L
    start_prm[1] += L
    def create_new_observer ():
        return psigobservers.LinearSystemLearner ( *(start_prm+[{'a': (r,end1), 'b': (r,end2)}]), **gen_kwargs )

else:
    raise IOError, "Invalid observer model: %s" % ( options.observer, )

# Estimation related parameters
ana_kwargs = { "nafc": options.gen_nafc,
        "sigmoid":     options.ana_sigmoid,
        "core":        options.ana_core }

if not options.ana_nafc is None:
    ana_kwargs["nafc"] = options.ana_nafc

# Create observer
OO = psigobservers.Observer ( *gen_prm, **gen_kwargs )

# Create stimulus levels
if options.fixed_levels is not None:
    # Use the stimulus levels that were given on the command line
    message = "'fixed-levels' must be a sequence of numbers, of length nblocks."
    try:
        x = eval(options.fixed_levels)
    except Exception:
        raise ValueError(message + options.fixed_levels)
    if not operator.isSequenceType(x) or False in [operator.isNumberType(i) for i in x]:
        raise ValueError(message + options.fixed_levels)
    if not len(x) == options.nblocks:
        raise ValueError("Argument mismatch: You gave "+str(len(x))+" levels on the command line"+\
                " but specified "+str(options.nblocks)+" blocks.")
    print "levels:", x
else:
    # Create stimuli (roughly based on results from Hills PhD thesis in chapter 5)
    Fx = mgrid[.1:.99:1j*(options.nblocks)]
    # Fx -= Fx.mean()
    # Fx /= Fx.std()
    # Fx *= 0.22    # this corresponds to Hills PhD thesis (sigma_p ~ 0.11 in 2AFC)
    # if options.ana_nafc == 1:
    #     Fx += 0.5
    # else:
    #     Fx += 0.6   # This corresponds to Hills PhD thesis (p_bar ~ 0.8 in 2AFC)
    # Fx = array( Fx.tolist()+[.99] )
    # Fx = clip(Fx,.001,.999)
    x = OO.getlevels(Fx)
    y = mgrid[0.001:0.999:100j]
    print x,Fx

# Priors
constraints = ["unconstrained","unconstrained","Uniform(0,.1)"]
# priors      = ["Gauss(4,.1)", "Gamma(1,4)","Beta(2,50)"]   # Hilft auch nicht so viel
priors      = ["Gauss(0,100)", "Gamma(1.01,2000)","Beta(2,50)"]
if options.ana_nafc < 2:
    constraints += ["Uniform(0,.1)"]
    priors += ["Beta(1,10)"]
print priors

# Organize output
if len(options.outputfile) > 0:
    if os.path.exists( options.outputfile ):
        sys.stderr.write ( "Output file %s exists." %(options.outputfile,) )
        while True:
            decision = raw_input ( " Overwrite? [Y/N] " )
            if decision.upper() == "Y":
                break
            elif decision.upper() == "N":
                sys.stderr.write ("Terminating\n")
                sys.exit()
    outfile = open(options.outputfile,"w")
else:
    outfile = os.path.devnull

# Parse and set seed

def check_int(value):
    try:
        int(value)
        return True
    except ValueError:
        return False

if options.seed not in ["fixed", "time"] and not check_int(options.seed):
    raise ValueError("'seed' must be either 'fixed', 'time' or an integer value.")
elif options.seed == 'fixed':
    print "Seed is default."
elif options.seed == 'time':
    seed = int(time.time())
    print "Seed is time since epoch in seconds: '%d'" % seed
    pypsignifit.set_seed(seed)
else:
    seed = int(options.seed)
    print "Seeed is value given on command line: '%d'" % seed
    pypsignifit.set_seed(seed)

############################################################
#                                                          #
# Perform the simulation                                   #
#                                                          #
############################################################

def check_ci ( observer, inference ):
    true_thres = observer.thres
    ci = inference.getCI ( 1, (.025, .975) )
    if ci[0]<true_thres and ci[1]>true_thres:
        return 1
    else:
        return 0

def writelog ( f, Bnpr=None, Bpar=None, mcmc=None, mcmc_conv=1 ):
    ptile = pylab.prctile
    if Bnpr is None:
        outs = "run m.gen w.gen "
        outs += "m.npr.e m.npr.l m.npr.h w.npr.e w.npr.l w.npr.h d.npr d.npr.crit nu.npr rpd.npr rpd.npr.l rpd.npr.h rkd.npr rkd.npr.l rkd.npr.h infl.npr." \
                + " infl.npr.".join([str(x) for x in range(options.nblocks)]) + " "
        outs += "m.par.e m.par.l m.par.h w.par.e w.par.l w.par.h d.par d.par.crit nu.par rpd.par rpd.par.l rpd.par.h rkd.par rkd.par.l rkd.par.h infl.par." \
                + " infl.par.".join([str(x) for x in range(options.nblocks)]) + " "
        outs += "m.bay.e m.bay.l m.bay.h w.bay.e w.bay.l w.bay.h d.bay d.bay.p nu.bay rpd.bay rpd.bay.p rkd.bay rkd.bay.p conv.bay Rhat.0 Rhat.1 Rhat.2 infl.bay." \
                + " infl.bay.".join([str(x) for x in range(options.nblocks)])
        outs += "\n"
    else:
        outs = "%d %g %g " % (simulation, gen_prm[0], gen_prm[1] )
        # Bnpr
        outs += "%g %g %g " % (Bnpr.estimate[0],Bnpr.getCI(1,(.025,)),Bnpr.getCI(1,(.975))) # m.npr.e m.npr.l m.npr.h
        outs += "%g %g %g " % (Bnpr.estimate[1],ptile(Bnpr.mcestimates[:,1],2.5),ptile(Bnpr.mcestimates[:,1],97.5)) # w.npr.e w.npr.l w.npr.h
        outs += "%g %g %g " % (Bnpr.deviance, ptile(Bnpr.mcdeviance,95), psigcorrect.estimate_nu (Bnpr)[0]) # d.npr d.npr.crit nu.npr
        outs += "%g %g %g " % (Bnpr.Rpd,ptile(Bnpr.mcRpd,2.5),ptile(Bnpr.mcRpd,97.5)) # rpd.npr rpd.npr.l rpd.npr.h
        outs += "%g %g %g " % (Bnpr.Rkd,ptile(Bnpr.mcRkd,2.5),ptile(Bnpr.mcRkd,97.5)) # rkd.npr rkd.npr.l rkd.npr.h
        outs += ("%g "*options.nblocks) % tuple(Bnpr.infl)
        # Bpar
        outs += "%g %g %g " % (Bpar.estimate[0],Bpar.getCI(1,(.025,)),Bpar.getCI(1,(.975))) # m.par.e m.par.l m.par.h
        outs += "%g %g %g " % (Bpar.estimate[1],ptile(Bpar.mcestimates[:,1],2.5),ptile(Bpar.mcestimates[:,1],97.5)) # w.par.e w.par.l w.par.h
        outs += "%g %g %g " % (Bpar.deviance, ptile(Bpar.mcdeviance,95), psigcorrect.estimate_nu (Bpar)[0]) # d.par d.par.crit nu.par
        outs += "%g %g %g " % (Bpar.Rpd,ptile(Bpar.mcRpd,2.5),ptile(Bpar.mcRpd,97.5)) # rpd.par rpd.par.l rpd.par.h
        outs += "%g %g %g " % (Bpar.Rkd,ptile(Bpar.mcRkd,2.5),ptile(Bpar.mcRkd,97.5)) # rkd.par rkd.par.l rkd.par.h
        outs += ("%g "*options.nblocks) % tuple(Bpar.infl)
        # Bay
        outs += "%g %g %g " % (mcmc.estimate[0],mcmc.getCI(1,(.025,)),mcmc.getCI(1,(.975))) # m.bay.e m.bay.l m.bay.h
        outs += "%g %g %g " % (mcmc.estimate[1],ptile(mcmc.mcestimates[:,1],2.5),ptile(mcmc.mcestimates[:,1],97.5)) # m.bay.e m.bay.l m.bay.h
        outs += "%g %g %g " % (mcmc.deviance,mcmc.bayesian_p('deviance'), psigcorrect.estimate_nu (mcmc)[0]) # d.bay d.bay.p
        outs += "%g %g " % (mcmc.Rpd,mcmc.bayesian_p('Rpd')) # d.rpd d.rpd.p
        outs += "%g %g " % (mcmc.Rkd,mcmc.bayesian_p('Rkd')) # d.rkd d.rkd.p
        outs += "%d %g %g %g " % (mcmc_conv,mcmc.Rhat(0),mcmc.Rhat(1),mcmc.Rhat(2))
        outs += ("%g "*options.nblocks) % tuple(mcmc.infl)
        outs += "\n"
    f.write ( outs )
    return

writelog ( outfile )
if options.datareduce:
    outfile_reduced = open ( options.outputfile + "reduce","w" )
    writelog ( outfile_reduced )

count_npr = 0.
count_par = 0.
count_bay = 0.
not_converged = 0

sys.stderr.write("\n")
for simulation in xrange ( options.nsimulations ):
    sys.stderr.write ( "\nSimulation %d is running" % ( simulation, ) )
    O = create_new_observer ()
    # print "\nO=",O
    random.shuffle ( x )
    data = O.DoAnExperiment ( x, ntrials=options.blocksize )
    print "\ndata =",data
    print constraints
    Bnpr = pypsignifit.BootstrapInference ( data, sample=options.nbootstrap, priors=constraints, parametric=False, **ana_kwargs )
    print "Done npar"
    Bpar = pypsignifit.BootstrapInference ( data, sample=options.nbootstrap, priors=constraints, parametric=True,  **ana_kwargs )
    print "Done par"

    ####################
    # How to make sure that in the end ALL chains have converged?
    # We can give upper and lower limits for m and w from our sampling positions.
    # m cannot be outside the sampled range and w should not be wider than the sampled range (or twice that)
    mcmc = pypsignifit.BayesInference ( data, sample=True, priors=priors, **ana_kwargs )
    for prm in [0,1,2]:
        if not mcmc.geweke(prm)[2] is None:
            for j in mcmc.geweke(prm)[2]:
                mcmc.resample(j)
    N = mcmc.mcestimates.shape[0]
    mcmc.sample( start = mcmc.farstart )
    mcmc.sample( start = mcmc.farstart )
    for prm in [0,1,2]:
        if not mcmc.geweke(prm)[2] is None:
            for j in mcmc.geweke(prm)[2]:
                mcmc.resample(j)
    print "Rhat:  ",mcmc.Rhat(0),mcmc.Rhat(1),mcmc.Rhat(2)
    print "Geweke:",mcmc.geweke(0)[2],mcmc.geweke(1)[2],mcmc.geweke(2)[2]
    mcmc_conv = 1
    if mcmc.Rhat (0)>1.1 or mcmc.Rhat (1)>1.1 or mcmc.Rhat (2)>1.1:
        not_converged += 1
        mcmc_conv = 0
        # pypsignifit.ConvergenceMCMC(mcmc,0)
        # pypsignifit.ConvergenceMCMC(mcmc,1)
        # pypsignifit.ConvergenceMCMC(mcmc,2)
        # pypsignifit.GoodnessOfFit(mcmc)
        # pypsignifit.show()
        # sys.exit()
    # pypsignifit.ConvergenceMCMC(mcmc,1)
    else:
        count_bay += check_ci ( O, mcmc )

    count_npr += check_ci ( O, Bnpr )
    count_par += check_ci ( O, Bpar )
    # print count_bay, mcmc.estimate, pylab.prctile(mcmc.mcestimates[:,0], (2.5,97.5)), pylab.prctile(mcmc.mcestimates[:,1], (2.5,97.5))
    print count_bay, mcmc.getCI(1,(.025,.975))

    writelog ( outfile, Bnpr, Bpar, mcmc, mcmc_conv )

    if options.datareduce:
        data = array ( data )
        dataml = data.copy ()
        nu = psigcorrect.estimate_nu (Bpar)[0]
        print "==============", nu, "==============="
        dataml[:,1] = ( data[:,1] * nu ).astype("i")
        dataml[:,2] = ( data[:,2] * nu ).astype("i")
        print dataml
        try:
            Bnpr = pypsignifit.BootstrapInference ( dataml, sample=options.nbootstrap, priors=constraints, parametric=False, **ana_kwargs )
        except:
            sys.stderr ( "An error ocurred in simulation %d during nonparametric bootstrap\n" % ( simulation,) )
        print "Done npar"
        try:
            Bpar = pypsignifit.BootstrapInference ( dataml, sample=options.nbootstrap, priors=constraints, parametric=True,  **ana_kwargs )
        except:
            sys.stderr ( "An error ocurred in simulation %d during parametric bootstrap\n" % ( simulation,) )
        print "Done par"

        ####################
        # How to make sure that in the end ALL chains have converged?
        # We can give upper and lower limits for m and w from our sampling positions.
        # m cannot be outside the sampled range and w should not be wider than the sampled range (or twice that)
        datamcmc = data.copy ()
        # nu = psigcorrect.estimate_nu (mcmc)[0]
        print "==============", nu, "==============="
        datamcmc[:,1] = ( data[:,1] * nu ).astype("i")
        datamcmc[:,2] = ( data[:,2] * nu ).astype("i")
        print datamcmc
        try:
            mcmc = pypsignifit.BayesInference ( datamcmc, sample=True, priors=priors, **ana_kwargs )
            for prm in [0,1,2]:
                if not mcmc.geweke(prm)[2] is None:
                    for j in mcmc.geweke(prm)[2]:
                        mcmc.resample(j)
            N = mcmc.mcestimates.shape[0]
            mcmc.sample( start = mcmc.farstart )
            mcmc.sample( start = mcmc.farstart )
            for prm in [0,1,2]:
                if not mcmc.geweke(prm)[2] is None:
                    for j in mcmc.geweke(prm)[2]:
                        mcmc.resample(j)
        except:
            sys.stderr ( "An error ocurred in simulation %d during mcmc sampling\n" % ( simulation,) )
        print "Rhat:  ",mcmc.Rhat(0),mcmc.Rhat(1),mcmc.Rhat(2)
        print "Geweke:",mcmc.geweke(0)[2],mcmc.geweke(1)[2],mcmc.geweke(2)[2]
        mcmc_conv = 1
        if mcmc.Rhat (0)>1.1 or mcmc.Rhat (1)>1.1 or mcmc.Rhat (2)>1.1:
            not_converged += 1
            mcmc_conv = 0
        writelog ( outfile_reduced, Bnpr, Bpar, mcmc, mcmc_conv )

sys.stderr.write ( "\r"+50*" "+"\n" )

outfile.close()

print "Coverages:"
print "  nonparametric bootstrap:", count_npr/options.nsimulations
print "  parametric bootstrap:   ", count_par/options.nsimulations
print "  MCMC (bayesian):        ", count_bay/(options.nsimulations-not_converged)
print "  MCMC runs that did not converge:",not_converged

# pypsignifit.show()
